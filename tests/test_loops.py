import random
from fractions import Fraction

import pytest
from support import (
    HALF,
    clock,
    coin,
    first_certain_return,
    loop,
    mass,
    random_events,
    random_universe,
    swap,
)

from hyperprobability import (
    INFINITY,
    OMEGA,
    ONE,
    ZERO,
    OmegaFraction,
    Ordinal,
    Universe,
    frequency,
    hyperprobability,
    kac,
    refine,
    self_simulation,
    strange_loops,
)
from hyperprobability._exact import solve

SEEDS = range(60)


def mean_returns(universe: Universe, event: frozenset, states: tuple) -> dict:
    """The mean return time to ``event`` from each of its states in the closed class ``states``."""
    outside = [q for q in states if q not in event]
    hitting = {}
    if outside:
        # The mean hitting time of the event solves h = 1 + K h off the event.
        matrix = [[Fraction(q == r) - universe.step(q).get(r, 0) for r in outside] for q in outside]
        solution = solve(matrix, [[Fraction(1)] for _ in outside])
        hitting = {q: row[0] for q, row in zip(outside, solution)}
    return {
        q: 1 + sum(p * hitting.get(r, 0) for r, p in universe.step(q).items())
        for q in states
        if q in event
    }


def subdivide(coarse: Universe, cost: int) -> Universe:
    """``coarse`` with every step spread over ``cost`` steps.

    State ``(s, i)`` is the ``i``-th evaluation of a step from ``s``. Each limit
    rule is carried over by projecting an attractor onto its phase-0 states and
    sending it to its coarse exit law at phase 0.
    """
    kernel = {}
    for s in coarse.states:
        for i in range(cost - 1):
            kernel[(s, i)] = {(s, i + 1): 1}
        kernel[(s, cost - 1)] = {(t, 0): p for t, p in coarse.step(s).items()}
    rules = []
    for k in range(coarse.depth):
        exits = {a.members: a.exit for a in coarse.attractors(k)}

        def rule(members, exits=exits):
            projected = frozenset(s for s, i in members if i == 0)
            return {(t, 0): p for t, p in exits[projected].items()}

        rules.append(rule)
    return Universe(kernel, rules)


def test_strange_loops():
    u = swap()
    assert [(a.level, a.states) for a in strange_loops(u)] == [
        (0, ("a",)),
        (0, ("b",)),
        (1, ("a", "b")),
    ]
    assert [a.states for a in strange_loops(u, 1)] == [("a", "b")]
    assert strange_loops(u, 2) == ()
    assert strange_loops(clock()) == () and strange_loops(coin()) == ()
    assert [a.states for a in strange_loops(loop())] == [("a", "b")]
    with pytest.raises(ValueError):
        strange_loops(u, -1)


def test_kac_examples():
    (bound,) = kac(clock(), "tick")
    assert bound.attractor.states == ("tick", "tock") and bound.event == frozenset({"tick"})
    assert bound.frequency == HALF and bound.guarantee == 2 and bound.ordinal_probability == HALF
    assert bound.holds and bound.tight

    (bound,) = kac(coin(), "H")
    assert bound.frequency == HALF and bound.guarantee == OMEGA
    assert bound.ordinal_probability == OmegaFraction(1, OMEGA)
    assert bound.holds and not bound.tight

    (bound,) = kac(loop(), "out")
    assert bound.attractor.states == ("out",) and bound.frequency == 1 and bound.guarantee == 1
    assert bound.tight
    (bound,) = kac(loop(), "a")
    assert bound.attractor.states == ("a", "b") and bound.frequency == HALF and bound.tight
    assert kac(loop(), []) == ()


@pytest.mark.parametrize("seed", SEEDS)
def test_the_kac_bridge(seed):
    """SPEC Theorem 6.2, against Kac's lemma and exhaustive search."""
    u = random_universe(seed, size=random.Random(seed).randint(1, 4))
    for event in random_events(u, seed):
        bounds = kac(u, event)
        assert [b.attractor for b in bounds] == [a for a in u.attractors() if a.members & event]
        stopped = u.stop_at(event)
        for bound in bounds:
            states = bound.attractor.states
            assert bound.event == frozenset(q for q in states if q in event)
            assert bound.frequency == sum(bound.attractor.law[q] for q in bound.event)
            assert all(frequency(u, event, q) == bound.frequency for q in states)
            returns = mean_returns(u, event, states)
            # Kac's lemma: from the stationary law on the event, the mean return time is 1/π_A(E).
            assert sum(bound.attractor.law[q] * m for q, m in returns.items()) == 1
            assert bound.guarantee == max(first_certain_return(u, event, q) for q in bound.event)
            assert bound.holds
            if bound.guarantee.is_finite:
                n = bound.guarantee.to_int()
                assert max(returns.values()) <= n and bound.frequency >= Fraction(1, n)
                # Tight exactly when no return comes before step n.
                early = n > 1 and any(
                    mass(stopped.law(u.step(q), n - 2), event) for q in bound.event
                )
                assert bound.tight == (not early)
            else:
                assert bound.guarantee == OMEGA and not bound.tight
                assert bound.ordinal_probability.is_infinitesimal


def test_refine():
    assert refine(OMEGA, 2) == OMEGA
    assert refine(2, OMEGA) == OMEGA * 2
    assert refine(3, 4) == 12
    assert refine(OMEGA + 3, 2) == OMEGA + 6
    assert refine(ZERO, 3) == 0
    assert refine(INFINITY, 5) is INFINITY
    with pytest.raises(ValueError):
        refine(OMEGA, 0)
    for guarantee, cost in ((1.5, 2), (OMEGA, 1.5), ("ω", 2), (True, 2)):
        with pytest.raises(TypeError):
            refine(guarantee, cost)


@pytest.mark.parametrize("seed", range(40))
def test_refine_is_the_guarantee_of_a_subdivided_run(seed):
    """SPEC Proposition 8.1: spreading each step over ``c`` evaluations turns ``ℍ`` into
    ``c·ℍ``."""
    coarse = random_universe(seed)
    for cost in (1, 2, 3):
        fine = subdivide(coarse, cost)
        for event in random_events(coarse, seed):
            marked = [(s, 0) for s in event]
            for s in coarse.states:
                expected = refine(hyperprobability(coarse, event, s), cost)
                assert hyperprobability(fine, marked, (s, 0)) == expected


def test_self_simulation():
    assert self_simulation(lambda a: 1 + a) == OMEGA
    assert self_simulation(lambda a: OMEGA + a) == Ordinal.omega_power(2)
    assert self_simulation(lambda a: OMEGA * 2 + a) == Ordinal.omega_power(2)
    assert self_simulation(lambda a: a) == 0
    assert self_simulation(lambda a: a, start=5) == 5
    assert self_simulation(lambda a: refine(a, 2), start=1) == OMEGA
    assert self_simulation(lambda a: max(a, OMEGA * 3)) == OMEGA * 3
    with pytest.raises(ValueError):
        self_simulation(lambda a: OMEGA if a == 0 else ONE)
