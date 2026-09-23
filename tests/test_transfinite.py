import random
from fractions import Fraction

import pytest
from support import (
    chain,
    coin,
    loop,
    mass,
    ordinals,
    random_events,
    random_universe,
    swap,
)

from hyperprobability import (
    INFINITY,
    OMEGA,
    Ordinal,
    TransfiniteLaw,
    Universe,
    hyperprobability,
)

SEEDS = range(60)


def split_limit(stage: Ordinal) -> tuple[int, Ordinal]:
    k = next(i for i, c in enumerate(stage.coefficients) if c)
    return k, Ordinal(tuple(c - (i == k) for i, c in enumerate(stage.coefficients)))


def deterministic_universe(seed: int) -> Universe:
    """A deterministic universe whose every limit rule names a state."""
    rng = random.Random(seed)
    states = list(range(rng.randint(1, 5)))
    step = {state: rng.choice(states) for state in states}
    rules = []
    for k in range(rng.randint(0, 2)):

        def rule(members, key=f"{seed}/{k}"):
            return random.Random(f"{key}/{sorted(members)}").choice(states)

        rules.append(rule)
    return Universe.deterministic(step, rules)


def test_the_coin():
    law = TransfiniteLaw(coin(), "H", "T")
    assert law.exactly(0) == 0
    for n in range(1, 12):
        assert law.exactly(n) == Fraction(1, 2**n)
    assert law.before(OMEGA) == 1 and law.exactly(OMEGA) == 0 and law.never() == 0
    assert law.guarantee() == OMEGA
    assert law.event == frozenset({"H"}) and law.initial == {"T": 1}
    assert law.universe.states == ("H", "T")
    assert repr(law) == "TransfiniteLaw(event={'H'}, initial={'T': Fraction(1, 1)})"


def test_a_strange_loop_makes_the_law_jump():
    law = TransfiniteLaw(loop(), "out", "a")
    assert law.before(OMEGA) == 0
    assert law.exactly(OMEGA) == 1
    ((attractor, weight),) = law.loops_into(OMEGA).items()
    assert attractor.states == ("a", "b") and attractor.level == 0 and weight == 1
    assert law.loops_into(5) == {} and law.loops_into(0) == {}


def test_the_exit_chain_lands_three_steps_after_the_limit():
    law = TransfiniteLaw(chain(), "out", "a")
    assert law.at_most(OMEGA + 2) == 0
    assert law.exactly(OMEGA) == 0 and law.loops_into(OMEGA) == {}
    assert law.exactly(OMEGA + 3) == 1


def test_the_swap_jumps_only_at_omega_squared():
    law = TransfiniteLaw(swap(), "out", "a")
    assert all(law.at_most(stage) == 0 for stage in ordinals(2, 4))
    assert law.exactly(Ordinal.omega_power(2)) == 1
    ((attractor, weight),) = law.loops_into(Ordinal.omega_power(2)).items()
    assert attractor.level == 1 and weight == 1


@pytest.mark.parametrize("seed", SEEDS)
def test_the_law_only_grows(seed):
    u = random_universe(seed)
    for event in random_events(u, seed):
        for q in u.states:
            law = TransfiniteLaw(u, event, q)
            previous = Fraction(0)
            for stage in ordinals(u.depth + 1, 3):
                current = law.at_most(stage)
                assert previous <= law.before(stage) <= current <= 1
                assert law.exactly(stage) == current - law.before(stage)
                previous = current


@pytest.mark.parametrize("seed", SEEDS)
def test_limit_jumps_come_from_strange_loops(seed):
    """SPEC Theorem 5.2."""
    u = random_universe(seed)
    for event in random_events(u, seed):
        for q in u.states:
            law = TransfiniteLaw(u, event, q)
            for stage in ordinals(u.depth + 2, 2):
                if stage.is_limit:
                    jumps = law.loops_into(stage)
                    assert law.exactly(stage) == sum(jumps.values(), Fraction(0))
                    for attractor in jumps:
                        assert not attractor.members & event and not attractor.collapses


@pytest.mark.parametrize("seed", SEEDS)
def test_collapsing_limits_never_jump(seed):
    u = random_universe(seed, depth=0)
    taller = Universe({q: u.step(q) for q in u.states}, [None, None])
    for event in random_events(u, seed):
        for q in u.states:
            law = TransfiniteLaw(taller, event, q)
            for stage in ordinals(3, 2):
                if stage.is_limit:
                    assert law.exactly(stage) == 0 and law.loops_into(stage) == {}


@pytest.mark.parametrize("seed", SEEDS)
def test_without_strange_loops_nothing_happens_after_the_finite_stages(seed):
    """SPEC Corollary 5.3: the law is complete at ``ω`` and guarantees stop there."""
    u = random_universe(seed, depth=0)
    taller = Universe({q: u.step(q) for q in u.states}, [None, None])
    for event in random_events(u, seed):
        for q in u.states:
            law = TransfiniteLaw(taller, event, q)
            finite = law.before(OMEGA)
            assert all(law.at_most(stage) == finite for stage in ordinals(3, 2) if stage >= OMEGA)
            guarantee = law.guarantee()
            assert guarantee is INFINITY or guarantee <= OMEGA


@pytest.mark.parametrize("seed", SEEDS)
def test_before_a_limit_in_deterministic_universes(seed):
    """SPEC Proposition 5.1: a deterministic run meets the event within n steps of its level, or never."""
    u = deterministic_universe(seed)
    n = len(u)
    for event in random_events(u, seed):
        for q in u.states:
            law = TransfiniteLaw(u, event, q)
            stages = list(ordinals(u.depth + 1, 2)) + [Ordinal.omega_power(u.depth + 1)]
            for stage in stages:
                if stage.is_limit:
                    k, gamma = split_limit(stage)
                    assert law.before(stage) == law.at_most(gamma + Ordinal.omega_power(k - 1) * n)


@pytest.mark.parametrize("seed", SEEDS)
def test_before_a_limit_is_squeezed_by_the_run_below_it(seed):
    """SPEC Proposition 5.1: before the ``m``-th step of the level below, the run has met
    the event or is still transient, or it never meets it."""
    u = random_universe(seed)
    for event in random_events(u, seed):
        stopped = u.stop_at(event)
        for q in u.states:
            law = TransfiniteLaw(u, event, q)
            for stage in list(ordinals(u.depth + 1, 1))[1:] + [Ordinal.omega_power(u.depth + 1)]:
                if not stage.is_limit:
                    continue
                k, gamma = split_limit(stage)
                settled = {s for a in stopped.attractors(k - 1) for s in a}
                for m in range(0, 25, 4):
                    below = stopped.law(q, gamma + Ordinal.omega_power(k - 1) * m)
                    transient = sum(p for s, p in below.items() if s not in settled)
                    assert mass(below, event) <= law.before(stage) <= mass(below, event) + transient


@pytest.mark.parametrize("seed", SEEDS)
def test_never_is_what_the_last_rule_leaves(seed):
    u = random_universe(seed)
    limit = Ordinal.omega_power(u.depth + 1)
    later = Ordinal.omega_power(u.depth + 3) * 2 + OMEGA + 5
    for event in random_events(u, seed):
        for q in u.states:
            law = TransfiniteLaw(u, event, q)
            assert law.never() == 1 - law.at_most(limit) == 1 - law.at_most(later)
            assert law.guarantee() == hyperprobability(u, event, q)
