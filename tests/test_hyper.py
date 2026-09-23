import pickle
import random
from fractions import Fraction

import pytest
from support import (
    HALF,
    chain,
    clock,
    coin,
    dead_end,
    first_certain,
    first_certain_return,
    fork,
    loop,
    random_events,
    random_universe,
    swap,
)

from hyperprobability import (
    INFINITY,
    OMEGA,
    ONE,
    ZERO,
    Infinity,
    Ordinal,
    TransfiniteLaw,
    Universe,
    hyperprobabilities,
    hyperprobability,
    return_hyperprobability,
)

SEEDS = range(80)
N = Ordinal.from_int


@pytest.mark.parametrize(
    "universe, event, start, expected",
    [
        (clock(), "tock", "tick", N(1)),
        (clock(), "tick", "tick", ZERO),
        (coin(), "H", "T", OMEGA),
        (loop(), "out", "a", OMEGA),
        (chain(), "out", "a", OMEGA + 3),
        (chain(), "out", "x1", N(3)),
        (fork(), "out", "a", OMEGA + 2),
        (swap(), "out", "a", Ordinal.omega_power(2)),
        (swap(), "b", "a", OMEGA),
        (dead_end(), "out", "a", INFINITY),
    ],
)
def test_examples(universe, event, start, expected):
    assert hyperprobability(universe, event, start) == expected
    assert first_certain(universe, frozenset({event}), start) == expected


def test_return_guarantees():
    assert return_hyperprobability(clock(), "tick", "tick") == N(2)
    assert return_hyperprobability(coin(), "H", "H") == OMEGA
    assert return_hyperprobability(loop(), "out", "a") == OMEGA
    assert return_hyperprobability(dead_end(), "a", "a") == N(1)
    assert return_hyperprobability(dead_end(), "out", "a") is INFINITY


def test_every_state_at_once():
    assert hyperprobabilities(swap(), "out") == {
        "a": Ordinal.omega_power(2),
        "b": Ordinal.omega_power(2),
        "out": ZERO,
    }
    assert hyperprobabilities(coin(), []) == {"H": INFINITY, "T": INFINITY}


def test_an_initial_law_takes_the_worst_state():
    u = chain()
    assert hyperprobability(u, "out", {"a": HALF, "x2": HALF}) == OMEGA + 3
    assert hyperprobability(u, "out", {"x1": HALF, "x2": HALF}) == N(3)
    assert hyperprobability(u, "out", {"x2": HALF, "out": HALF}) == N(2)


def test_infinity():
    assert Infinity() is INFINITY
    assert pickle.loads(pickle.dumps(INFINITY)) is INFINITY
    assert str(INFINITY) == "∞" and repr(INFINITY) == "INFINITY"
    assert INFINITY == INFINITY and INFINITY != OMEGA and INFINITY != 3
    for value in (ZERO, OMEGA, Ordinal.omega_power(9), 7):
        assert value < INFINITY and INFINITY > value and value <= INFINITY and INFINITY >= value
        assert not INFINITY < value and not INFINITY <= value
    assert INFINITY <= INFINITY and INFINITY >= INFINITY and not INFINITY < INFINITY
    assert max([OMEGA, INFINITY, ONE]) is INFINITY
    assert INFINITY + ONE is INFINITY and ONE + INFINITY is INFINITY
    assert len({INFINITY, Infinity()}) == 1
    for other in ("a", 1.5, True):
        with pytest.raises(TypeError):
            INFINITY < other


@pytest.mark.parametrize("block", range(40))
def test_the_recursion_finds_the_least_certain_stage(block):
    """SPEC Theorem 4.4 against exhaustive search over stages, in 2000 universes.

    So many because some mistakes show up only rarely: taking the best rather
    than the worst state after a limit that forks changes about one universe
    in three hundred.
    """
    for seed in range(50 * block, 50 * block + 50):
        u = random_universe(seed)
        for event in random_events(u, seed):
            values = hyperprobabilities(u, event)
            for q in u.states:
                assert values[q] == first_certain(u, event, q), (seed, sorted(event), q)


@pytest.mark.parametrize("block", range(10))
def test_return_guarantees_against_search(block):
    """SPEC Proposition 4.8."""
    for seed in range(50 * block, 50 * block + 50):
        u = random_universe(seed)
        for event in random_events(u, seed):
            for q in u.states:
                expected = first_certain_return(u, event, q)
                assert return_hyperprobability(u, event, q) == expected, (seed, sorted(event), q)


@pytest.mark.parametrize("seed", SEEDS)
def test_initial_laws_take_the_worst_state_of_their_support(seed):
    """SPEC Proposition 4.3."""
    u = random_universe(seed)
    rng = random.Random(seed)
    for event in random_events(u, seed):
        states = rng.sample(list(u.states), rng.randint(1, len(u)))
        law = {q: Fraction(1, len(states)) for q in states}
        expected = max(hyperprobability(u, event, q) for q in states)
        assert hyperprobability(u, event, law) == expected == first_certain(u, event, law)


@pytest.mark.parametrize("seed", SEEDS)
def test_guarantees_are_monotone_in_the_event(seed):
    """SPEC Proposition 4.6."""
    u = random_universe(seed)
    events = random_events(u, seed, count=5)
    for small in events:
        for large in events:
            if small <= large:
                for q in u.states:
                    assert hyperprobability(u, large, q) <= hyperprobability(u, small, q)


@pytest.mark.parametrize("seed", SEEDS)
def test_guarantees_are_infinite_exactly_when_the_event_can_be_missed(seed):
    u = random_universe(seed)
    for event in random_events(u, seed):
        for q in u.states:
            never = TransfiniteLaw(u, event, q).never()
            assert (hyperprobability(u, event, q) is INFINITY) == (never > 0)


@pytest.mark.parametrize("seed", SEEDS)
def test_guarantees_stay_below_the_first_collapsed_level(seed):
    """SPEC Corollary 4.5: at most ``ω**(depth+1)``, with every coefficient below the
    number of states."""
    u = random_universe(seed)
    bound = Ordinal.omega_power(u.depth + 1)
    for event in random_events(u, seed):
        for value in hyperprobabilities(u, event).values():
            assert value is INFINITY or value <= bound
            if value is not INFINITY:
                assert all(c < len(u) for c in value.coefficients)


@pytest.mark.parametrize("seed", SEEDS)
def test_guarantees_depend_only_on_what_is_possible(seed):
    """SPEC Proposition 4.9: reweighting every row and every limit target, keeping which
    states they can reach, changes no guarantee."""
    u = random_universe(seed)
    rng = random.Random(seed)

    def reweight(law: dict) -> dict:
        weights = {state: rng.randint(1, 9) for state in law}
        total = sum(weights.values())
        return {state: Fraction(weight, total) for state, weight in weights.items()}

    kernel = {q: reweight(u.step(q)) for q in u.states}
    rules = [
        {a.members: reweight(dict(a.exit)) for a in u.attractors(k) if not a.collapses}
        for k in range(u.depth)
    ]
    v = Universe(kernel, rules)
    for event in random_events(u, seed):
        assert hyperprobabilities(v, event) == hyperprobabilities(u, event)


@pytest.mark.parametrize("seed", SEEDS)
def test_collapsing_levels_change_nothing(seed):
    """Once a guarantee exists, higher collapsing levels keep it (SPEC Theorem 4.4)."""
    u = random_universe(seed)
    kernel = {q: u.step(q) for q in u.states}
    rules = [u._rules[k] for k in range(u.depth)]
    taller = Universe(kernel, rules + [None, None])
    for event in random_events(u, seed):
        assert hyperprobabilities(taller, event) == hyperprobabilities(u, event)
