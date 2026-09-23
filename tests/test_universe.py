import random
from fractions import Fraction

import pytest
from support import HALF, clock, coin, loop, random_events, random_universe, swap

from hyperprobability import OMEGA, Ordinal, Universe

SEEDS = range(60)


def random_ordinal(rng: random.Random, exponents: int) -> Ordinal:
    return Ordinal(tuple(rng.randint(0, 3) for _ in range(rng.randint(0, exponents))))


def random_limit(rng: random.Random, exponents: int) -> Ordinal:
    while True:
        alpha = Ordinal((0,) + tuple(rng.randint(0, 2) for _ in range(exponents - 1)))
        if alpha:
            return alpha


@pytest.mark.parametrize(
    "kernel, error",
    [
        ([("a", {"a": 1})], TypeError),
        ({}, ValueError),
        ({"a": ["a"]}, TypeError),
        ({"a": {"a": HALF}}, ValueError),
        ({"a": {"a": 2, "b": -1}, "b": {"b": 1}}, ValueError),
        ({"a": {"b": 1}}, ValueError),
        ({"a": {"a": 1.0}}, TypeError),
        ({"a": {"a": True}}, TypeError),
        ({"a": {"a": "one"}}, ValueError),
    ],
)
def test_a_kernel_must_be_contained_and_exact(kernel, error):
    with pytest.raises(error):
        Universe(kernel)


def test_rows_accept_exact_strings():
    u = Universe({"a": {"a": "1/3", "b": "2/3"}, "b": {"a": "0.5", "b": "1/2"}})
    assert u.step("a") == {"a": Fraction(1, 3), "b": Fraction(2, 3)}
    assert u.step("b") == {"a": HALF, "b": HALF}


def test_limit_rules_are_validated():
    step = {"a": "b", "b": "a"}
    with pytest.raises(TypeError):
        Universe.deterministic(step, {frozenset(step): "a"})
    with pytest.raises(TypeError):
        Universe.deterministic(step, [{("a", "b"): "a"}])
    with pytest.raises(ValueError, match="not a level-0 attractor"):
        Universe.deterministic(step, [{frozenset({"a"}): "a"}])
    with pytest.raises(TypeError):
        Universe.deterministic(step, [5])
    with pytest.raises(TypeError):
        Universe.deterministic(step, [{frozenset(step): 3.5}])
    with pytest.raises(ValueError):
        Universe.deterministic(step, [{frozenset(step): {"a": HALF}}])
    with pytest.raises(TypeError):
        Universe.deterministic(["a"])


def test_accessors():
    u = clock()
    assert u.states == ("tick", "tock")
    assert len(u) == 2 and "tick" in u and "noon" not in u and [] not in u
    assert u.depth == 0
    assert u.step("tick") == {"tock": 1}
    with pytest.raises(ValueError):
        u.step("noon")
    assert u.kernel() == {"tick": {"tock": 1}, "tock": {"tick": 1}}
    both = {"tick": HALF, "tock": HALF}
    assert u.kernel(1) == {"tick": both, "tock": both}
    with pytest.raises(ValueError):
        u.kernel(-1)
    with pytest.raises(TypeError):
        u.kernel(True)
    assert repr(u) == "Universe(states=('tick', 'tock'), depth=0)"


def test_attractors():
    (cycle,) = clock().attractors()
    assert (cycle.level, cycle.states, cycle.period, cycle.kind) == (
        0,
        ("tick", "tock"),
        2,
        "periodic_cycle",
    )
    assert cycle.deterministic and cycle.collapses
    assert dict(cycle.law) == dict(cycle.exit) == {"tick": HALF, "tock": HALF}
    assert "tick" in cycle and "noon" not in cycle and [] not in cycle
    assert list(cycle) == ["tick", "tock"] and len(cycle) == 2
    assert cycle.members == frozenset({"tick", "tock"})
    assert repr(cycle) == "Attractor(level=0, states=('tick', 'tock'), kind='periodic_cycle')"

    (toss,) = coin().attractors()
    assert (toss.kind, toss.period, toss.deterministic) == ("recurrent_class", 1, False)

    u = swap()
    first, second, out = u.attractors(0)
    assert first.kind == "fixed_point" and dict(first.exit) == {"b": 1} and not first.collapses
    assert out.collapses
    pair, rest = u.attractors(1)
    assert pair.states == ("a", "b") and pair.period == 2 and dict(pair.exit) == {"out": 1}
    assert [a.states for a in u.attractors(2)] == [("out",)]
    assert u.attractors(0) == u.attractors(0)
    assert {a: 1 for a in u.attractors(0)}[first] == 1


def test_laws_at_transfinite_stages():
    u = swap()
    assert u.law("a") == {"a": 1}
    assert u.law("a", 5) == {"a": 1}
    assert u.law("a", OMEGA) == {"b": 1}
    assert u.law("a", OMEGA * 2) == {"a": 1}
    assert u.law("a", OMEGA * 3 + 4) == {"b": 1}
    assert u.law("a", Ordinal.omega_power(2)) == {"out": 1}
    assert u.law({"a": HALF, "b": HALF}, OMEGA) == {"a": HALF, "b": HALF}
    assert (
        clock().law("tick", OMEGA) == clock().law("tick", OMEGA + 1) == {"tick": HALF, "tock": HALF}
    )
    with pytest.raises(TypeError):
        u.law("a", 1.5)
    with pytest.raises(ValueError):
        u.law("a", -1)
    with pytest.raises(ValueError):
        u.law("nowhere")
    with pytest.raises(ValueError):
        u.law({"a": HALF})


def test_the_law_can_move_after_the_kernels_settle():
    """SPEC Proposition 2.4: past the last rule the kernels stop changing, but the law need not."""
    u = Universe.deterministic(
        {"a": "a", "b": "x", "x": "x"}, [{frozenset({"a"}): "b", frozenset({"x"}): "a"}]
    )
    assert u.kernel(2) == u.kernel(3)
    square = Ordinal.omega_power(2)
    assert u.law("a", square) == u.law("a", square + OMEGA) == {"a": HALF, "b": HALF}
    assert u.law("a", square + 1) == {"a": HALF, "x": HALF}


def test_collapse_names_the_attractors_that_hold_the_tail():
    u = swap()
    ((fixed, weight),) = u.collapse("a", OMEGA).items()
    assert fixed.states == ("a",) and weight == 1
    ((pair, weight),) = u.collapse("a", Ordinal.omega_power(2)).items()
    assert pair.level == 1 and pair.states == ("a", "b") and weight == 1
    for stage in (0, 3, OMEGA + 1):
        with pytest.raises(ValueError):
            u.collapse("a", stage)


def test_in_uni_ignores_strange_loops():
    assert loop().in_uni("a") == {"a": HALF, "b": HALF}
    assert loop().law("a", OMEGA) == {"out": 1}
    assert coin().in_uni("H") == {"H": HALF, "T": HALF}


def test_events_and_initial_laws_are_read_strictly():
    u = clock()
    assert u.stop_at("tick").step("tick") == {"tick": 1}
    assert u.stop_at(lambda state: state.startswith("to")).step("tock") == {"tock": 1}
    assert u.stop_at(["tick", "tock"]).kernel() == {"tick": {"tick": 1}, "tock": {"tock": 1}}
    with pytest.raises(ValueError):
        u.stop_at(["noon"])
    with pytest.raises(TypeError):
        u.stop_at(5)


@pytest.mark.parametrize("seed", SEEDS)
def test_higher_kernels_absorb_lower_ones(seed):
    """``K_j K_k = K_k`` for ``j < k``, and the levels stop changing past the last rule."""
    u = random_universe(seed)
    top = u.depth + 2
    kernels = [u.kernel(k) for k in range(top + 1)]
    for rows in kernels:
        assert all(sum(row.values()) == 1 for row in rows.values())
    for k in range(1, top + 1):
        for j in range(k):
            for q in u.states:
                after: dict = {}
                for r, p in kernels[j][q].items():
                    for t, s in kernels[k][r].items():
                        after[t] = after.get(t, 0) + p * s
                assert {t: v for t, v in after.items() if v} == kernels[k][q]
    assert u.kernel(u.depth + 1) == u.kernel(u.depth + 5)


@pytest.mark.parametrize("seed", SEEDS)
def test_composition(seed):
    """SPEC Theorem 2.3: the law at ``α + β`` is the law at ``β`` from the law at ``α``."""
    u = random_universe(seed)
    rng = random.Random(seed)
    for _ in range(8):
        alpha = random_ordinal(rng, u.depth + 2)
        beta = random_ordinal(rng, u.depth + 2)
        start = rng.choice(u.states)
        assert u.law(u.law(start, alpha), beta) == u.law(start, alpha + beta)


@pytest.mark.parametrize("seed", SEEDS)
def test_limit_laws_mix_the_exits_of_the_collapsing_attractors(seed):
    u = random_universe(seed)
    rng = random.Random(seed)
    for _ in range(5):
        stage = random_limit(rng, u.depth + 2)
        start = rng.choice(u.states)
        weights = u.collapse(start, stage)
        assert sum(weights.values()) == 1
        mixture: dict = {}
        for attractor, weight in weights.items():
            for state, p in attractor.exit.items():
                mixture[state] = mixture.get(state, 0) + weight * p
        assert mixture == u.law(start, stage)


@pytest.mark.parametrize("seed", SEEDS)
def test_the_tail_settles_where_the_run_below_the_limit_is_heading(seed):
    """At ``γ + ω**k`` the weight of an attractor lies between the chance that the level
    ``k - 1`` run from ``γ`` is already in it and that chance plus the transient mass."""
    u = random_universe(seed)
    rng = random.Random(seed)
    for _ in range(5):
        stage = random_limit(rng, u.depth + 2)
        k = next(i for i, c in enumerate(stage.coefficients) if c)
        gamma = Ordinal(tuple(c - (i == k) for i, c in enumerate(stage.coefficients)))
        start = rng.choice(u.states)
        weights = u.collapse(start, stage)
        attractors = u.attractors(k - 1)
        settled = {q for a in attractors for q in a}
        for m in range(0, 30, 6):
            below = u.law(start, gamma + Ordinal.omega_power(k - 1) * m)
            transient = sum(p for q, p in below.items() if q not in settled)
            for attractor in attractors:
                inside = sum(p for q, p in below.items() if q in attractor)
                assert inside <= weights.get(attractor, 0) <= inside + transient


@pytest.mark.parametrize("seed", SEEDS)
def test_stopped_attractors_that_miss_the_event_are_attractors(seed):
    """SPEC Lemma 4.2, with the same rows, stationary laws and exits."""
    u = random_universe(seed)
    for event in random_events(u, seed):
        stopped = u.stop_at(event)
        for k in range(u.depth + 2):
            originals = {a.members: a for a in u.attractors(k)}
            rows, stopped_rows = u.kernel(k), stopped.kernel(k)
            for attractor in stopped.attractors(k):
                if attractor.members & event:
                    assert len(attractor) == 1 and dict(attractor.exit) == dict(attractor.law)
                    continue
                original = originals[attractor.members]
                assert dict(original.law) == dict(attractor.law)
                assert dict(original.exit) == dict(attractor.exit)
                assert all(rows[q] == stopped_rows[q] for q in attractor)


def test_an_attractor_that_misses_the_event_need_not_survive_stopping():
    """The converse of SPEC Lemma 4.2(b) fails (SPEC Theorem 5.2)."""
    u = Universe.deterministic({"r": "e", "e": "s", "s": "s"}, [{frozenset({"s"}): "r"}])
    assert [a.states for a in u.attractors(1)] == [("r",)]
    assert [a.states for a in u.stop_at("e").attractors(1)] == [("e",)]


@pytest.mark.parametrize("seed", SEEDS)
def test_in_uni_is_the_collapsed_limit(seed):
    u = random_universe(seed, depth=0)
    for q in u.states:
        law = u.in_uni(q)
        assert law == u.law(q, OMEGA)
        assert u.law(law, 1) == law
