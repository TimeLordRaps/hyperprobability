import pytest
from support import HALF, clock, coin, loop, mass, random_events, random_universe

from hyperprobability import OMEGA, TransfiniteLaw, frequency, probability, reach_probability

SEEDS = range(60)


def test_examples():
    u = coin()
    assert probability(u, "H", "T") == 0
    assert probability(u, "H", "T", 1) == probability(u, "H", "T", OMEGA) == HALF
    assert reach_probability(u, "H", "T") == 1 and frequency(u, "H", "T") == HALF
    c = clock()
    assert probability(c, "tick", "tick", 3) == 0 and probability(c, "tick", "tick", 4) == 1
    assert probability(c, "tick", "tick", OMEGA) == frequency(c, "tick", "tock") == HALF
    u = loop()
    assert probability(u, "out", "a", OMEGA) == 1
    assert reach_probability(u, "out", "a") == frequency(u, "out", "a") == 0


@pytest.mark.parametrize("seed", SEEDS)
def test_finite_stages_take_ordinary_steps(seed):
    u = random_universe(seed)
    for event in random_events(u, seed):
        for q in u.states:
            law = {q: 1}
            for n in range(8):
                assert probability(u, event, q, n) == mass(law, event)
                after: dict = {}
                for r, p in law.items():
                    for t, s in u.step(r).items():
                        after[t] = after.get(t, 0) + p * s
                law = after


@pytest.mark.parametrize("seed", SEEDS)
def test_the_in_universe_law_is_the_long_run_average(seed):
    """The Cesaro limit of the step laws is the one family of laws that is invariant,
    harmonic, and inside its own attractor for a start in an attractor."""
    u = random_universe(seed)
    home = {q: a.members for a in u.attractors() for q in a}
    laws = {q: u.in_uni(q) for q in u.states}
    for q, law in laws.items():
        assert sum(law.values()) == 1
        assert u.law(law, 1) == law
        mixed: dict = {}
        for r, p in u.step(q).items():
            for t, s in laws[r].items():
                mixed[t] = mixed.get(t, 0) + p * s
        assert mixed == law
        if q in home:
            assert set(law) <= home[q]
        for event in random_events(u, seed):
            assert frequency(u, event, q) == mass(law, event)


@pytest.mark.parametrize("seed", SEEDS)
def test_reaching_an_event_is_the_law_before_omega(seed):
    u = random_universe(seed)
    for event in random_events(u, seed):
        for q in u.states:
            law = TransfiniteLaw(u, event, q)
            reach = reach_probability(u, event, q)
            assert reach == law.before(OMEGA) == law.at_most(OMEGA) - law.exactly(OMEGA)
            assert all(law.at_most(n) <= reach for n in range(10))
