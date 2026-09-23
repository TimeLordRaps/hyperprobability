import pickle
import random
from fractions import Fraction

import pytest
from support import HALF, chain, clock, coin, dead_end, loop, ordinals, swap

from hyperprobability import (
    OMEGA,
    ONE,
    ZERO,
    OmegaFraction,
    Ordinal,
    Universe,
    frequency,
    kac,
    ordinal_probability,
    reach_probability,
)

W = OmegaFraction.from_ordinal(OMEGA)
# Past every root of the polynomials below, so each rational function has the sign there
# that it has at infinity.
LARGE = 10**12


def value_at(x: OmegaFraction, n: int) -> Fraction:
    top = sum(c * n**i for i, c in enumerate(x.numerator))
    bottom = sum(c * n**i for i, c in enumerate(x.denominator))
    return Fraction(top, bottom)


def random_element(rng: random.Random) -> OmegaFraction:
    def poly() -> list[int]:
        return [rng.randint(-5, 5) for _ in range(rng.randint(0, 3))]

    while True:
        bottom = poly()
        if any(bottom):
            return OmegaFraction(poly(), bottom)


@pytest.mark.parametrize("seed", range(60))
def test_arithmetic_agrees_with_evaluation(seed):
    """Evaluation at a large integer is a field homomorphism wherever it is defined."""
    rng = random.Random(seed)
    x, y = random_element(rng), random_element(rng)
    for n in (LARGE, LARGE + 7):
        assert value_at(x + y, n) == value_at(x, n) + value_at(y, n)
        assert value_at(x - y, n) == value_at(x, n) - value_at(y, n)
        assert value_at(x * y, n) == value_at(x, n) * value_at(y, n)
        assert value_at(-x, n) == -value_at(x, n)
        if y:
            assert value_at(x / y, n) == value_at(x, n) / value_at(y, n)
            assert value_at(y**-2, n) == value_at(y, n) ** -2
    assert (x < y) == (value_at(y - x, LARGE) > 0)
    assert (x == y) == (value_at(x, LARGE) == value_at(y, LARGE))


@pytest.mark.parametrize("seed", range(60))
def test_an_ordered_field(seed):
    rng = random.Random(seed)
    x, y, z = (random_element(rng) for _ in range(3))
    zero, one = OmegaFraction(), OmegaFraction(1)
    assert x + y == y + x and x * y == y * x
    assert (x + y) + z == x + (y + z) and (x * y) * z == x * (y * z)
    assert x * (y + z) == x * y + x * z
    assert x + zero == x and x * one == x and x - x == zero and x + (-x) == 0
    if x:
        assert x * x.reciprocal() == one == x / x
    assert [x < y, x == y, x > y].count(True) == 1
    assert (x <= y) == (x < y or x == y) and (x >= y) == (y <= x)
    if x < y:
        assert x + z < y + z
        if z > 0:
            assert x * z < y * z
    if x > 0 and y > 0:
        assert x * y > 0 and x + y > 0
    assert abs(x) >= 0 and abs(x) * abs(y) == abs(x * y)


def test_iota_embeds_the_ordinals_with_their_natural_operations():
    alphas = list(ordinals(3, 2))
    for alpha in alphas:
        a = OmegaFraction.from_ordinal(alpha)
        assert a.is_ordinal and a.as_ordinal() == alpha
        assert a == alpha and alpha == a and hash(a) == hash(alpha)
        assert str(a) == str(alpha)
        for beta in alphas:
            b = OmegaFraction.from_ordinal(beta)
            assert (a < b) == (alpha < beta) and (a < beta) == (alpha < b)
            assert a + b == alpha.natural_add(beta)
            assert a * b == alpha.natural_mul(beta)


def test_iota_does_not_preserve_ordinary_addition():
    assert ONE + OMEGA == OMEGA
    assert 1 + W > W and OmegaFraction(1) + OMEGA > OMEGA
    assert OMEGA + ONE == W + 1
    assert W - OMEGA == 0


def test_infinite_and_infinitesimal():
    assert W > 10**100 and -W < -(10**100)
    assert 0 < 1 / W < Fraction(1, 10**100)
    assert sorted([W, 1, HALF, 1 / W, -W, 0]) == [-W, 0, 1 / W, HALF, 1, W]
    x = OmegaFraction((5, 2), (1, 4))
    assert x.numerator == (5, 2) and x.denominator == (1, 4)
    assert x.degree == 0 and not x.is_infinite and not x.is_infinitesimal
    assert x.standard_part() == HALF
    assert W.degree == 1 and W.is_infinite and not W.is_infinitesimal
    assert (1 / W).degree == -1 and (1 / W).is_infinitesimal and (1 / W).standard_part() == 0
    assert OmegaFraction().is_infinitesimal and OmegaFraction().standard_part() == 0


def test_lowest_terms():
    assert OmegaFraction((HALF, Fraction(1, 3)), 1) == OmegaFraction((3, 2), 6)
    assert OmegaFraction((HALF, Fraction(1, 3)), 1).numerator == (3, 2)
    assert OmegaFraction(["1/2", "1/3"], 1).denominator == (6,)
    assert OmegaFraction((2, 4), (6, 2)).numerator == (1, 2)
    assert OmegaFraction((1, 1), (1, 1)) == 1
    assert OmegaFraction((-1, 0, 1), (1, 1)) == W - 1
    negative = OmegaFraction(1, (0, -1))
    assert negative == -(1 / W) and negative.numerator == (-1,) and negative.denominator == (0, 1)
    assert OmegaFraction(OMEGA * 2 + 1) == 2 * W + 1 and OmegaFraction(1, OMEGA) == 1 / W
    assert OmegaFraction(0, (1, 2)).denominator == (1,)
    assert W**3 == Ordinal.omega_power(3) and W**0 == 1 and W**-1 == 1 / W
    assert (W + 1) ** 2 == W**2 + 2 * W + 1 and (2 * W) ** -2 == OmegaFraction(1, (0, 0, 4))
    assert W + HALF == OmegaFraction((1, 2), 2) and HALF * W == OmegaFraction((0, 1), 2)
    assert HALF / W == OmegaFraction(1, (0, 2)) and OMEGA * W == W**2 == W * OMEGA


@pytest.mark.parametrize(
    "value, text",
    [
        (1 / W, "1/ω"),
        (1 / (2 * W + 3), "1/(ω*2 + 3)"),
        ((2 * W + 1) / (2 * W**2), "(ω*2 + 1)/(ω^2*2)"),
        (OmegaFraction(-3, 4), "-3/4"),
        (1 - W, "-ω + 1"),
        (OmegaFraction(), "0"),
        (3 * W**2 / (2 * W + 2), "ω^2*3/(ω*2 + 2)"),
        ((5 - 6 * W**2) / (15 * W + 15), "(-ω^2*6 + 5)/(ω*15 + 15)"),
        (W**-2, "1/ω^2"),
        (W / 3, "ω/3"),
    ],
)
def test_printing(value, text):
    assert str(value) == text
    assert eval(repr(value), {"OmegaFraction": OmegaFraction}) == value
    assert pickle.loads(pickle.dumps(value)) == value


def test_repr():
    assert repr(1 / W) == "OmegaFraction((1,), (0, 1))"
    assert repr(OmegaFraction()) == "OmegaFraction((), (1,))"


def test_equal_values_hash_alike():
    assert OmegaFraction(3) == 3 == Ordinal.from_int(3) and hash(OmegaFraction(3)) == hash(3)
    assert OmegaFraction(1, 2) == HALF and hash(OmegaFraction(1, 2)) == hash(HALF)
    assert OmegaFraction(-3) == -3 and hash(OmegaFraction(-3)) == hash(-3)
    assert W == OMEGA and hash(W) == hash(OMEGA)
    assert len({OmegaFraction(2), 2, Fraction(2)}) == 1
    assert len({OmegaFraction(2), Ordinal.from_int(2)}) == 1
    assert OmegaFraction() == 0 == ZERO and hash(OmegaFraction()) == hash(0)
    assert W != 3 and W != "ω" and W != 1.0


def test_errors():
    for make in (
        lambda: OmegaFraction(1, 0),
        lambda: OmegaFraction(1, [0, 0]),
        lambda: OmegaFraction().reciprocal(),
        lambda: W / 0,
        lambda: OmegaFraction() ** -1,
    ):
        with pytest.raises(ZeroDivisionError):
            make()
    for bad in (1.5, True, object(), OmegaFraction(1), [1.5], complex(1, 0)):
        with pytest.raises(TypeError):
            OmegaFraction(bad)
        with pytest.raises(TypeError):
            OmegaFraction(1, bad)
    with pytest.raises(ValueError):
        OmegaFraction("one")
    for operation in (
        lambda: W < 1.5,
        lambda: W + 1.5,
        lambda: W * True,
        lambda: W**0.5,
        lambda: W**W,
        lambda: W < "a",
    ):
        with pytest.raises(TypeError):
            operation()
    for operation in (
        lambda: (-W).as_ordinal(),
        lambda: (1 / W).as_ordinal(),
        lambda: OmegaFraction().degree,
        lambda: W.standard_part(),
        lambda: OmegaFraction.from_ordinal(-1),
    ):
        with pytest.raises(ValueError):
            operation()


@pytest.mark.parametrize(
    "universe, event, start, expected",
    [
        (clock(), "tick", "tick", OmegaFraction(1, 2)),
        (coin(), "H", "H", 1 / W),
        (loop(), "out", "a", 1 / W),
        (chain(), "out", "a", 1 / (W + 3)),
        (swap(), "out", "a", W**-2),
        (dead_end(), "out", "a", OmegaFraction()),
    ],
)
def test_ordinal_probability(universe, event, start, expected):
    assert ordinal_probability(universe, event, start) == expected


def test_ordinal_probability_sees_what_in_universe_probability_cannot():
    u = loop()
    p = ordinal_probability(u, "out", "a")
    assert 0 < p < Fraction(1, 10**100) and p.standard_part() == 0
    assert reach_probability(u, "out", "a") == frequency(u, "out", "a") == 0


def test_ordinal_probability_is_not_a_measure():
    """SPEC Proposition 7.5: neither additive nor maxitive, and only a law spread over an
    attractor's part of the event is held below its frequency."""
    u = clock()
    tick, tock = ordinal_probability(u, "tick", "tick"), ordinal_probability(u, "tock", "tick")
    assert tick + tock == Fraction(3, 2) and ordinal_probability(u, ["tick", "tock"], "tick") == 1
    assert ordinal_probability(u, [], "tick") == 0
    u = coin()
    heads, tails = ordinal_probability(u, "H", "H"), ordinal_probability(u, "T", "H")
    assert heads == tails == 1 / W and ordinal_probability(u, ["H", "T"], "H") == 1
    u = Universe({"e1": {"e2": 1}, "e2": {"x": 1}, "x": {"e1": HALF, "x": HALF}})
    event = ["e1", "e2"]
    assert ordinal_probability(u, event, "e1") == 1 > frequency(u, event, "e1") == HALF
    (bound,) = kac(u, event)
    spread = ordinal_probability(u, event, {"e1": HALF, "e2": HALF})
    assert spread == bound.ordinal_probability == 1 / W and bound.holds and not bound.tight
