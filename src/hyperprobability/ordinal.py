"""Ordinal probability: exact fractions of ordinals in the ordered field ℚ(ω).

Hyperprobabilities are ordinals below ``ω**ω``. In Cantor normal form such an
ordinal is a polynomial in ``ω`` with natural-number coefficients, and its
natural (Hessenberg) sum and product are exactly polynomial sum and product.
The map ``ι`` that reads an ordinal as its polynomial therefore embeds the
ordinals in ℚ(ω), the field of rational functions in ``ω`` ordered by their
behaviour as ``ω`` grows past every rational.

Ordinal probability lives in that field: ``p_μ(E) = 1/ι(ℍ⁺_μ(E))``, one over
the number of evaluations that surely produce the event. It is built from a
hyperprobability, which must exist first; that is the sense in which ordinal
probability is the higher-order notion.

``ι`` preserves order, natural sum and natural product. It does not preserve
ordinary ordinal addition: ``1 + ω == ω`` as ordinals, but ``1 + ι(ω) > ι(ω)``
in the field. See SPEC.md §7.
"""

from __future__ import annotations

import math
import numbers
from fractions import Fraction

from ordinatics.ordinals import Ordinal

from hyperprobability.hyper import INFINITY, return_hyperprobability
from hyperprobability.universe import Universe, _stage, exact

# Coefficients of 1, ω, ω**2, ...; rationals while computing, integers once reduced.
Poly = tuple[Fraction, ...]
IntPoly = tuple[int, ...]

_ZERO = Fraction(0)
_ONE = Fraction(1)


def _trim(values: list[Fraction]) -> Poly:
    while values and not values[-1]:
        values.pop()
    return tuple(values)


def _add(a: Poly, b: Poly) -> Poly:
    if len(a) < len(b):
        a, b = b, a
    values = list(a)
    for i, c in enumerate(b):
        values[i] += c
    return _trim(values)


def _neg(a: Poly) -> Poly:
    return tuple(-c for c in a)


def _mul(a: Poly, b: Poly) -> Poly:
    if not a or not b:
        return ()
    values = [_ZERO] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                values[i + j] += x * y
    return _trim(values)


def _scale(a: Poly, factor: Fraction) -> Poly:
    return _trim([c * factor for c in a])


def _divmod(a: Poly, b: Poly) -> tuple[Poly, Poly]:
    """Division with remainder over ℚ: ``a == q*b + r`` with ``deg r < deg b``."""
    shift = len(a) - len(b)
    if shift < 0:
        return (), a
    remainder = list(a)
    quotient = [_ZERO] * (shift + 1)
    for s in range(shift, -1, -1):
        coefficient = remainder[s + len(b) - 1] / b[-1]
        quotient[s] = coefficient
        if coefficient:
            for i, c in enumerate(b):
                remainder[s + i] -= coefficient * c
    return _trim(quotient), _trim(remainder[: len(b) - 1])


def _gcd(a: Poly, b: Poly) -> Poly:
    """The monic greatest common divisor of two polynomials, not both zero."""
    while b:
        a, b = b, _divmod(a, b)[1]
    return _scale(a, 1 / a[-1])


def _reduce(top: Poly, bottom: Poly) -> tuple[IntPoly, IntPoly]:
    """Lowest terms: coprime integer polynomials, jointly primitive, the denominator's lead positive."""
    top = tuple(Fraction(c) for c in top)
    bottom = tuple(Fraction(c) for c in bottom)
    if not bottom:
        raise ZeroDivisionError("OmegaFraction with a zero denominator")
    if not top:
        return (), (1,)
    common = _gcd(top, bottom)
    if len(common) > 1:
        top = _divmod(top, common)[0]
        bottom = _divmod(bottom, common)[0]
    parts = top + bottom
    scale = math.lcm(*(c.denominator for c in parts))
    content = math.gcd(*(int(c * scale) for c in parts))
    factor = Fraction(scale, content if bottom[-1] > 0 else -content)
    return tuple(int(c * factor) for c in top), tuple(int(c * factor) for c in bottom)


def _coefficients(value: object) -> Poly:
    if isinstance(value, OmegaFraction):
        raise TypeError("divide OmegaFractions with / rather than nesting them")
    if isinstance(value, Ordinal):
        return tuple(Fraction(c) for c in value.coefficients)
    if isinstance(value, (numbers.Number, str)):
        return _trim([exact(value)])
    try:
        items = list(value)  # type: ignore[call-overload]
    except TypeError:
        raise TypeError(
            "each part of an OmegaFraction is an exact rational, an Ordinal, "
            "or a sequence of exact rational coefficients"
        ) from None
    return _trim([exact(item) for item in items])


def _term(coefficient: int, power: int) -> str:
    # The notation of ordinatics: ω^2*3 + ω*2 + 1.
    if power == 0:
        return str(coefficient)
    base = "ω" if power == 1 else f"ω^{power}"
    return base if coefficient == 1 else f"{base}*{coefficient}"


def _terms(poly: IntPoly) -> int:
    return sum(1 for c in poly if c)


def _format(poly: IntPoly) -> str:
    parts: list[str] = []
    for power in range(len(poly) - 1, -1, -1):
        c = poly[power]
        if not c:
            continue
        text = _term(abs(c), power)
        if parts:
            parts.append(("- " if c < 0 else "+ ") + text)
        else:
            parts.append(("-" if c < 0 else "") + text)
    return " ".join(parts) if parts else "0"


def _show(poly: IntPoly) -> str:
    return "(" + ", ".join(str(c) for c in poly) + ("," if len(poly) == 1 else "") + ")"


class OmegaFraction:
    """An exact element of ℚ(ω): a rational function of ``ω``, ordered at infinity.

    ``OmegaFraction(numerator, denominator)`` reads each part as an exact
    rational, an :class:`~ordinatics.ordinals.Ordinal` (through ``ι``), or a
    sequence of exact rationals giving the coefficients of ``1, ω, ω**2, ...``.
    The value is kept in lowest terms: integer polynomials with no common
    factor, no common divisor of all their coefficients, and a positive
    leading coefficient below. Since ``ω`` exceeds every rational, ``x < y``
    exactly when ``y - x`` has a positive leading coefficient. So ``1/ω`` is a
    positive infinitesimal and ``ω`` is infinite.

    Arithmetic and comparison accept ints, fractions and Ordinals, the latter
    through ``ι``. An OmegaFraction equals an Ordinal ``α`` exactly when it is
    ``ι(α)``, and then it hashes like ``α``.
    """

    __slots__ = ("_numerator", "_denominator")

    def __init__(self, numerator: object = 0, denominator: object = 1) -> None:
        self._numerator, self._denominator = _reduce(
            _coefficients(numerator), _coefficients(denominator)
        )

    @classmethod
    def _from_polys(cls, top: Poly, bottom: Poly) -> OmegaFraction:
        result = object.__new__(cls)
        result._numerator, result._denominator = _reduce(top, bottom)
        return result

    @classmethod
    def from_ordinal(cls, alpha: object) -> OmegaFraction:
        """``ι(α)``: an Ordinal or natural number read as the polynomial of its normal form."""
        result = object.__new__(cls)
        result._numerator, result._denominator = _stage(alpha).coefficients, (1,)
        return result

    @property
    def numerator(self) -> IntPoly:
        """Integer coefficients of the numerator, lowest power of ``ω`` first."""
        return self._numerator

    @property
    def denominator(self) -> IntPoly:
        """Integer coefficients of the denominator, lowest power of ``ω`` first."""
        return self._denominator

    @property
    def is_ordinal(self) -> bool:
        """Whether this is ``ι(α)`` for an ordinal ``α`` below ``ω**ω``."""
        return self._denominator == (1,) and all(c >= 0 for c in self._numerator)

    def as_ordinal(self) -> Ordinal:
        """The ordinal ``α`` with ``ι(α)`` equal to this element."""
        if not self.is_ordinal:
            raise ValueError(f"{self} is not the image of an ordinal")
        return Ordinal(self._numerator)

    @property
    def degree(self) -> int:
        """How fast this grows with ``ω``: positive if infinite, negative if infinitesimal."""
        if not self._numerator:
            raise ValueError("zero has no degree")
        return len(self._numerator) - len(self._denominator)

    @property
    def is_infinite(self) -> bool:
        """Whether this exceeds every rational in absolute value."""
        return bool(self._numerator) and self.degree > 0

    @property
    def is_infinitesimal(self) -> bool:
        """Whether this is below every positive rational in absolute value (zero included)."""
        return not self._numerator or self.degree < 0

    def standard_part(self) -> Fraction:
        """The rational that a finite element approaches as ``ω`` grows without bound."""
        if self.is_infinite:
            raise ValueError(f"{self} is infinite and has no standard part")
        if self.is_infinitesimal:
            return _ZERO
        return Fraction(self._numerator[-1], self._denominator[-1])

    def reciprocal(self) -> OmegaFraction:
        """``1/self``."""
        if not self._numerator:
            raise ZeroDivisionError("zero has no reciprocal")
        return OmegaFraction._from_polys(self._denominator, self._numerator)

    def _sign(self) -> int:
        if not self._numerator:
            return 0
        return 1 if self._numerator[-1] > 0 else -1

    def _compare(self, other: object) -> object:
        value = _coerce(other)
        if value is NotImplemented:
            return NotImplemented
        difference = _add(
            _mul(self._numerator, value._denominator),
            _neg(_mul(value._numerator, self._denominator)),
        )
        if not difference:
            return 0
        return 1 if difference[-1] > 0 else -1

    def __add__(self, other: object) -> OmegaFraction:
        value = _coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return OmegaFraction._from_polys(
            _add(
                _mul(self._numerator, value._denominator), _mul(value._numerator, self._denominator)
            ),
            _mul(self._denominator, value._denominator),
        )

    __radd__ = __add__

    def __neg__(self) -> OmegaFraction:
        return OmegaFraction._from_polys(_neg(self._numerator), self._denominator)

    def __pos__(self) -> OmegaFraction:
        return self

    def __abs__(self) -> OmegaFraction:
        return -self if self._sign() < 0 else self

    def __sub__(self, other: object) -> OmegaFraction:
        value = _coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return self + (-value)

    def __rsub__(self, other: object) -> OmegaFraction:
        value = _coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return value + (-self)

    def __mul__(self, other: object) -> OmegaFraction:
        value = _coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return OmegaFraction._from_polys(
            _mul(self._numerator, value._numerator), _mul(self._denominator, value._denominator)
        )

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> OmegaFraction:
        value = _coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return self * value.reciprocal()

    def __rtruediv__(self, other: object) -> OmegaFraction:
        value = _coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return value * self.reciprocal()

    def __pow__(self, exponent: object) -> OmegaFraction:
        if isinstance(exponent, bool) or not isinstance(exponent, numbers.Integral):
            return NotImplemented
        power = int(exponent)
        factor = self if power >= 0 else self.reciprocal()
        power = abs(power)
        result = OmegaFraction(1)
        while power:
            if power & 1:
                result = result * factor
            power >>= 1
            if power:
                factor = factor * factor
        return result

    def __eq__(self, other: object) -> bool:
        value = _coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return self._numerator == value._numerator and self._denominator == value._denominator

    def __lt__(self, other: object) -> bool:
        c = self._compare(other)
        return c if c is NotImplemented else c < 0  # type: ignore[operator,return-value]

    def __le__(self, other: object) -> bool:
        c = self._compare(other)
        return c if c is NotImplemented else c <= 0  # type: ignore[operator,return-value]

    def __gt__(self, other: object) -> bool:
        c = self._compare(other)
        return c if c is NotImplemented else c > 0  # type: ignore[operator,return-value]

    def __ge__(self, other: object) -> bool:
        c = self._compare(other)
        return c if c is NotImplemented else c >= 0  # type: ignore[operator,return-value]

    def __hash__(self) -> int:
        # Equal values hash alike: an ordinal image like its Ordinal, a
        # rational constant like its Fraction.
        if self.is_ordinal:
            return hash(self.as_ordinal())
        if len(self._numerator) == 1 and len(self._denominator) == 1:
            return hash(Fraction(self._numerator[0], self._denominator[0]))
        return hash((OmegaFraction, self._numerator, self._denominator))

    def __bool__(self) -> bool:
        return bool(self._numerator)

    def __reduce__(self) -> tuple:
        return (OmegaFraction, (self._numerator, self._denominator))

    def __str__(self) -> str:
        top = _format(self._numerator)
        if self._denominator == (1,):
            return top
        bottom = _format(self._denominator)
        if _terms(self._numerator) > 1:
            top = f"({top})"
        if _terms(self._denominator) > 1 or "*" in bottom:
            bottom = f"({bottom})"
        return f"{top}/{bottom}"

    def __repr__(self) -> str:
        return f"OmegaFraction({_show(self._numerator)}, {_show(self._denominator)})"


def _coerce(value: object) -> OmegaFraction:
    if isinstance(value, OmegaFraction):
        return value
    if isinstance(value, Ordinal):
        return OmegaFraction.from_ordinal(value)
    if isinstance(value, bool) or not isinstance(value, numbers.Rational):
        return NotImplemented
    return OmegaFraction._from_polys(_trim([Fraction(value.numerator, value.denominator)]), (_ONE,))


def ordinal_probability(universe: Universe, event: object, initial: object) -> OmegaFraction:
    """``p_μ(E) = 1/ι(ℍ⁺_μ(E))``, or zero when no stage guarantees ``event``.

    It is one over the number of evaluations that surely produce ``event``,
    not counting the start. For a law spread over all of an attractor's part
    of the event, it never exceeds the attractor's in-universe frequency of
    the event (the Kac bridge, SPEC Theorem 6.2). A single start can: one
    whose next step is always in the event has ordinal probability 1.
    """
    guarantee = return_hyperprobability(universe, event, initial)
    if guarantee is INFINITY:
        return OmegaFraction()
    return OmegaFraction.from_ordinal(guarantee).reciprocal()
