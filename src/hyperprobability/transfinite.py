"""Transfinite probability: the law of the first stage at which an event happens.

For an event ``E`` and an initial law ``μ``, ``τ_E`` is the first ordinal
stage at which the run is in ``E``. Its law lives on the ordinals below
``ω**ω`` and one more point, *never*. It grows by ordinary steps at successor
stages. At a limit stage it jumps only when the run's tail settles in an
attractor whose limit rule carries it into ``E``: a strange loop. Default
rules collapse each attractor to its own in-universe law, so they never jump
(SPEC Theorem 5.2). The law collapses out through these attractors, and
``ℍ_μ(E)`` is the first stage at which it is complete.
"""

from __future__ import annotations

from fractions import Fraction

from ordinatics.ordinals import Ordinal

from hyperprobability._exact import Row
from hyperprobability.hyper import Guarantee, _guarantees
from hyperprobability.universe import Attractor, Universe, _last_limit, _stage


class TransfiniteLaw:
    """The law of ``τ_E``, the first stage at which the run from ``initial`` is in ``event``."""

    def __init__(self, universe: Universe, event: object, initial: object) -> None:
        self._universe = universe
        self._marked = universe._event(event)
        self._stopped = universe._stop_at(self._marked)
        self._initial = universe._vector(initial)

    @property
    def universe(self) -> Universe:
        return self._universe

    @property
    def event(self) -> frozenset:
        """The states of the event."""
        return self._universe._members(sorted(self._marked))

    @property
    def initial(self) -> dict:
        """The initial law."""
        return self._universe._named(self._initial)

    def at_most(self, stage: object) -> Fraction:
        """``P(τ_E ≤ stage)``."""
        return self._mass(self._stopped._push(self._initial, _stage(stage)))

    def before(self, stage: object) -> Fraction:
        """``P(τ_E < stage)``.

        At a limit ``λ = γ + ω**k`` this is the probability that the level
        ``k - 1`` run from stage ``γ`` meets the event (SPEC Proposition 5.1).
        """
        alpha = _stage(stage)
        if not alpha:
            return Fraction(0)
        if alpha.is_successor:
            return self.at_most(alpha.predecessor)
        k, gamma = _last_limit(alpha)
        vector = self._stopped._push(self._initial, gamma)
        data = self._stopped._level(k - 1)
        caught = [a for a, members in enumerate(data.classes) if members[0] in self._marked]
        return sum(
            (p * data.absorption[q][a] for q, p in vector.items() for a in caught), Fraction(0)
        )

    def exactly(self, stage: object) -> Fraction:
        """``P(τ_E = stage)``."""
        return self.at_most(stage) - self.before(stage)

    def never(self) -> Fraction:
        """The probability that the run is in the event at no stage below ``ω**ω``."""
        return 1 - self.before(Ordinal.omega_power(self._universe.depth + 1))

    def guarantee(self) -> Guarantee:
        """``ℍ_μ(E)``, the least stage with ``P(τ_E ≤ stage) = 1``."""
        values = _guarantees(self._universe, self._marked)
        return max(values[q] for q in self._initial)

    def loops_into(self, stage: object) -> dict[Attractor, Fraction]:
        """The strange loops through which the law jumps into the event at ``stage``.

        For a limit ``stage = γ + ω**k`` these are the level ``k - 1``
        attractors of the universe that miss the event, each with the
        probability that the run's tail settles there and its limit lands in
        the event. They sum to ``exactly(stage)`` (SPEC Theorem 5.2). Zero and
        successor stages have none.
        """
        alpha = _stage(stage)
        if not alpha.is_limit:
            return {}
        k, _ = _last_limit(alpha)
        originals = {a.members: a for a in self._universe.attractors(k - 1)}
        event = self.event
        result = {}
        for attractor, weight in self._stopped._collapse(self._initial, alpha).items():
            if attractor.members & event:
                continue
            into = sum((p for state, p in attractor.exit.items() if state in event), Fraction(0))
            if into:
                result[originals[attractor.members]] = weight * into
        return result

    def _mass(self, vector: Row) -> Fraction:
        return sum((p for q, p in vector.items() if q in self._marked), Fraction(0))

    def __repr__(self) -> str:
        return f"TransfiniteLaw(event={set(self.event)!r}, initial={self.initial!r})"
