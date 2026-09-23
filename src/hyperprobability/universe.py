"""Contained universes and their transfinite runs.

A :class:`Universe` is a finite set of states with an exact rational Markov
kernel whose steps never leave the set; that is what *contained* means here.
A run is indexed by the ordinals below ``ω**ω``:

* a successor stage ``α + 1`` takes one step of the kernel from stage ``α``;
* a limit stage ``γ + ω**k`` looks at the run's tail along the stages
  ``γ + ω**(k-1)·n``. Almost surely that tail settles in an *attractor*, a
  closed class of the level ``k - 1`` kernel, and the level ``k`` limit rule
  says where the run goes next.

The default rule *collapses* an attractor to its stationary law, the ordinary
in-universe probability of its states. Any other rule makes the attractor a
*strange loop*: the limit of the loop sends the run somewhere else, possibly
back to where the loop began. Either way the state at a limit depends only on
the run before it, so every run is causal.

Stage ``α = ω**d·c_d + ... + ω·c_1 + c_0`` has kernel
``K_d**c_d ... K_1**c_1 K_0**c_0``, applied left to right. ``K_0`` is the step
kernel and ``K_k`` sends a state to the law at stage ``ω**k`` of a run that
starts there. See SPEC.md §1-2.
"""

from __future__ import annotations

import numbers
from collections.abc import Callable, Hashable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Any, Optional, Union

from ordinatics.ordinals import Ordinal

from hyperprobability._exact import Row, absorption, closed_classes, period, push_power, stationary

# A state is any hashable value. A target, where a limit rule sends an
# attractor, is None (collapse), a state, or a law. A limit rule is None, a
# mapping from attractors (frozensets of states) to targets, or a callable.
State = Hashable
Target = Any
LimitRule = Union[None, Mapping[frozenset, Target], Callable[[frozenset], Target]]

_ONE = Fraction(1)


def exact(value: object) -> Fraction:
    """Return ``value`` as an exact :class:`~fractions.Fraction`.

    Integers, fractions, other :class:`numbers.Rational` values and strings such
    as ``"1/3"`` or ``"0.25"`` are accepted. Floats, complex numbers and booleans
    are rejected: every probability here is an exact rational.
    """
    if isinstance(value, bool):
        raise TypeError("a probability must be an exact rational, not bool")
    if isinstance(value, numbers.Rational):
        return Fraction(value.numerator, value.denominator)
    if isinstance(value, str):
        try:
            return Fraction(value)
        except ValueError:
            raise ValueError(f"{value!r} is not an exact rational") from None
    raise TypeError(f"a probability must be an exact rational, not {type(value).__name__}")


def _count(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise TypeError(f"{name} must be a nonnegative integer")
    result = int(value)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _stage(stage: object) -> Ordinal:
    """Read a stage: an :class:`~ordinatics.ordinals.Ordinal` or a nonnegative integer."""
    if isinstance(stage, Ordinal):
        return stage
    try:
        return Ordinal.from_int(stage)
    except TypeError:
        raise TypeError("a stage is an Ordinal or an exact nonnegative integer") from None


def _last_limit(alpha: Ordinal) -> tuple[int, Ordinal]:
    """Split a limit stage ``λ`` as ``γ + ω**k`` and return ``(k, γ)``."""
    k = next(i for i, c in enumerate(alpha.coefficients) if c)
    values = list(alpha.coefficients)
    values[k] -= 1
    return k, Ordinal(tuple(values))


@dataclass(frozen=True, repr=False)
class Attractor:
    """A closed class of the level-``level`` kernel: where runs of that level settle.

    ``law`` is the class's stationary law, the in-universe probability of each
    of its states over a long run. ``exit`` is the law of the next limit stage
    for a run whose tail settles here. The attractor *collapses* when ``exit``
    is ``law``; otherwise it is a strange loop.
    """

    level: int
    states: tuple[State, ...]
    period: int
    law: Mapping[State, Fraction]
    exit: Mapping[State, Fraction]
    deterministic: bool

    @property
    def members(self) -> frozenset:
        """The states of the attractor, as the key a limit rule uses for it."""
        return frozenset(self.states)

    @property
    def collapses(self) -> bool:
        """Whether the limit of this attractor resolves to its own in-universe law."""
        return dict(self.exit) == dict(self.law)

    @property
    def kind(self) -> str:
        """``"fixed_point"``, ``"periodic_cycle"`` or ``"recurrent_class"``.

        The first two are the attractor kinds of :mod:`ordinatics.dynamics`; the
        third is a class whose runs branch.
        """
        if len(self.states) == 1:
            return "fixed_point"
        return "periodic_cycle" if self.deterministic else "recurrent_class"

    def __contains__(self, state: object) -> bool:
        try:
            return state in self.members
        except TypeError:
            return False

    def __iter__(self) -> Iterator[State]:
        return iter(self.states)

    def __len__(self) -> int:
        return len(self.states)

    def __hash__(self) -> int:
        return hash((Attractor, self.level, self.states))

    def __repr__(self) -> str:
        return f"Attractor(level={self.level}, states={self.states!r}, kind={self.kind!r})"


@dataclass(frozen=True)
class _Level:
    """One level's kernel with its attractors, absorption and exit laws, by state index."""

    rows: list[Row]
    succ: list[list[int]]
    classes: list[list[int]]
    absorption: list[list[Fraction]]
    stationary: list[Row]
    exits: list[Row]


class _Stop:
    """A limit rule of a stopped universe: stopped states stay, other attractors keep their rule."""

    def __init__(self, stops: frozenset, rule: Optional[Callable[[frozenset], Target]]) -> None:
        self.stops = stops
        self.rule = rule

    def __call__(self, members: frozenset) -> Target:
        if len(members) == 1:
            (state,) = members
            if state in self.stops:
                return state
        return None if self.rule is None else self.rule(members)


class Universe:
    """A contained universe: finitely many states, an exact kernel and limit rules.

    ``kernel`` maps each state to a mapping from successor states to exact
    probabilities; its keys, in order, are the states. Every row must be
    nonnegative, sum to exactly one and step only to states of the universe.

    ``limits[k - 1]`` is the level-``k`` limit rule. At a limit stage
    ``γ + ω**k`` it receives the level ``k - 1`` attractor in which the run's
    tail settled and names the target of the limit. A rule is ``None`` (every
    attractor collapses), a mapping from attractors to targets, keyed by
    frozensets of states, or a callable taking such a frozenset. A target is
    ``None`` (collapse to the attractor's stationary law), a state, or a law.
    Attractors a mapping leaves out collapse, and so does every level past the
    last rule.
    """

    def __init__(
        self,
        kernel: Mapping[State, Mapping[State, object]],
        limits: Sequence[LimitRule] = (),
    ) -> None:
        if not isinstance(kernel, Mapping):
            raise TypeError("kernel must map each state to its row of successor probabilities")
        if not kernel:
            raise ValueError("a universe needs at least one state")
        if isinstance(limits, (Mapping, str)) or not isinstance(limits, Sequence):
            raise TypeError("limits must be a sequence of limit rules, one per level")
        self._states: tuple[State, ...] = tuple(kernel)
        self._index = {state: i for i, state in enumerate(self._states)}
        self._rows = [self._row(state, kernel[state]) for state in self._states]
        rules: list[Optional[Callable[[frozenset], Target]]] = []
        tables: list[Optional[dict]] = []
        for level, rule in enumerate(limits, 1):
            if rule is None:
                rules.append(None)
                tables.append(None)
            elif isinstance(rule, Mapping):
                table = dict(rule)
                for key in table:
                    if not isinstance(key, frozenset):
                        raise TypeError(
                            f"limit rule {level} must key attractors by frozensets of states, "
                            f"not {key!r}"
                        )
                rules.append(table.get)
                tables.append(table)
            elif callable(rule):
                rules.append(rule)
                tables.append(None)
            else:
                raise TypeError(f"limit rule {level} must be None, a mapping or a callable")
        self._rules = tuple(rules)
        self._levels: list[_Level] = []
        self._stopped: dict[frozenset[int], Universe] = {}
        self._memo: dict[Any, Any] = {}
        for level in range(self.depth + 1):
            self._level(level)
        for level, table in enumerate(tables):
            if table is None:
                continue
            known = {self._members(members) for members in self._level(level).classes}
            for key in table:
                if key not in known:
                    names = ", ".join(
                        self._show(members) for members in sorted(known, key=self._order)
                    )
                    raise ValueError(
                        f"limit rule {level + 1} names {self._show(key)}, which is not a "
                        f"level-{level} attractor; the level-{level} attractors are {names}"
                    )

    @classmethod
    def deterministic(
        cls, step: Mapping[State, State], limits: Sequence[LimitRule] = ()
    ) -> Universe:
        """The universe in which each state moves to ``step[state]`` with certainty."""
        if not isinstance(step, Mapping):
            raise TypeError("step must map each state to its successor")
        return cls({state: {successor: _ONE} for state, successor in step.items()}, limits)

    @property
    def states(self) -> tuple[State, ...]:
        """The states, in the order of the kernel's keys."""
        return self._states

    @property
    def depth(self) -> int:
        """The number of limit rules; every level past it collapses."""
        return len(self._rules)

    def step(self, state: State) -> dict[State, Fraction]:
        """The law of the stage after ``state``: the kernel's row."""
        if not self._is_state(state):
            raise ValueError(f"{state!r} is not a state of this universe")
        return self._named(self._rows[self._index[state]])

    def kernel(self, level: int = 0) -> dict[State, dict[State, Fraction]]:
        """The level-``level`` kernel ``K_level``: each state's law ``ω**level`` stages later."""
        rows = self._level(_count(level, "level")).rows
        return {state: self._named(rows[i]) for i, state in enumerate(self._states)}

    def attractors(self, level: int = 0) -> tuple[Attractor, ...]:
        """The closed classes of ``K_level``, ordered by their first state."""
        level = _count(level, "level")
        data = self._level(level)
        return tuple(
            Attractor(
                level=level,
                states=tuple(self._states[i] for i in members),
                period=period(members, data.succ),
                law=MappingProxyType(self._named(law)),
                exit=MappingProxyType(self._named(exit_law)),
                deterministic=all(len(data.rows[i]) == 1 for i in members),
            )
            for members, law, exit_law in zip(data.classes, data.stationary, data.exits)
        )

    def law(self, initial: object, stage: object = 0) -> dict[State, Fraction]:
        """The law of the run at ``stage`` (an Ordinal or an integer).

        ``initial`` is a state or a law on states. This is transfinite
        probability: by SPEC Theorem 2.3, ``law(law(μ, α), β) == law(μ, α + β)``.
        """
        return self._named(self._push(self._vector(initial), _stage(stage)))

    def collapse(self, initial: object, stage: object) -> dict[Attractor, Fraction]:
        """How the law at a limit stage collapses out through attractors.

        At ``λ = γ + ω**k`` the run's tail has almost surely settled in a
        level ``k - 1`` attractor. This returns each attractor that can hold
        the tail with the probability that it does. The law at ``λ`` is the
        mixture of their ``exit`` laws with these weights (SPEC §2).
        """
        alpha = _stage(stage)
        if not alpha.is_limit:
            raise ValueError(f"{alpha} is not a limit stage")
        return self._collapse(self._vector(initial), alpha)

    def in_uni(self, initial: object) -> dict[State, Fraction]:
        """The in-universe law: the long-run average of the run's step laws.

        Each attractor contributes its stationary law, weighted by the
        probability that the run settles there. This is ordinary probability,
        the law every default limit rule collapses to.
        """
        return self._named(self._in_uni(self._vector(initial)))

    def stop_at(self, event: object) -> Universe:
        """The universe in which every state of ``event`` is absorbing, at every level.

        Its runs agree with this universe's runs until they first meet the
        event, and stay where they met it. Every attractor of the stopped
        universe that misses the event is an attractor of this universe and
        keeps its limit rule (SPEC Lemma 4.2).
        """
        return self._stop_at(self._event(event))

    def __contains__(self, state: object) -> bool:
        return self._is_state(state)

    def __len__(self) -> int:
        return len(self._states)

    def __repr__(self) -> str:
        return f"Universe(states={self._states!r}, depth={self.depth})"

    # Internals shared with the rest of the package.

    def _stop_at(self, marked: frozenset[int]) -> Universe:
        stopped = self._stopped.get(marked)
        if stopped is None:
            stops = frozenset(self._states[i] for i in marked)
            kernel = {
                state: {state: _ONE} if i in marked else self._named(self._rows[i])
                for i, state in enumerate(self._states)
            }
            stopped = Universe(kernel, [_Stop(stops, rule) for rule in self._rules])
            self._stopped[marked] = stopped
        return stopped

    def _level(self, level: int) -> _Level:
        # Past the last rule every level collapses, and K_k = K_(depth+1) for
        # all k > depth: the Cesaro limit of an idempotent kernel is itself.
        level = min(level, self.depth + 1)
        while len(self._levels) <= level:
            k = len(self._levels)
            rows = self._rows if k == 0 else self._limit_rows(self._levels[-1])
            self._levels.append(self._build(k, rows))
        return self._levels[level]

    def _build(self, k: int, rows: list[Row]) -> _Level:
        succ = [sorted(row) for row in rows]
        classes = closed_classes(len(rows), succ)
        laws = [stationary(rows, members) for members in classes]
        rule = self._rules[k] if k < self.depth else None
        exits = []
        for members, law in zip(classes, laws):
            target = None if rule is None else rule(self._members(members))
            exits.append(law if target is None else self._target(target, k + 1, members))
        return _Level(rows, succ, classes, absorption(rows, classes), laws, exits)

    @staticmethod
    def _limit_rows(level: _Level) -> list[Row]:
        rows = []
        for weights in level.absorption:
            row: Row = {}
            for a, weight in enumerate(weights):
                if weight:
                    for t, p in level.exits[a].items():
                        row[t] = row.get(t, 0) + weight * p
            rows.append({t: p for t, p in row.items() if p})
        return rows

    def _target(self, target: Target, level: int, members: Sequence[int]) -> Row:
        if self._is_state(target):
            return {self._index[target]: _ONE}
        if isinstance(target, Mapping):
            return self._law_row(target, f"the level-{level} limit of {self._show(members)}")
        raise TypeError(
            f"limit rule {level} sends {self._show(members)} to {target!r}, which is neither "
            "a state of this universe nor a law on its states"
        )

    def _push(self, vector: Row, alpha: Ordinal) -> Row:
        for k in range(len(alpha.coefficients) - 1, -1, -1):
            count = alpha.coefficients[k]
            if count:
                vector = push_power(vector, self._level(k).rows, count)
        return vector

    def _collapse(self, vector: Row, alpha: Ordinal) -> dict[Attractor, Fraction]:
        k, before = _last_limit(alpha)
        vector = self._push(vector, before)
        weights = self._level(k - 1).absorption
        result = {}
        for a, attractor in enumerate(self.attractors(k - 1)):
            weight = sum((p * weights[q][a] for q, p in vector.items()), Fraction(0))
            if weight:
                result[attractor] = weight
        return result

    def _in_uni(self, vector: Row) -> Row:
        data = self._level(0)
        result: Row = {}
        for q, p in vector.items():
            for a, weight in enumerate(data.absorption[q]):
                if weight:
                    for t, s in data.stationary[a].items():
                        result[t] = result.get(t, 0) + p * weight * s
        return {t: v for t, v in result.items() if v}

    def _row(self, state: State, row: object) -> Row:
        if not isinstance(row, Mapping):
            raise TypeError(f"the row of {state!r} must map successor states to probabilities")
        return self._law_row(row, f"the row of {state!r}")

    def _law_row(self, law: Mapping, what: str) -> Row:
        result: Row = {}
        total = Fraction(0)
        for state, value in law.items():
            if not self._is_state(state):
                raise ValueError(
                    f"{what} puts mass on {state!r}, which is not a state of this universe"
                )
            p = exact(value)
            if p < 0:
                raise ValueError(f"{what} gives {state!r} the negative probability {p}")
            total += p
            if p:
                result[self._index[state]] = p
        if total != 1:
            raise ValueError(f"{what} sums to {total}, not 1")
        return result

    def _event(self, event: object) -> frozenset[int]:
        if self._is_state(event):
            return frozenset((self._index[event],))
        if callable(event):
            return frozenset(i for i, state in enumerate(self._states) if event(state))
        try:
            members = iter(event)  # type: ignore[call-overload]
        except TypeError:
            raise TypeError("an event is a state, an iterable of states or a predicate") from None
        result = set()
        for state in members:
            if not self._is_state(state):
                raise ValueError(f"{state!r} is not a state of this universe")
            result.add(self._index[state])
        return frozenset(result)

    def _vector(self, initial: object) -> Row:
        if self._is_state(initial):
            return {self._index[initial]: _ONE}
        if isinstance(initial, Mapping):
            return self._law_row(initial, "the initial law")
        raise ValueError(f"{initial!r} is not a state of this universe or a law on its states")

    def _is_state(self, value: object) -> bool:
        try:
            return value in self._index
        except TypeError:
            return False

    def _named(self, row: Row) -> dict[State, Fraction]:
        return {self._states[i]: row[i] for i in sorted(row)}

    def _members(self, indices: Sequence[int]) -> frozenset:
        return frozenset(self._states[i] for i in indices)

    def _order(self, members: frozenset) -> list[int]:
        return sorted(self._index[state] for state in members if self._is_state(state))

    def _show(self, members: object) -> str:
        if isinstance(members, frozenset):
            known = [s for s in members if self._is_state(s)]
            states = sorted(known, key=self._index.__getitem__) + [
                s for s in members if not self._is_state(s)
            ]
        else:
            states = [self._states[i] for i in members]  # type: ignore[union-attr]
        return "{" + ", ".join(repr(state) for state in states) + "}"
