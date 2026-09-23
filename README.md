# hyperprobability

Ordinary probability says how likely an event is inside a universe.
**Hyperprobability** counts how many evaluations of a simulation of that universe
are guaranteed to produce the event. The count is a transfinite ordinal.

Take a run that swaps between two states forever, where the limit of the swapping is a
way out. No finite stage of the run reaches the exit, and the run spends no time there,
so every in-universe probability of the exit is 0. Yet the exit is certain by stage ω,
the first stage after all the finite ones. Its hyperprobability is ω, and its ordinal
probability is the positive infinitesimal 1/ω.

Everything is exact. There are finitely many states, every probability is a rational
number, and the stages are the ordinals below ω^ω from
[ordinatics](https://github.com/TimeLordRaps/ordinatics).

## Install

```bash
pip install git+https://github.com/TimeLordRaps/hyperprobability
```

Python 3.10 or later. The only dependency is ordinatics 0.3, which pip installs from
PyPI.

## Example

```python
from hyperprobability import (
    OMEGA,
    TransfiniteLaw,
    Universe,
    frequency,
    hyperprobability,
    ordinal_probability,
    reach_probability,
)

# a and b swap forever. The limit of the swapping sends the run to "out".
loop = Universe.deterministic(
    {"a": "b", "b": "a", "out": "out"},
    [{frozenset({"a", "b"}): "out"}],
)

# No finite stage reaches "out", and the run spends no time there,
print(reach_probability(loop, "out", "a"))  # 0
print(frequency(loop, "out", "a"))  # 0
# but "out" is certain by stage ω.
print(hyperprobability(loop, "out", "a"))  # ω
print(ordinal_probability(loop, "out", "a"))  # 1/ω

# The whole jump at ω comes through one strange loop.
for attractor, p in TransfiniteLaw(loop, "out", "a").loops_into(OMEGA).items():
    print(attractor.states, p)  # ('a', 'b') 1
```

A fair coin shows heads eventually, but no finite number of tosses is sure to:

```python
toss = {"H": "1/2", "T": "1/2"}
coin = Universe({"H": toss, "T": toss})

print(frequency(coin, "H", "T"))  # 1/2
print(hyperprobability(coin, "H", "T"))  # ω
print(ordinal_probability(coin, "H", "H"))  # 1/ω

first_heads = TransfiniteLaw(coin, "H", "T")
print(*(first_heads.exactly(n) for n in range(1, 5)))  # 1/2 1/4 1/8 1/16
print(first_heads.at_most(OMEGA))  # 1
```

## Concepts

- **Universe.** A finite set of states, an exact rational step kernel, and limit rules
  at finitely many levels (`Universe`, `Universe.deterministic`). A run is indexed by
  the ordinals below ω^ω. A successor stage takes one step. At a limit stage the run's
  tail has settled in an *attractor*, a closed class of the level below, and that
  level's rule says where the run goes next. By default an attractor *collapses* to its
  stationary law. Any other target makes it a **strange loop** (`strange_loops`).
- **In-universe probability.** The chance of the event at a stage (`probability`), the
  chance of reaching it at a finite stage (`reach_probability`), and the long-run
  fraction of time spent in it (`frequency`).
- **Transfinite probability.** `TransfiniteLaw` is the law of the first stage at which
  the run is in the event. It grows by ordinary steps, and jumps at a limit stage only
  through strange loops (`loops_into`).
- **Hyperprobability.** `hyperprobability(universe, event, start)` is the least stage by
  which the event is certain, or `INFINITY` if the run can miss it.
  `return_hyperprobability` does not count the start. Unless it is `INFINITY`, a
  guarantee is at most ω^(d+1), where d is the number of levels with rules, and each of
  its coefficients is less than the number of states. Hyperprobability depends only on
  which transitions are possible, not on their weights.
- **Ordinal probability.** `ordinal_probability` is 1/ℍ⁺ in ℚ(ω), the ordered field of
  rational functions in ω (`OmegaFraction`). It is 1/N for a guarantee of N steps, an
  infinitesimal for an infinite guarantee, and 0 when the event can be missed. It is
  monotone, but it is not a measure.
- **The Kac bridge.** `kac(universe, event)` looks at each level-0 attractor that meets
  the event. There the in-universe frequency of the event is at least 1/ℍ⁺, where ℍ⁺
  is the worst guarantee of a return to the event. Equality holds exactly when every
  return takes the same number of steps.
- **Simulations of simulations.** `refine(ℍ, c)` = c·ℍ is the guarantee when each step
  costs c evaluations of a finer simulation. `self_simulation` finds the least guarantee
  that pays for running itself: h·ω when the cost is h and then the run again.

## Examples

| Universe | Event | Start | ℍ | ordinal probability | chance of reaching it at a finite stage | frequency |
|---|---|---|---|---|---|---|
| clock | tock | tick | 1 | 1 | 1 | 1/2 |
| coin | H | T | ω | 1/ω | 1 | 1/2 |
| loop | out | a | ω | 1/ω | 0 | 0 |
| chain | out | a | ω + 3 | 1/(ω + 3) | 0 | 0 |
| swap | out | a | ω² | 1/ω² | 0 | 0 |
| dead end | out | a | ∞ | 0 | 0 | 0 |

`examples/coin_and_clock.py` and `examples/strange_loop.py` print these and more. SPEC
Example 3.4 defines the universes.

## Specification

[SPEC.md](SPEC.md) gives the definitions, proves what it can, and labels every claim as
DEFINITION, PROVED, KNOWN, CHECKED, INTERPRETATION, TRANSPORT, CONJECTURE or OPEN. Its
appendix maps each result to the tests that check it. For example, the level recursion
that computes hyperprobability is checked against brute-force search over stages in
2000 random universes. No proof is machine-checked yet.

## Family

hyperprobability builds on [ordinatics](https://github.com/TimeLordRaps/ordinatics).
It sits beside [hypermath](https://github.com/TimeLordRaps/hypermath),
[hyperlogic](https://github.com/TimeLordRaps/hyperlogic),
[hyperethics](https://github.com/TimeLordRaps/hyperethics),
[hyperphysics](https://github.com/TimeLordRaps/hyperphysics) and
[hyperstratum](https://github.com/TimeLordRaps/hyperstratum). SPEC §9 describes how
they relate.

## Development

```bash
python -m pip install -e ".[dev]"
```

```bash
python -m pytest
```

```bash
ruff check .
```

```bash
ruff format --check .
```

## Citation

See [CITATION.cff](CITATION.cff).

## License

Apache-2.0. See [LICENSE](LICENSE).
