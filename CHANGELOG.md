# Changelog

## 0.1.0 - Unreleased

Initial release:

- Contained universes: finite states, exact rational step kernels, and limit rules at
  finitely many levels. Runs are indexed by the ordinals below `omega**omega`.
- Attractors at every level, with their period, kind, stationary law and exit law.
  Strange loops are attractors whose limit is not their own stationary law.
- In-universe probability: the law at a stage, the chance of reaching an event, and
  the long-run frequency of an event.
- Transfinite probability: the law of the first stage at which an event happens, whose
  jumps at limit stages are traced to strange loops.
- Hyperprobability: the least stage by which an event is certain, computed exactly by a
  level recursion. Return guarantees are also available.
- Ordinal probability in the ordered field `Q(omega)`, and the Kac bridge to
  in-universe frequency.
- Subdivided runs (`refine`) and self-simulation (`self_simulation`).
- SPEC.md, with definitions, proofs and a status label on every claim. The tests check
  the results against brute-force search, Kac's lemma and exact identities.
