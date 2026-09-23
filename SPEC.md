# Hyperprobability: specification

Version 0.1.0 (unreleased). This document defines hyperprobability, transfinite
probability and ordinal probability for contained universes. It proves what it can and
labels every claim. The Python package in `src/hyperprobability` implements these
definitions exactly, and the appendix lists the tests behind each result.

## 0. Reading and conventions

### 0.1 What the request asks for

The repository started from Tyler Roost's description of the idea. Each phrase below is
read one way (INTERPRETATION), and the section that makes the reading precise is named.

| Phrase | Reading here | Where |
|---|---|---|
| probability is in-uni | ordinary probability inside a contained universe: the chance of an event at a stage, and its long-run frequency | Def 2.6, `inuni.py` |
| minimum guaranteed simulation evaluations that result in that event, measured in transfinite ordinals | ℍ_μ(E), the least ordinal stage α with P(τ_E ≤ α) = 1 | Def 4.1, `hyper.py` |
| transfinite probability | the law of τ_E, the first stage at which the run is in E, over the ordinal stages below ω^ω | Def 4.1, §5, `transfinite.py` |
| collapses out through causal strange loop attractors | at a limit stage the law is a mixture over the attractors the run can settle in, and the law of τ_E jumps only through strange loops | Prop 2.5, Thm 5.2, `loops.py` |
| ordinal probability, a higher-order concept that hyperprobability explains | p = 1/ι(ℍ⁺) in the ordered field ℚ(ω), defined from hyperprobability | §7, `ordinal.py` |
| possibility dynamics | ℍ depends only on which transitions are possible, not on their weights | Prop 4.9 |
| contained, constraint, causally representable universe | finite states, exact rational laws that never leave them, and limit stages that depend only on the past | Def 1.1, Remark 1.3 |

A *simulation evaluation* is one stage of a run: a step of the kernel at a successor
stage, or one application of a limit rule at a limit stage. It is not hypermath's
simulation relation (§9).

### 0.2 Claim labels

| Label | Meaning |
|---|---|
| DEFINITION | a stipulation |
| PROVED | proved in this document; "(sketch)" marks a proof whose routine steps are left out |
| KNOWN | a published result, cited |
| CHECKED | verified by the test suite on exact instances, not proved |
| INTERPRETATION | a reading of a phrase or a construction, not a mathematical claim |
| TRANSPORT | a form borrowed from another field, declared in hyperphysics' format |
| CONJECTURE | believed and unproved |
| OPEN | a question without a proposed answer |

No proof here is machine-checked.

### 0.3 Notation

- S is a finite nonempty set of states and n = |S|.
- A *law* is a probability vector on S with rational entries. A *kernel* is a stochastic
  matrix on S. Laws are row vectors: (μK)(t) = Σ_s μ(s)K(s, t). δ_q is the point mass at
  q, supp μ is the support of μ, and μ(E) = Σ_{s∈E} μ(s).
- Ordinals below ω^ω are written in Cantor normal form, α = ω^m·c_m + ⋯ + ω·c_1 + c_0
  with c_i ∈ ℕ. The code prints this as `ω^2*3 + ω*2 + 1`.
- α + β and α·β are ordinal addition and multiplication. They are not commutative:
  1 + ω = ω < ω + 1, and 2·ω = ω < ω·2 = ω + ω.
- α ⊕ β and α ⊗ β are Hessenberg's natural sum and product, which add and multiply
  Cantor normal forms like polynomials in ω.
- ∞ is a symbol above every ordinal, and ∞ + α = α + ∞ = ∞. In code it is `INFINITY`.

## 1. Contained universes

**Definition 1.1 (universe).** DEFINITION. A universe U = (S, K_0, L_1, …, L_d) has:

- a finite nonempty set S of states;
- a step kernel K_0 with exact rational entries;
- a depth d ≥ 0 and, for 1 ≤ k ≤ d, a limit rule L_k.

L_k assigns a *target* to each level-(k−1) attractor (Def 1.2). A target is either
*collapse* or a law on S; a state q stands for δ_q. For k > d, L_k collapses every
attractor.

In code this is `Universe(kernel, limits)`, where `limits[k-1]` is L_k. A rule is
`None`, a mapping keyed by frozensets of states, or a callable. A mapping collapses
every attractor it leaves out. `Universe.deterministic(step, limits)` builds a universe
whose rows are point masses.

**Definition 1.2 (levels).** DEFINITION. By recursion on k ≥ 0:

- 𝒜_k is the set of closed communicating classes of K_k, the *level-k attractors*.
- For A ∈ 𝒜_k, π_A is the unique stationary law of K_k carried by A.
- h_k(q, A) is the probability that the K_k-chain started at q is eventually in A.
- The *exit law* is e_A = π_A if L_{k+1}(A) collapses, and e_A = L_{k+1}(A) otherwise.
- K_{k+1}(q, ·) = Σ_{A∈𝒜_k} h_k(q, A)·e_A.

Code: `Universe.kernel(k)` and `Universe.attractors(k)`, with `Attractor.law` = π_A and
`Attractor.exit` = e_A.

**Remark 1.3 (contained, causal, representable).** INTERPRETATION.

- *Contained*: every row of every K_k and every target is a law on S, so no run ever
  leaves S. The constructor rejects rows and targets that put mass outside S or do not
  sum to exactly 1. Containment is the constraint.
- *Causal*: the state at a limit stage depends only on the run below it, through the
  attractor its tail settled in (Remark 2.8). Nothing at stage α depends on a later
  stage.
- *Representable*: S is finite and every probability is an exact rational, so the law
  at every stage below ω^ω is computed exactly (Def 2.1), without floating point or
  sampling. `exact` rejects floats.

## 2. Transfinite runs

**Definition 2.1 (law at a stage).** DEFINITION. For α = ω^m·c_m + ⋯ + ω·c_1 + c_0, let

    K_α = K_m^{c_m} ⋯ K_1^{c_1} K_0^{c_0},

with the factors applied left to right, and K_0 = I at stage 0. The law at stage α of
the run from μ is μK_α. Code: `Universe.law(μ, α)`.

**Lemma 2.2.** PROVED.
(a) Every K_k is stochastic.
(b) K_j K_k = K_k for all j < k.

*Proof.*
(a) A finite chain is absorbed into its closed classes with probability 1, so
Σ_A h_k(q, A) = 1 [Kemeny–Snell 1960; Norris 1997]. Each e_A is a law.
(b) Each h_k(·, A) is K_k-harmonic: Σ_r K_k(q, r)·h_k(r, A) = h_k(q, A). For transient q
this is first-step analysis; for q in a closed class B, both sides equal [B = A]. So
K_k K_{k+1} = K_{k+1}. By induction on k − j,
K_j K_k = K_j (K_{k−1} K_k) = (K_j K_{k−1}) K_k = K_{k−1} K_k = K_k. ∎

**Theorem 2.3 (composition).** PROVED. For α, β < ω^ω, K_{α+β} = K_α K_β. So the law at
β of a run started from the law at α is the law at α + β.

*Proof.* The case β = 0 is trivial. Otherwise let ω^m be the leading power of β, with
coefficient b_m ≥ 1. Ordinal addition drops the terms of α below ω^m and adds the
coefficients at ω^m:

    α + β = ω^M·a_M + ⋯ + ω^m·(a_m + b_m) + ω^{m−1}·b_{m−1} + ⋯ + b_0.

In K_α K_β, the block K_{m−1}^{a_{m−1}} ⋯ K_0^{a_0} is followed by K_m, which absorbs it
by Lemma 2.2(b). What remains is K_{α+β}. ∎

**Proposition 2.4 (stabilization).** PROVED. K_{d+1} is the Cesàro limit of K_d and is
idempotent, and K_k = K_{d+1} for every k > d.

*Proof.* L_{d+1} collapses every attractor, so K_{d+1}(q, ·) = Σ_A h_d(q, A)·π_A. That
is the Cesàro limit lim_N (1/N) Σ_{i<N} K_d^i, which is idempotent [Kemeny–Snell 1960].
The Cesàro limit of an idempotent kernel is the kernel itself, so K_{d+2} = K_{d+1}. The
rest follows by induction. ∎

The kernels stabilize, but the laws need not (CHECKED). Take depth 1 and S = {a, b, x}:

- K_0: a → a, b → x, x → x;
- L_1: {a} ↦ b and {x} ↦ a.

Then K_1 sends a → b, b → a and x → a. Its one attractor {a, b} collapses, so K_2 is
{a: ½, b: ½} from every state. From a:

| stage | law |
|---|---|
| ω² | {a: ½, b: ½} |
| ω² + 1 | {a: ½, x: ½} |
| ω² + ω | {a: ½, b: ½} |

What does stop moving at ω^{d+1} is the probability of having met an event
(Lemma 4.2(d)).

**Proposition 2.5 (limit laws mix exits).** PROVED. Let λ = γ + ω^k with k ≥ 1, and let
ν_m be the law at γ + ω^{k−1}·m. Let T be the set of states outside every level-(k−1)
attractor.

(a) The law at λ is Σ_{A∈𝒜_{k−1}} w_A·e_A. Here w_A = Σ_q ν_0(q)·h_{k−1}(q, A) is the
probability that the level-(k−1) run from stage γ settles in A. Code:
`Universe.collapse(μ, λ)` returns the weights w_A.

(b) For every m ≥ 0: ν_m(A) ≤ w_A ≤ ν_m(A) + ν_m(T).

*Proof.* (a) The law at λ is ν_0 K_k (Theorem 2.3); expand K_k. (b) Harmonicity of
h_{k−1}(·, A) gives w_A = Σ_q ν_m(q)·h_{k−1}(q, A). A state in A contributes 1, a state
in another attractor contributes 0, and a state in T contributes at most 1. ∎

**Definition 2.6 (in-universe law).** DEFINITION. The in-universe law is
in_uni(μ) = Σ_{A∈𝒜_0} (Σ_q μ(q)·h_0(q, A))·π_A. The in-universe frequency of E is
in_uni(μ)(E). Code: `Universe.in_uni`, `frequency`, `probability`, `reach_probability`.

**Proposition 2.7.** KNOWN (a, b); PROVED (c, d).

(a) in_uni(μ) is the Cesàro limit lim_N (1/N) Σ_{i<N} μK_0^i.
(b) Almost surely the run settles in some A ∈ 𝒜_0, and the fraction of stages i < N at
which it is in E then tends to π_A(E) (the ergodic theorem for finite irreducible
chains [Norris 1997]). So in_uni(μ)(E) is the expected long-run frequency of E.
(c) If L_1 collapses every attractor, the law at ω is in_uni(μ).
(d) in_uni(μ)·K_0 = in_uni(μ).

Strange loops do not change in_uni, which depends only on K_0.

*Proof of (c) and (d).* (c) K_1 is then the Cesàro limit of K_0. (d) Each π_A is
K_0-invariant. ∎

**Remark 2.8 (pathwise runs).** PROVED (sketch). A random run (X_α)_{α<ω^ω} is built
along the countable well-order ω^ω:

- X_0 ~ μ;
- X_{α+1} ~ K_0(X_α, ·);
- at λ = γ + ω^k, the tail (X_{γ+ω^{k−1}·m})_{m∈ℕ} is a K_{k−1}-chain (by induction on k),
  so it is eventually in one attractor A almost surely. Then X_λ ~ e_A, drawn
  independently of everything else given A.

By Ionescu-Tulcea's theorem, transfinite induction gives a process whose law at α is
μK_α. A collapsing attractor sends its tail to π_A, forgetting where in A the tail was.
Only laws are used below; this construction is the intended reading of "the run".

## 3. Attractors and strange loops

**Definition 3.1 (attractor).** DEFINITION. A level-k attractor A ∈ 𝒜_k carries:

- its *period*: the gcd of the cycle lengths of the support graph of K_k on A;
- its *kind*: a fixed point if |A| = 1, a periodic cycle if every row of K_k on A is a
  point mass, and a recurrent class otherwise;
- its stationary law π_A and its exit law e_A.

The first two kinds are the attractor kinds of `ordinatics.dynamics`. Code: `Attractor`.

**Definition 3.2 (strange loop).** DEFINITION. An attractor A is a *strange loop* if
e_A ≠ π_A: the limit of a run that settles in A does not resolve to A's in-universe law.
Otherwise A *collapses*. Only attractors of levels k < d can be strange loops. Code:
`strange_loops`, `Attractor.collapses`.

**Remark 3.3 (Hofstadter).** TRANSPORT, declared in hyperphysics' five-item format:

1. *Source, by name*: Hofstadter's strange loop [Hofstadter 1979, 2007].
2. *Target form*: a closed class of K_k whose level-(k+1) limit sends the run somewhere
   other than its own stationary law, possibly back into the class.
3. *What transports*: going around a loop at one level lands you somewhere unexpected
   at the level above. The limit of the loop makes a move that no pass around the loop
   makes.
4. *What does not transport*: self-reference, symbols and representation, and every
   claim about minds. Nothing in a universe refers to itself.
5. *Whether the target constrains a parameter*: yes. The target fixes the number of
   levels crossed, exactly one (from k to k + 1), which the source leaves informal.

The "strange attractor" of chaotic dynamics is unrelated.

**Example 3.4 (the running examples).** DEFINITION. Tests and examples use these
universes (`tests/support.py`):

| Universe | Step kernel K_0 | Limit rules | Strange loops |
|---|---|---|---|
| clock | tick → tock → tick | none | none |
| coin | every state → {H: ½, T: ½} | none | none |
| loop | a → b → a, out → out | L_1: {a, b} ↦ out | {a, b} at level 0 |
| chain | a → b → a, x1 → x2 → x3 → out → out | L_1: {a, b} ↦ x1 | {a, b} at level 0 |
| fork | a → b → a, x → out, y1 → y2 → out, out → out | L_1: {a, b} ↦ {x: ½, y1: ½} | {a, b} at level 0 |
| swap | a → a, b → b, out → out | L_1: {a} ↦ b, {b} ↦ a; L_2: {a, b} ↦ out | {a}, {b} at level 0; {a, b} at level 1 |
| dead end | a → a, out → out | none | none |

In the swap, K_1 sends a → b and b → a, so {a, b} is a level-1 attractor.

## 4. Hyperprobability

**Definition 4.1.** DEFINITION. Let E ⊆ S be an event and μ an initial law.

- The *stopped universe* U^E has step kernel K^E_0(q) = δ_q for q ∈ E and
  K^E_0(q) = K_0(q) otherwise. Its rules are L^E_k({e}) = e for e ∈ E, and
  L^E_k(A) = L_k(A) for every other level-(k−1) attractor A of U^E. Such an A is an
  attractor of U by Lemma 4.2(b) at level k − 1, so L_k(A) is defined. Code:
  `Universe.stop_at`.
- F_μ(α) = (μK^E_α)(E) is the probability that the run has met E by stage α, written
  P_μ(τ_E ≤ α).
- P_μ(τ_E < α) = sup_{β<α} F_μ(β), and 0 at α = 0.
- P_μ(τ_E = α) = F_μ(α) − P_μ(τ_E < α).
- P_μ(never) = 1 − F_μ(ω^{d+1}), the probability that the run is in E at no stage below
  ω^ω (Lemma 4.2(d)). By Theorem 5.2 this also equals 1 − P_μ(τ_E < ω^{d+1}).
- The *hyperprobability* of E is ℍ_μ(E) = min{α < ω^ω : F_μ(α) = 1}, or ∞ if there is
  no such α.

Code: `TransfiniteLaw(U, E, μ)` with `at_most`, `before`, `exactly`, `never` and
`guarantee`; `hyperprobability(U, E, μ)`; `hyperprobabilities(U, E)` for every state
at once. An event is a state, an iterable of states, or a predicate.

Pathwise (INTERPRETATION), τ_E is the first stage at which the run of Remark 2.8 is in
E. The stopped run agrees with the run until τ_E and stays in E afterwards. The
definitions above are the law-level version, and they are what the code computes.

**Lemma 4.2.** PROVED. Every e ∈ E is absorbing at every level: K^E_k(e) = δ_e. For
every k ≥ 0:

(a) if (νK^E_k)(E) = 0 for a law ν, then νK^E_k = νK_k;
(b) every attractor of K^E_k other than the singletons {e} with e ∈ E misses E, is an
attractor of K_k, and has the same rows, stationary law and exit;
(c) F_μ is weakly increasing;
(d) F_μ is constant on [ω^{d+1}, ω^ω).

*Proof.* **Absorption.** By induction on k: {e} is a closed class of K^E_{k−1} with
exit δ_e.

**(a) and (b), by joint induction on k.**

*Level 0.* The hypothesis of (a) forces ν(E) = 0, and off E the two kernels have the
same rows. For (b), a closed class of K^E_0 containing e ∈ E is {e}, because e is
absorbing. Any other class A misses E, so its rows are rows of K_0. Hence A is closed and
irreducible for K_0, with the same stationary law, and L^E_1(A) = L_1(A).

*Level k ≥ 1.* Assume (a) and (b) at level k − 1.

For (a), write K^E_k(q, ·) = Σ_A h^E_{k−1}(q, A)·e^E_A. The hypothesis gives every
singleton {e} ⊆ E zero weight, so the K^E_{k−1}-chain from ν never enters E. Every
state r it can reach therefore has (δ_r K^E_{k−1})(E) = 0, and (a) at k − 1 gives
K^E_{k−1}(r) = K_{k−1}(r). So the two chains from ν agree on everything they reach. They
reach the same attractors with the same absorption probabilities, and by (b) at k − 1
those attractors have the same stationary laws and exits. Hence νK^E_k = νK_k.

For (b), a class containing some e ∈ E is {e}. Let A be any other class. For r ∈ A,
K^E_k(r) is carried by A, which misses E, so K^E_k(r) = K_k(r) by (a). Hence A is a
closed irreducible class of K_k with the same stationary law, and its exit is
L^E_{k+1}(A) = L_{k+1}(A).

**(c).** If α ≤ β, write β = α + δ. By Theorem 2.3,
F_μ(β) = Σ_s (μK^E_α)(s)·(δ_s K^E_δ)(E), and δ_s K^E_δ = δ_s for s ∈ E. So
F_μ(β) ≥ F_μ(α).

**(d).** U^E has depth d. By Lemma 2.2(b) and Proposition 2.4 applied to U^E,
K^E_j K^E_{d+1} = K^E_{d+1} for every j. So f(s) = K^E_{d+1}(s, E) satisfies
K^E_j f = f for every j: f is harmonic at every level. Since f ≥ 0, the set
Z = {s : f(s) = 0} is closed under every K^E_j, and so is E.

The law ν = μK^E_{d+1} is a mixture of stationary laws of attractors of K^E_d. The
singletons {e} lie in E. Every other attractor misses E and is closed, so f = 0 on it.
Hence ν is carried by E ∪ Z.

Every α ≥ ω^{d+1} equals ω^{d+1} + δ for some δ. The kernel K^E_δ keeps E's mass in E
and Z's mass in Z, and Z is disjoint from E, so F_μ(α) = (νK^E_δ)(E) = ν(E) =
F_μ(ω^{d+1}). ∎

**Proposition 4.3 (initial laws).** PROVED. ℍ_μ(E) = max_{q∈supp μ} ℍ_q(E).

*Proof.* F_μ = Σ_q μ(q)·F_q equals 1 exactly when every F_q with μ(q) > 0 does. By
Lemma 4.2(c), {α : F_q(α) = 1} is the final segment starting at ℍ_q(E). ∎

**Theorem 4.4 (the level recursion).** PROVED. Let R_{−1}(q) = 0 for q ∈ E and ∞
otherwise. For 0 ≤ k ≤ d, let Supp_c(q) = supp(δ_q (K^E_k)^c). Call c *good* for q if
R_{k−1}(r) < ω^k for every r ∈ Supp_c(q) \ E. Define R_k(q) as the first case that
applies:

- if some c is good, then R_k(q) = ω^k·c + max{R_{k−1}(r) : r ∈ Supp_c(q) \ E} for the
  least good c, with max ∅ = 0;
- if no attractor of K^E_k that misses E is reachable from q, then R_k(q) = ω^{k+1};
- otherwise R_k(q) = ∞.

Then:

(i) the least good c, when there is one, is at most n − 1;
(ii) R_k(q) is the least α < ω^{k+1} with F_q(α) = 1 if there is one; otherwise it is
ω^{k+1} if P_q(τ_E < ω^{k+1}) = 1; otherwise it is ∞;
(iii) ℍ_q(E) = R_d(q).

This recursion is what the code computes. `_resolve` in `hyper.py` is one step of it.

*Proof.* **(ii), by induction on k.** The case k = −1 holds because stage 0 is certain
exactly when q ∈ E. Let k ≥ 0.

*The least certain stage below ω^{k+1}.* Every α < ω^{k+1} is uniquely ω^k·c + β with
β < ω^k. By Theorem 2.3,

    F_q(ω^k·c + β) = Σ_r (δ_q (K^E_k)^c)(r)·F_r(β),

which is 1 exactly when F_r(β) = 1 for every r ∈ Supp_c(q). This holds for every r ∈ E.
For r ∉ E and β < ω^k, F_r(β) = 1 exactly when r has a certain stage below ω^k and the
least one is at most β. By induction, that is R_{k−1}(r) ≤ β: R_{k−1}(r) is that least
stage when it lies below ω^k, and it is ω^k or ∞ otherwise. So ω^k·c + β is certain
exactly when c is good and β ≥ max{R_{k−1}(r) : r ∈ Supp_c(q) \ E}. The least certain
stage below ω^{k+1} is therefore the first case of the definition.

*The second and third cases.* The stages ω^k·c are cofinal below ω^{k+1}, and F is
increasing. So P_q(τ_E < ω^{k+1}) = lim_c (δ_q (K^E_k)^c)(E), the probability that the
K^E_k-chain from q is absorbed in E (Prop 5.1(a)). This is 1 exactly when no attractor
that misses E is reachable from q: a finite chain is absorbed into its closed classes
almost surely, and each reachable closed class gets positive probability.

**(i).** Call r *bad* if r ∉ E and R_{k−1}(r) ≥ ω^k. Suppose c is good. For each
r ∈ Supp_c(q) \ E, F_r(R_{k−1}(r)) = 1, so F_r(ω^k) = 1 by Lemma 4.2(c), and K^E_k(r) is
carried by E. Since E is absorbing, Supp_{c'}(q) ⊆ E for every c' > c. So every
K^E_k-path from q to a bad state is shorter than c.

A path from q to a bad state b of length at least n − 1 would either repeat a state or
visit all n states. In the first case its cycle can be traversed again and again. In
the second case b's row points back to some state on the path, which closes a cycle
through b. Either way there would be arbitrarily long paths from q to b. So every such
path is shorter than n − 1, and c = n − 1 is good.

**(iii).**
- If R_d(q) < ω^{d+1}, then R_d(q) is the least certain stage, which is ℍ_q(E).
- If R_d(q) = ω^{d+1}, no smaller stage is certain, and F_q(ω^{d+1}) ≥
  P_q(τ_E < ω^{d+1}) = 1.
- If R_d(q) = ∞, every level-d attractor collapses, so the ones that miss E add nothing
  at ω^{d+1} (Theorem 5.2). Hence F_q(ω^{d+1}) = P_q(τ_E < ω^{d+1}) < 1, and F_q stays
  there by Lemma 4.2(d). ∎

Theorem 5.2 and Proposition 5.1(a) do not depend on Theorem 4.4.

In the loop, R_0(a) = ∞ and R_1(a) = ω. At level 0, P_a(τ < ω) = 0 even though
F_a(ω) = 1. The jump at ω is found one level up, as the stage ω·1 + 0. This is why R_k
records P(τ_E < ω^{k+1}) = 1 and not F(ω^{k+1}) = 1.

**Corollary 4.5 (shape).** PROVED. ℍ_q(E) is ∞, or ω^{d+1}, or
ω^d·c_d + ⋯ + ω·c_1 + c_0 with every c_i ≤ n − 1.

*Proof.* By Theorem 4.4, R_k(q) is ∞, or ω^{k+1}, or ω^k·c + R_{k−1}(r) with c ≤ n − 1
and R_{k−1}(r) < ω^k. Induct on k. ∎

In particular a guarantee within finitely many steps is a guarantee within n − 1 steps.
The oracle in `tests/support.py` searches coefficients up to 2^n + 1 and shares no code
with the recursion.

**Proposition 4.6 (monotone in the event).** PROVED. If E ⊆ E', then
ℍ_μ(E') ≤ ℍ_μ(E).

*Proof.* **Claim:** for any universe V and any G ⊆ S,
supp K^{V^G}_k(q) ⊆ supp K^V_k(q) ∪ G for all k and q. We prove it by induction on k.
At k = 0 it is immediate. For the step, write K^{V^G}_{k+1}(q) = Σ_A h(q, A)·e_A.

- Singletons in G have exits in G.
- Any other attractor A of V^G misses G. By Lemma 4.2(b), A is an attractor of V with
  the same exit. A K^{V^G}_k-path from q to A never enters G, since G is absorbing. By
  the induction hypothesis each of its steps is also a step of K^V_k, so A is reachable
  under K^V_k and supp e_A ⊆ supp K^V_{k+1}(q).

Composing, the claim extends to every stage kernel K_α.

**Applying the claim.** Take V = U^E and G = E'. Stopping at E and then at E' ⊇ E is
stopping at E', so V^G = U^{E'}. If F^E_q(α) = 1, then supp(δ_q K^{U^E}_α) ⊆ E. So
supp(δ_q K^{U^{E'}}_α) ⊆ E ∪ E' = E', that is, F^{E'}_q(α) = 1. ∎

**Definition 4.7 (return hyperprobability).** DEFINITION. For β < ω^ω, let
P_μ(τ⁺_E ≤ 1 + β) = F_{μK_0}(β): the probability that the run after its first step has
met E. Then ℍ⁺_μ(E) = min{α ≥ 1 : P_μ(τ⁺_E ≤ α) = 1}, or ∞ if there is no such α. The
start does not count. Every α ≥ 1 is 1 + β for exactly one β. Code:
`return_hyperprobability`.

**Proposition 4.8.** PROVED. ℍ⁺_μ(E) = 1 + max_{r∈supp μK_0} ℍ_r(E), with ordinal
addition, so 1 + ω = ω and 1 + ∞ = ∞.

*Proof.* β ↦ 1 + β is strictly increasing onto [1, ω^ω). Apply Proposition 4.3 to
μK_0. ∎

**Proposition 4.9 (possibility dynamics).** PROVED. ℍ_q(E) depends only on E and on
the support graphs of K_0, …, K_d. Equivalently, it depends only on E, the support
graph of K_0 and the supports of the exits of the strange loops. Reweighting rows and
targets without changing their supports changes no hyperprobability.

*Proof.* The recursion of Theorem 4.4 reads only these things:

- the event E;
- the supports of K^E_0, …, K^E_d;
- their closed classes and reachability, both of which the supports determine.

supp K^E_0 is determined by supp K_0 and E. supp K^E_{k+1}(q) is the union of supp e^E_A
over the attractors A of K^E_k reachable from q. Here e^E_{e} = δ_e, and for any other A,
e^E_A = e_A (Lemma 4.2(b)). Its support is supp K_{k+1}(r) for any r ∈ A, because
h_k(r, A) = 1. ∎

INTERPRETATION: hyperprobability sees which evaluations are possible. In-universe
probability sees how likely they are. This is the sense in which hyperprobability is a
*possibility dynamics*.

## 5. Transfinite probability and the collapse

**Proposition 5.1 (just before a limit).** PROVED. Let λ = γ + ω^k with k ≥ 1. Let
ν_m = μK^E_γ (K^E_{k−1})^m be the stopped law at γ + ω^{k−1}·m, and let T be the set of
states outside every attractor of K^E_{k−1}.

(a) P_μ(τ_E < λ) = Σ_q ν_0(q)·Σ_{e∈E} h^E_{k−1}(q, {e}) = lim_m F_μ(γ + ω^{k−1}·m).
(b) For every m, F_μ(γ + ω^{k−1}·m) ≤ P_μ(τ_E < λ) ≤ F_μ(γ + ω^{k−1}·m) + ν_m(T).
(c) Suppose every row of K^E_{k−1} is a point mass. This holds, for example, when K_0
is deterministic and L_1, …, L_{k−1} send every attractor to a single state. Then the run
meets E within n − 1 steps of level k − 1 or never, so
P_μ(τ_E < λ) = F_μ(γ + ω^{k−1}·m) for every m ≥ n − 1.

Code: `TransfiniteLaw.before` uses (a).

*Proof.*
(a) The stages γ + ω^{k−1}·m are cofinal in λ and F is increasing. The K^E_{k−1}-chain
from ν_0 is in E at time m exactly when it has been absorbed at some e ∈ E by then, and
these probabilities increase to the absorption probabilities.
(b) After time m, only the mass in T can still be absorbed in E, since the other
attractors miss E and are closed.
(c) The first n positions of a deterministic path either meet E, or repeat a state and
so cycle outside E forever. If they are n distinct states outside E, then E is empty. ∎

**Theorem 5.2 (limit jumps come from strange loops).** PROVED. With λ and ν_0 as above,

    P_μ(τ_E = λ) = Σ_A w_A·e_A(E),

summed over the attractors A of K^E_{k−1} that miss E. By Lemma 4.2(b) these are
level-(k−1) attractors of U, with the same exits. Here w_A = Σ_q ν_0(q)·h^E_{k−1}(q, A)
is the probability that the stopped run's tail settles in A. A collapsing A contributes
π_A(E) = 0. So the law of τ_E can jump at a limit stage only through strange loops whose
exits reach E.

The converse of Lemma 4.2(b) fails, so the sum cannot run over every attractor of U that
misses E. Let K_0 be r → e → s → s, let E = {e}, and let L_1({s}) = r. Then {r} is a
level-1 attractor of U, but in U^E the state r steps to e and is absorbed there.

Code: `TransfiniteLaw.loops_into`.

*Proof.* F_μ(λ) = (ν_0 K^E_k)(E) = Σ_{A∈𝒜^E_{k−1}} w_A·e^E_A(E). The singletons in E
contribute Σ_e w_{e}, which is P_μ(τ_E < λ) by Proposition 5.1(a). The other attractors
are attractors of U that miss E, with the same exits, by Lemma 4.2(b). ∎

**Corollary 5.3 (no strange loops).** PROVED. If no attractor at any level is a strange
loop, then F_μ is constant on [ω, ω^ω) and ℍ_μ(E) ∈ ℕ ∪ {ω, ∞}.

*Proof.* Every attractor collapses, so K_1 is the Cesàro limit of K_0, and every later
level equals K_1 (Proposition 2.4). The same holds in U^E, whose attractors are
singletons in E or attractors of U (Lemma 4.2(b)). The argument of Lemma 4.2(d), with
d + 1 replaced by 1, shows F_μ is constant from ω. So ℍ is finite, or ω, or ∞. ∎

**Remark 5.4 (the collapse).** INTERPRETATION. Tyler's phrase was "transfinite
probability which collapses out through causal strange loop attractors from
hyperprobability". Here it reads as follows:

- At every limit stage the run's law is a mixture over the attractors its tail can
  settle in (Proposition 2.5).
- A collapsing attractor resolves to its in-universe law π_A. This is where the
  transfinite reduces to ordinary probability.
- The law of τ_E, which is transfinite probability, grows by ordinary steps at
  successor stages. It jumps at a limit only through strange loops (Theorem 5.2).
- It is complete exactly from the stage ℍ_μ(E) on. So hyperprobability is where
  transfinite probability finishes collapsing.

## 6. The Kac bridge

**Lemma 6.1 (Kac).** KNOWN. Let P be an irreducible chain on a finite set with
stationary law π, and let B be a set with π(B) > 0. Let π_B = π(· | B) and let τ⁺_B be
the first return time to B. Then E_{π_B}[τ⁺_B] = 1/π(B) [Kac 1947].

**Theorem 6.2 (the Kac bridge).** PROVED. Let A ∈ 𝒜_0 meet E, and let
ℍ⁺_A(E) = max_{q∈A∩E} ℍ⁺_q(E). In ℚ(ω) (§7),

    π_A(E) ≥ 1/ι(ℍ⁺_A(E)).

Equality holds exactly when ℍ⁺_A(E) = N is finite and every return from A ∩ E to E
takes exactly N steps.

Code: `kac(U, E)` returns a `KacBound` for each such A.

*Proof.* By Proposition 6.3, ℍ⁺_A(E) ∈ ℕ_{≥1} ∪ {ω}. A stationary law of an irreducible
class has full support, so π_A(E) is a positive rational.

- If ℍ⁺_A(E) = ω, then 1/ω is a positive infinitesimal, below every positive rational,
  and the inequality is strict.
- If ℍ⁺_A(E) = N, then τ⁺_E ≤ N almost surely from every q ∈ A ∩ E. Kac's lemma for the
  chain on A with B = A ∩ E gives 1/π_A(E) = E_{π_B}[τ⁺_E] ≤ N. Equality holds exactly
  when τ⁺_E = N almost surely under π_B, which charges every state of A ∩ E. ∎

The bound needs the worst state of A ∩ E. A single start can do better than the
frequency (Proposition 7.5(f)).

**Proposition 6.3.** PROVED. Let q lie in a level-0 attractor A that meets E. Then
ℍ_q(E) ≤ ω and ℍ⁺_q(E) ∈ ℕ_{≥1} ∪ {ω}. Both are finite exactly when the finite-stage
hitting time of E is bounded.

*Proof.* At finite stages the run from q stays in A. A is finite and irreducible, so the
run hits A ∩ E almost surely, P_q(τ_E < ω) = 1 and F_q(ω) = 1. The statement about ℍ⁺
follows by Proposition 4.8. ∎

## 7. Ordinal probability

**Definition 7.1 (ℚ(ω)).** DEFINITION. ℚ(ω) is the field of rational functions in an
indeterminate ω with rational coefficients. It is ordered at infinity: f > 0 when
f(x) > 0 for all large enough real x. Equivalently, the leading coefficients of the
numerator and the denominator have the same sign. Code: `OmegaFraction`, kept in lowest
terms with integer coefficients.

**Proposition 7.2 (ι).** PROVED, using KNOWN facts about natural operations. The map

    ι(ω^m·c_m + ⋯ + c_0) = c_m·ω^m + ⋯ + c_0

embeds (ordinals < ω^ω, <, ⊕, ⊗) into ℚ(ω). It preserves order, ι(α ⊕ β) = ι(α) + ι(β),
and ι(α ⊗ β) = ι(α)·ι(β). Its image is the semiring ℕ[ω]. ℚ(ω) is the fraction field of
ℤ[ω], the ring of differences of ℕ[ω]. The order of ℚ(ω) is the only field order that
extends the order of the ordinals. So ℚ(ω) is the canonical ordered field generated by
the ordinals below ω^ω with their natural operations.

*Proof.*
- Cantor normal forms compare lexicographically from the top, which is the sign at
  infinity of their difference.
- The natural sum adds Cantor normal forms coefficientwise. The natural product
  multiplies them as polynomials, with ω^i ⊗ ω^j = ω^{i+j} [Hessenberg 1906;
  Conway 1976].
- The order of ℤ[ω] is fixed by that of ℕ[ω], since a − b > 0 exactly when a > b.
- An order on an integral domain extends to its fraction field in exactly one way:
  x/y > 0 exactly when xy > 0. ∎

**Remark 7.3.** PROVED. ι does not preserve ordinal addition: 1 + ω = ω, but
ι(1) + ι(ω) = 1 + ω > ω. Guarantees are computed with ordinal addition, which respects
the order in which stages happen. They are compared and inverted in ℚ(ω), which measures
sizes.

**Definition 7.4 (ordinal probability).** DEFINITION.

    p_μ(E) = 1/ι(ℍ⁺_μ(E)) if ℍ⁺_μ(E) < ∞, and p_μ(E) = 0 otherwise.

It is one over the number of evaluations that surely produce E, not counting the start.
ℍ⁺ is used rather than ℍ so that p is defined when μ starts in E, where ℍ = 0. Code:
`ordinal_probability`.

INTERPRETATION: ordinal probability is built from hyperprobability and the field ℚ(ω),
which is the sense in which it is higher order. Hyperprobability is what gives it its
meaning.

**Proposition 7.5.** PROVED; (d)–(f) are also CHECKED.

(a) 0 ≤ p ≤ 1. p = 1 exactly when μK_0 is carried by E, and p = 0 exactly when
ℍ⁺ = ∞.
(b) p is either 1/N with N ∈ ℕ_{≥1}, or a positive infinitesimal (ℍ⁺ ≥ ω), or 0. Its
standard part is 1/N or 0.
(c) p is monotone: E ⊆ E' implies p_μ(E) ≤ p_μ(E'). Also p(∅) = 0 and p(S) = 1.
(d) p is not additive. From tick in the clock, p(tick) + p(tock) = ½ + 1 = 3/2, while
p({tick, tock}) = 1.
(e) p is not maxitive. In the coin, p(H) = p(T) = 1/ω, while p({H, T}) = 1.
(f) The Kac bound holds only for starts spread over the event. If supp μ = A ∩ E for a
level-0 attractor A, then p_μ(E) ≤ π_A(E), by Theorem 6.2 with Propositions 4.3
and 4.8. A single start can exceed the frequency. Take S = {e1, e2, x} with
e1 → e2 → x and x → {e1: ½, x: ½}, and E = {e1, e2}. Then p_{e1}(E) = 1 > ½ = π(E), but
for any law spread over {e1, e2}, p = 1/ω.

*Proof.*
(a) ℍ⁺ ≥ 1, and ℍ⁺ = 1 exactly when F_{μK_0}(0) = 1.
(b) ι(α) is infinite exactly when α ≥ ω.
(c) Monotonicity is Proposition 4.6 through Proposition 4.8: 1 + · is monotone, and
x ↦ 1/x reverses order on positive elements. ℍ(∅) = ∞, and μK_0 is always carried by S.
(d)–(f) Direct computation: in the clock, ℍ⁺(tick) = 2 and ℍ⁺(tock) = 1 from tick; in
the coin, ℍ⁺(H) = 1 + ω = ω; in (f), ℍ⁺_{e1} = 1 + ℍ_{e2} = 1, ℍ⁺_{e2} = 1 + ℍ_x = ω, and
π = (¼, ¼, ½). ∎

So p is a normalized monotone set function, but it is neither a measure nor a
possibility measure.

**Remark 7.6 (surreals).** KNOWN. ℚ(ω) embeds into Conway's surreal numbers by sending ω
to the surreal ω, under which 1/ω is ε. Surreal sum and product restrict to the natural
operations on ordinals [Conway 1976]. Nothing here needs more than ℚ(ω).

**Example 7.7.** CHECKED. The running examples. Each line runs from the given start,
and "frequency" is the in-universe frequency of the event.

| Universe | Event | Start | ℍ | ℍ⁺ | p | P(τ < ω) | frequency | strange loops |
|---|---|---|---|---|---|---|---|---|
| clock | tick | tick | 0 | 2 | 1/2 | 1 | 1/2 | none |
| clock | tock | tick | 1 | 1 | 1 | 1 | 1/2 | none |
| coin | H | T | ω | ω | 1/ω | 1 | 1/2 | none |
| loop | out | a | ω | ω | 1/ω | 0 | 0 | {a, b} at level 0 |
| chain | out | a | ω + 3 | ω + 3 | 1/(ω + 3) | 0 | 0 | {a, b} at level 0 |
| fork | out | a | ω + 2 | ω + 2 | 1/(ω + 2) | 0 | 0 | {a, b} at level 0 |
| swap | out | a | ω² | ω² | 1/ω² | 0 | 0 | {a}, {b} at level 0; {a, b} at level 1 |
| dead end | out | a | ∞ | ∞ | 0 | 0 | 0 | none; P(never) = 1 |

In the loop, in-universe probability gives 0 twice: the run never reaches "out" at a
finite stage, and it spends no time there. Hyperprobability still says the event is
certain by stage ω, and ordinal probability records this as the positive infinitesimal
1/ω.

## 8. Simulations of simulations

**Proposition 8.1 (subdivision).** PROVED (sketch); CHECKED. Fix a finite c ≥ 1 and
build U_c as follows.

- Its states are S × {0, …, c − 1}.
- Its steps are (s, i) → (s, i + 1) for i < c − 1, and (s, c − 1) → (t, 0) with
  probability K_0(s, t).
- Each limit rule sends an attractor to the exit of its phase-0 projection, placed at
  phase 0.

Then for all s and E, and every β = c·α + i with i < c,

    F^{U_c}_{(s,0)}(β; E × {0}) = F^U_s(α; E),

and therefore

    ℍ^{U_c}_{(s,0)}(E × {0}) = c·ℍ^U_s(E), with c·∞ = ∞.

Since c·(ω^m·a_m + ⋯ + ω·a_1 + a_0) = ω^m·a_m + ⋯ + ω·a_1 + c·a_0, a finite slowdown
multiplies only the finite part of a guarantee: c·ω = ω. Code: `refine(ℍ, c)`. The test
builds U_c exactly as above.

*Proof (sketch).* Every β < ω^ω is c·α + i for a unique α and i < c (ordinal division).
Stop U_c at G = E × {0} and start at phase 0.

- Phase-(>0) copies of states of E are never reached.
- c steps from phase 0 are one step of K^E_0, landing at phase 0.
- By induction on k ≤ d, the attractors reachable at level k are the singletons in G
  and the lifts of the attractors of U^E: A × {0, …, c − 1} at level 0, and A × {0} at
  higher levels. They have the same absorption probabilities and exits placed at
  phase 0.

So the stopped law at c·α, for α < ω^{d+1}, is the stopped law of U^E at α placed at
phase 0. The next i < c steps move only mass outside G. At ω^{d+1} both sides are the
limits along ω^d·m (Proposition 5.1(a) and Theorem 5.2; the level-d attractors
collapse). By Lemma 4.2(d), both are constant from there on. ∎

**Remark 8.2 (infinite cost).** INTERPRETATION. `refine(g, cost)` = cost·g for any
ordinal cost ≥ 1. With cost ω, each step of the run needs ω evaluations: 2 becomes ω·2,
and ω becomes ω². This is not constructed as a universe here.

**Proposition 8.3 (self-simulation).** PROVED. For ordinals h ≥ 1 and α, h + α = α
exactly when α ≥ h·ω. So the least α ≥ s with h + α = α is max(s, h·ω). That is ω for
h = 1, and ω² for h = ω or h = ω·2.

*Proof.*
- If α ≥ h·ω, write α = h·ω + γ. Then h + h·ω = h·(1 + ω) = h·ω, so h + α = α.
- If α < h·ω = sup_n h·n, then α < h·n for some n. Suppose h + α = α. By induction
  h·n + α = α, but h·n + α ≥ h·n > α. ∎

INTERPRETATION: a run that must first spend h evaluations and then run itself again
needs a guarantee that absorbs h, and the least one is h·ω. `self_simulation(cost)`
searches for the least fixed point with `ordinatics.calculus.least_fixed_point`. Its
result is a checked fixed point, but its passage through limits is heuristic.

## 9. The family

hyperprobability sits beside Tyler Roost's other public repositories. The links below
are stated as relationships, not as imports, except for ordinatics.

- **ordinatics**, a dependency (ordinatics ≥ 0.3, < 0.4). It provides exact ordinals
  below ω^ω (`Ordinal`, `OMEGA`, the natural sum and product) and
  `calculus.least_fixed_point`. The attractor kinds of Def 3.1 are those of
  `ordinatics.dynamics`. This package adds probability to ordinatics' ordinals: runs
  indexed by them, and events guaranteed at them.
- **hypermath** derives a formal universe from □ under application. Its simulation
  relation `==` relates derivations that reproduce each other's paths. "Simulation
  evaluation" here means one stage of a run of a universe, which is a different notion.
  No claim here uses hypermath's `==`.
- **hyperlogic** and **hyperethics** each derive their structure from one primitive
  operation: `turn` under circuit closure, and `create` under immanence. hyperprobability
  does not derive itself from a primitive. It rests on ordinatics and on the theory of
  finite Markov chains.
- **hyperphysics** states physical laws once, and borrowers cite them by name with a
  declared transport. No physical law is borrowed here. The one borrowed form, the
  strange loop, is declared in hyperphysics' format (Remark 3.3).
- **hyperstratum** defines a hyperfield as "distributed potential/substance across
  recursive form-levels". INTERPRETATION: a universe with its tower of levels
  K_0, K_1, … fits that description. The levels are recursive forms (each is built from
  the attractors of the one below), the law at each stage distributes probability over
  the states, and hyperprobability is read off the whole tower. This is a reading, not a
  formal link.

**Naming.** OPEN. hyperlogic records the family's convention: foundations import
nothing, grounded packages import inward, and a package grounded outside itself carries
a `grounded-` prefix. hyperprobability imports ordinatics, so by that convention its
name could be `grounded-hyperprobability`. The name is Tyler's decision.

## 10. Limits and open questions

- **Stages below ω^ω only.** Longer runs would need rules at infinitely many levels and
  ordinals beyond ordinatics' `Ordinal`. OPEN.
- **Finite state spaces only.** On a countable space a run can escape to infinity or be
  null recurrent, and absorption into closed classes and stationary laws can fail. OPEN.
- **Stage-independent rules.** L_k(A) depends on A but not on γ in γ + ω^k.
  Stage-dependent rules would break the kernel form of Theorem 2.3. OPEN.
- **p is not a measure** (Proposition 7.5). Whether it supports conditioning or
  independence is OPEN.
- **Quantiles.** Let ℍ^δ_μ(E) = min{α : F_μ(α) ≥ 1 − δ}. Then ℍ = sup_{δ>0} ℍ^δ.
  PROVED: ℍ^δ ≤ ℍ for every δ. If β < ℍ, then F(β) < 1, and any δ < 1 − F(β) gives
  ℍ^δ > β. If ℍ = ∞, F stays below 1 from ω^{d+1} on (Lemma 4.2(d)), so ℍ^δ = ∞ for small
  δ. A quantitative theory of ℍ^δ is OPEN: how fast it grows as δ → 0, beyond geometric
  tails like the coin's.
- **Metastability.** CONJECTURE (informal). The level hierarchy resembles the hierarchy
  of cycles and time scales of metastable chains [Freidlin–Wentzell 1984;
  Olivieri–Vares 2005]. A universe whose level-k rules model the exits of level-k cycles
  should put the mass of τ_E at stages with leading term ω^k exactly for events first
  reached on the k-th time scale. ℍ alone does not measure time scales: the coin has
  ℍ = ω and mixes in one step. So the conjecture concerns where the transfinite law puts
  its mass, not ℍ.
- **Axioms and formalization.** An axiomatic characterization of ℍ, and machine-checked
  (for example Lean) proofs of §2–§8, are OPEN.

## References

- M. Kac, On the notion of recurrence in discrete stochastic processes, *Bull. Amer.
  Math. Soc.* 53 (1947), 1002–1010.
- G. Hessenberg, *Grundbegriffe der Mengenlehre*, Göttingen, 1906.
- J. H. Conway, *On Numbers and Games*, Academic Press, 1976.
- J. G. Kemeny and J. L. Snell, *Finite Markov Chains*, Van Nostrand, 1960.
- J. R. Norris, *Markov Chains*, Cambridge University Press, 1997.
- M. I. Freidlin and A. D. Wentzell, *Random Perturbations of Dynamical Systems*,
  Springer, 1984.
- E. Olivieri and M. E. Vares, *Large Deviations and Metastability*, Cambridge
  University Press, 2005.
- D. R. Hofstadter, *Gödel, Escher, Bach*, Basic Books, 1979; *I Am a Strange Loop*,
  Basic Books, 2007.
- Sibling repositories: https://github.com/TimeLordRaps/ordinatics,
  https://github.com/TimeLordRaps/hypermath, https://github.com/TimeLordRaps/hyperlogic,
  https://github.com/TimeLordRaps/hyperethics,
  https://github.com/TimeLordRaps/hyperphysics,
  https://github.com/TimeLordRaps/hyperstratum.

## Appendix: results and their tests

All tests use exact rationals. The random universes of `tests/support.py` have up to six
states, up to three levels of rules, and deterministic, random and absorbing rows. The
oracles `first_certain` and `first_certain_return` search stages directly and share no
code with the recursion.

| Result | Tests |
|---|---|
| Lemma 2.2 | `test_universe.py::test_higher_kernels_absorb_lower_ones`, `test_exact.py::test_absorption_is_harmonic_and_complete` |
| Theorem 2.3 | `test_universe.py::test_composition` |
| Proposition 2.4 | `test_universe.py::test_higher_kernels_absorb_lower_ones`, `test_the_law_can_move_after_the_kernels_settle` |
| Proposition 2.5 | `test_universe.py::test_limit_laws_mix_the_exits_of_the_collapsing_attractors`, `test_the_tail_settles_where_the_run_below_the_limit_is_heading`, `test_collapse_names_the_attractors_that_hold_the_tail` |
| Proposition 2.7 | `test_inuni.py::test_the_in_universe_law_is_the_long_run_average`, `test_universe.py::test_in_uni_is_the_collapsed_limit`, `test_in_uni_ignores_strange_loops` |
| Def 3.1, 3.2, Example 3.4 | `test_universe.py::test_attractors`, `test_loops.py::test_strange_loops`, `test_universe.py::test_laws_at_transfinite_stages` |
| Lemma 4.2 | `test_universe.py::test_stopped_attractors_that_miss_the_event_are_attractors`, `test_transfinite.py::test_the_law_only_grows`, `test_never_is_what_the_last_rule_leaves` |
| Proposition 4.3 | `test_hyper.py::test_initial_laws_take_the_worst_state_of_their_support` |
| Theorem 4.4 | `test_hyper.py::test_the_recursion_finds_the_least_certain_stage` (2000 universes), `test_collapsing_levels_change_nothing`, `test_guarantees_are_infinite_exactly_when_the_event_can_be_missed` |
| Corollary 4.5 | `test_hyper.py::test_guarantees_stay_below_the_first_collapsed_level` |
| Proposition 4.6 | `test_hyper.py::test_guarantees_are_monotone_in_the_event` |
| Proposition 4.8 | `test_hyper.py::test_return_guarantees_against_search` |
| Proposition 4.9 | `test_hyper.py::test_guarantees_depend_only_on_what_is_possible` |
| Proposition 5.1 | `test_transfinite.py::test_before_a_limit_in_deterministic_universes`, `test_before_a_limit_is_squeezed_by_the_run_below_it`, `test_inuni.py::test_reaching_an_event_is_the_law_before_omega` |
| Theorem 5.2 | `test_transfinite.py::test_limit_jumps_come_from_strange_loops`, `test_collapsing_limits_never_jump`, `test_universe.py::test_an_attractor_that_misses_the_event_need_not_survive_stopping` |
| Corollary 5.3 | `test_transfinite.py::test_without_strange_loops_nothing_happens_after_the_finite_stages` |
| Theorem 6.2, Proposition 6.3 | `test_loops.py::test_the_kac_bridge` (against Kac's lemma and search), `test_kac_examples` |
| Proposition 7.2, Remark 7.3 | `test_ordinal.py::test_iota_embeds_the_ordinals_with_their_natural_operations`, `test_iota_does_not_preserve_ordinary_addition`, `test_an_ordered_field`, `test_arithmetic_agrees_with_evaluation` |
| Proposition 7.5 | `test_ordinal.py::test_ordinal_probability_is_not_a_measure` |
| Example 7.7 | `test_hyper.py::test_examples`, `test_ordinal.py::test_ordinal_probability`, `test_examples.py` |
| Proposition 8.1 | `test_loops.py::test_refine_is_the_guarantee_of_a_subdivided_run`, `test_refine` |
| Proposition 8.3 | `test_loops.py::test_self_simulation` |
