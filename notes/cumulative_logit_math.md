# Cumulative-logit PACS — exact likelihood with capture-rate correction

This note replaces the current "stack-and-treat-as-binary" approximation in
`R/PACS_test_cumulative.R` with a proper cumulative-logit formulation that
handles the per-cell capture rate `q_i` consistently with the existing binary
model in `R/PACS_test_logit.R`.

The implemented `method = "exact"` path uses the unpenalized Option B
likelihood derived here. Bias reduction for this curved model is deferred
until its adjusted scores have been derived explicitly.

## 1. Notation

- Cells indexed by `i = 1, ..., n`, with covariate row `x_i ∈ R^p` and
  per-cell capture rate `q_i ∈ (0, 1]`.
- Latent (uncaptured) ordinal response `Y_i ∈ {0, 1, ..., T}`. In snATAC,
  `Y_i` is the count of accessible fragments at a peak before sequencing
  thinning.
- Observed response `M_i ∈ {0, 1, ..., T}` after capture by the sequencing
  protocol.
- Threshold parameters `α = (α_1, ..., α_T)` with the proportional-odds
  ordering `α_1 ≥ α_2 ≥ ... ≥ α_T`.
- Slope `β ∈ R^p`. We collect `θ = (α^T, β^T)^T`.
- Logistic link `σ(η) = 1 / (1 + e^{-η})`.

For each `t = 1, ..., T` define

```
p_{it} = Pr(Y_i ≥ t | x_i) = σ(α_t + x_i^T β),
```

with conventions `p_{i,0} ≡ 1` and `p_{i,T+1} ≡ 0`. The category
probabilities are

```
π_{ik} = p_{ik} − p_{i,k+1},   k = 0, 1, ..., T.
```

Useful derivatives:

```
∂p_{it}/∂β  = p_{it}(1 − p_{it}) x_i ≡ u_{it} x_i,
∂p_{it}/∂α_s = u_{it} · 1[s = t].
```

## 2. Capture model

The binary PACS model uses `Pr(M = 1) = q_i p_i`, `Pr(M = 0) = 1 − q_i p_i`.
There are two natural extensions to the ordinal case; both reduce to the
binary model at `T = 1`.

### Option A — fragment-level thinning (matches biology)

Each underlying fragment is independently captured with probability `q_i`,
i.e. `M_i | Y_i = k ~ Binomial(k, q_i)`. Then

```
Pr(M_i = m | x_i, q_i)
   = Σ_{k = m}^{T} C(k, m) q_i^m (1 − q_i)^{k − m} · π_{ik},
```

where `C(k, m) = k! / (m! (k − m)!)`.

### Option B — cell-level all-or-nothing capture (computationally simpler)

With probability `q_i` the cell is observed perfectly and `M_i = Y_i`; with
probability `1 − q_i` the cell drops out completely and `M_i = 0`. Then

```
Pr(M_i = 0 | x_i, q_i) = 1 − q_i p_{i1},
Pr(M_i = m | x_i, q_i) = q_i (p_{im} − p_{i,m+1}),   1 ≤ m ≤ T.
```

### Tradeoff

Option A is the natural generalisation of the binary capture model and the
correct statistical treatment of fragment-level dropout. Option B is more
tractable (no Binomial mixture), and at high counts where capture is the
limiting factor it understates the rate at which a high-`Y_i` cell can produce
a moderate `M_i`. Empirically, for `T = 2` the two differ only when
`p_{i2}` is appreciable; for snATAC peak-level data where counts ≥ 2 are rare,
Option B is usually adequate.

A subtle interpretive point: Option B forces "low observed counts" to mean
"truly less accessible" plus "truly inaccessible at all", with the dropout
mass parked entirely on `M = 0`. Option A spreads the dropout mass across
`m = 0, 1, ..., k − 1` proportionally to `Binomial(k, q_i)`, so a cell with
`Y_i = 3` that loses one fragment still contributes to `M_i = 2`, not to
`M_i = 0`.

We will implement **Option B first** as a drop-in replacement (parity check at
`T = 1` against existing binary PACS) and add Option A as a second pass.

## 3. Log-likelihood

Per cell, log-density `ℓ_i(θ) = log Pr(M_i = m_i | x_i, q_i, θ)`.

### 3.1 Option B

```
ℓ_i =  1[m_i = 0]   · log(1 − q_i p_{i1})
     + 1[m_i ≥ 1]   · [ log q_i + log(p_{i,m_i} − p_{i,m_i + 1}) ].
```

The `log q_i` term depends only on data and `q_i`; it is constant in `θ` for
fixed data and can be dropped from optimisation.

### 3.2 Option A

```
ℓ_i = log Σ_{k = m_i}^{T} w_{ik} (p_{ik} − p_{i,k+1}),
w_{ik} = C(k, m_i) q_i^{m_i} (1 − q_i)^{k − m_i}.
```

Define the posterior over the latent count given the observed count:

```
γ_{ik} = w_{ik} (p_{ik} − p_{i,k+1}) / Σ_{j ≥ m_i} w_{ij} (p_{ij} − p_{i,j+1}),
         k = m_i, ..., T,    Σ_k γ_{ik} = 1.
```

`γ_{ik}` will play the role of E-step posteriors when computing the score,
Hessian, and Fisher information.

## 4. Score equations

### 4.1 Option B

Let `Δ_{ik} = p_{ik} − p_{i,k+1}` and `D_{ik} = u_{ik} − u_{i,k+1}` (with
`u_{i,0} = u_{i,T+1} = 0`). For cells with `m_i ≥ 1`:

```
∂ℓ_i/∂β   = D_{i,m_i} / Δ_{i,m_i} · x_i,
∂ℓ_i/∂α_t = (u_{i,m_i}   / Δ_{i,m_i}) · 1[t = m_i]
            − (u_{i,m_i+1} / Δ_{i,m_i}) · 1[t = m_i + 1].
```

(The numerator differs in the two indicator branches — `u_{i,m_i}` for the
`t = m_i` term and `u_{i,m_i+1}` for the `t = m_i + 1` term — because
`∂p_{ik}/∂α_t = u_{ik}·1[k = t]`.)

For cells with `m_i = 0`:

```
∂ℓ_i/∂β   = − q_i u_{i1} / (1 − q_i p_{i1}) · x_i,
∂ℓ_i/∂α_1 = − q_i u_{i1} / (1 − q_i p_{i1}),
∂ℓ_i/∂α_t = 0   for t ≥ 2.
```

Edge case `m_i = T`: by convention `p_{i,T+1} = 0` and `u_{i,T+1} = 0`,
so `Δ_{i,T} = p_{iT}` and `D_{i,T} = u_{iT}`. The score reduces to
`(u_{iT}/p_{iT}) x_i = (1 − p_{iT}) x_i` for `β`, and `u_{iT}/p_{iT}` for
`α_T`.

At `T = 1` these collapse to the existing binary PACS score in
`loss_gradient`.

### 4.2 Option A

Let `r_{ik} = (u_{ik} − u_{i,k+1}) / (p_{ik} − p_{i,k+1})`. Then

```
∂ℓ_i/∂β   = ( Σ_{k ≥ m_i} γ_{ik} · r_{ik} ) · x_i,
∂ℓ_i/∂α_t = γ_{it} · u_{it} / Δ_{it}
            − γ_{i,t−1} · u_{it} / Δ_{i,t−1},
```

with the convention `γ_{i,m_i − 1} ≡ 0`.

## 5. Observed information / Hessian

Write `L_i = exp(ℓ_i)`. Then

```
−∂² log L_i / ∂θ ∂θ^T = − L_i^{-1} ∂² L_i / ∂θ ∂θ^T  +  S_i S_i^T,
```

where `S_i = ∂ℓ_i/∂θ`. The needed second derivatives of `L_i` involve only

```
∂u_{it}/∂η_{it} = u_{it} (1 − 2 p_{it}),
```

with `η_{it} = α_t + x_i^T β`. All terms are closed-form.

### 5.1 Fisher information under Option B (block structure)

We compute the *expected* information

```
I(θ) = E_M [ S(θ) S(θ)^T ] = Σ_i Σ_m Pr(M_i = m | x_i, q_i, θ) · s_{i,m} s_{i,m}^T,
```

where `s_{i,m}` is the per-cell score evaluated at `M_i = m`.

The category probabilities are `Pr(M_i = 0) = 1 − q_i p_{i1}` and
`Pr(M_i = m) = q_i Δ_{i,m}` for `m ≥ 1`. Define

```
ν_{i,0} = q_i^2 u_{i1}^2 / (1 − q_i p_{i1})           (contribution from m = 0),
ν_{i,k} = q_i / Δ_{i,k}                                (m = k contribution scale, k ≥ 1).
```

Under Option B the blocks of `I(θ)` are as follows.

**α-block (tridiagonal in `t`).** For `t = 1`:

```
I_{α_1, α_1} = Σ_i [ q_i u_{i1}^2 / Δ_{i,1}  +  ν_{i,0} ].
```

For `2 ≤ t ≤ T`:

```
I_{α_t, α_t} = Σ_i [ q_i u_{it}^2 / Δ_{i,t}  +  q_i u_{it}^2 / Δ_{i,t−1} ].
```

Off-diagonal entries (only adjacent `t` couple, because score for `α_t`
is non-zero only when `m_i ∈ {t − 1, t}`):

```
I_{α_t, α_{t+1}} = − Σ_i q_i u_{it} u_{i,t+1} / Δ_{i,t}    for t = 1, ..., T − 1,
I_{α_t, α_s}     = 0    for |t − s| ≥ 2.
```

The `m_i = 0` mass enters `I_{α_1, α_1}` via `ν_{i,0}` and nowhere else,
because `∂ℓ_i/∂α_t = 0` at `m_i = 0` for `t ≥ 2`.

**β-block:**

```
I_{β, β} = Σ_i x_i x_i^T · v_i,
v_i = Σ_{k = 1}^{T} q_i D_{ik}^2 / Δ_{i,k}  +  ν_{i,0}.
```

**α–β cross-block.** For `t = 1`:

```
I_{α_1, β} = Σ_i x_i · [ q_i u_{i1} D_{i,1} / Δ_{i,1}  +  ν_{i,0} ].
```

For `2 ≤ t ≤ T`:

```
I_{α_t, β} = Σ_i x_i · [ q_i u_{it} D_{i,t} / Δ_{i,t}
                       − q_i u_{it} D_{i,t−1} / Δ_{i,t−1} ].
```

(Sign comes from `∂ℓ_i/∂α_t |_{m_i = t − 1} = − u_{it}/Δ_{i,t−1}`.)

**Reduction at `T = 1`.** Then `Δ_{i,1} = p_{i1}`, `D_{i,1} = u_{i1}`, the
α-block is a scalar `α_1`, and

```
I_{β, β} |_{T = 1}
   = Σ_i x_i x_i^T · [ q_i u_{i1}^2 / p_{i1}  +  q_i^2 u_{i1}^2 / (1 − q_i p_{i1}) ]
   = Σ_i x_i x_i^T · q_i u_{i1}^2 · [ 1 / p_{i1} + q_i / (1 − q_i p_{i1}) ].
```

Combine the bracket as `[ (1 − q_i p_{i1}) + q_i p_{i1} ] / [ p_{i1} (1 − q_i p_{i1}) ]
= 1 / [ p_{i1} (1 − q_i p_{i1}) ]`, and substitute `u_{i1}^2 = p_{i1}^2 (1 − p_{i1})^2`, giving

```
I_{β, β} |_{T = 1} = Σ_i x_i x_i^T · q_i p_{i1} (1 − p_{i1})^2 / (1 − q_i p_{i1}),
```

which matches `wii = q_i p_i (1 − p_i)^2 / (1 − q_i p_i)` in `infor_mat`
exactly.

### 5.2 Information under Option A

For Option A we use Louis' identity to avoid expanding the
mixture-derivative algebra. Louis gives the *observed* information:

```
J_i^{obs}(θ) = E_{γ_i}[ J_i^{complete}(θ) ]
             − E_{γ_i}[ S_i^{complete}(θ) S_i^{complete}(θ)^T ]
             + S_i(θ) S_i(θ)^T,
```

where the "complete-data" score and information `S^{complete}, J^{complete}`
are the standard cumulative-logit quantities for the latent `Y_i`
(McCullagh, 1980), and the expectations are taken under the latent-class
posterior `γ_{ik}` from §3.2. Summed over `i` this gives the observed
information `J^{obs}(θ)`. It is useful for Newton steps, but evaluating it at
the current iterate does not turn it into expected Fisher information.

If Option A is implemented, its expected information is obtained directly
by enumerating the observed categories:

```
I_i(θ) = Σ_{m=0}^{T} Pr(M_i = m | x_i, q_i, θ) s_{i,m}(θ) s_{i,m}(θ)^T,
I(θ)   = Σ_i I_i(θ),
```

where `s_{i,m}` is the observed-data score from §4.2 evaluated with
`M_i = m`. Thus Louis information may be used for optimization while the
outer expectation above supplies Fisher information when Fisher information
is specifically required.

## 6. Bias reduction and current estimation policy

The exact path currently maximizes the **unpenalized** Option B likelihood.
It does not claim Firth mean-bias reduction, unbiased slope estimates, or a
Firth-corrected likelihood-ratio test.

An earlier implementation added the Jeffreys-type objective

```
ℓ(θ) + (1/2) log det I(θ).
```

That objective was not a derived Firth adjustment for this nonlinear/curved
capture-adjusted cumulative-logit model. Firth's equivalence between the
adjusted score and the Jeffreys penalty has parameterization restrictions
([Firth, 1993](https://doi.org/10.1093/biomet/80.1.27)). Established bias
reduction for cumulative-link models instead starts from the first-order
bias and expected information and derives the adjusted scores
([Kosmidis, 2014](https://doi.org/10.1111/rssb.12025)).

There is also a concrete parameterization problem. Let `G = ∂(α,β) /
∂(ã,β)`. Applying the determinant penalty after the nonlinear order map gives

```
(1/2) log det(G^T I_{α,β} G)
  = (1/2) log det(I_{α,β}) + log |det G|
  = (1/2) log det(I_{α,β}) + Σ_{s=2}^{T} ã_s.
```

The last term is an additional threshold-gap penalty, not an optimization
constant. That defect is separable from the bias-reduction question: computing
the determinant from `I_(α,β)` would remove the extra term in one line. A
Jeffreys penalty in `(α, β)` was **not evaluated and rejected** as an estimator
in this work; it was deliberately not substituted because it still would not
establish Firth mean-bias reduction for the curved capture-adjusted model.
Removing the penalty makes the current MLE objective invariant to the order
reparameterization. Any future bias-reduced estimator must state its target
parameterization and derive the appropriate adjusted scores before it is
exposed through the API.

The same limitation applies to the legacy binary capture-adjusted path:
`loss_grad_pen()` is historically described as Firth correction, but its
Jeffreys-type adjustment has not been shown here to remove first-order bias
when `Pr(M_i = 1) = q_i sigma(eta_i)`. It remains unchanged for backward
compatibility, not as a validated reference implementation of Firth bias
reduction. The exact Option B path therefore uses the unpenalized likelihood
for both estimation and inference rather than treating the binary penalty as
settled statistical ground truth.

## 7. Order-constraint / identifiability

`α_1 ≥ α_2 ≥ ... ≥ α_T` is required for valid probabilities. We
reparameterise as

```
α_1 = ã_1,
α_t = ã_1 − Σ_{s = 2}^{t} exp(ã_s),    t = 2, ..., T,
```

so `ã ∈ R^T` is unconstrained and `α_{t−1} − α_t = exp(ã_t) > 0` for
`t ≥ 2`, which strictly enforces the ordering.

Score and Hessian transform via the Jacobian `J = ∂α/∂ã`, which is
lower-triangular with column 1 of all ones, `J_{t,s} = − exp(ã_s)` for
`2 ≤ s ≤ t`, and zero above the diagonal:

```
score_ã = J^T · score_α,
H_ã     = J^T · H_α · J  +  Σ_t (∂J/∂ã)_t · score_α[t]    (Newton),
```

with `(∂²α_t/∂ã_s ∂ã_s) = − exp(ã_s)` along the diagonal of the curvature
correction. The correction vanishes at the optimum and can be dropped for
Fisher scoring.

The order map is used only to enforce valid probabilities during numerical
optimization. Because the implemented objective is the unpenalized
likelihood, maximizing in `(ã, β)` and mapping back gives the same MLE as
constrained maximization in `(α, β)` for an interior solution.

Initial values: warm-start with `β = 0`, set `α̂_t` by inverting the
empirical cumulative rate after capture correction, then map back to `ã`:

```
α̂_t = logit( clip( mean_i 1[m_i ≥ t] / mean_i q_i,  ε,  1 − ε ) ),
ã_1 = α̂_1,
ã_t = log( max( ε,  α̂_{t−1} − α̂_t ) )    for t = 2, ..., T.
```

This matches the spirit of the existing `par_initial = 0.05` warm start when
`T = 1`.

## 8. Reduction to the binary model

When `T = 1`:

- Option A and Option B coincide. (At `T = 1`, the only Binomial mixture
  weights are `w_{i,1} = q_i` for `m = 1` and `w_{i,0} = 1 − q_i p_{i1}`
  marginal mass at `m = 0`, recovering the Bernoulli model.)
- The α-block is a scalar `α_1`; concatenating `(α_1, β)` gives the
  existing PACS design `θ = (intercept, β)`.
- The score in §4.1 reduces to `loss_gradient` term-for-term (verified
  above for both `m = 0` and `m = 1`).
- The Fisher information in §5.1 reduces to `infor_mat` (verified by the
  algebra at the end of §5.1: `q_i p_{i1} (1 − p_{i1})^2 / (1 − q_i p_{i1})`).

This gives precise likelihood, score, and expected-information parity tests.
The exact optimizer itself is intentionally unpenalized and therefore does
not match the legacy binary Firth optimizer.

## 9. Likelihood-ratio test under the new model

The full and null models are fit with the cumulative likelihood from §3. For
converged interior fits, `compare_models_cumu` uses the ordinary LRT

```
2 [ ℓ_full(θ̂_full) − ℓ_null(θ̂_null) ]   ~   χ²_{df},
```

with `df` equal to the number of `β` components constrained to zero under
the null. Both models re-estimate all `α_t` from the same `T`-dimensional
threshold space, so the thresholds contribute zero to the df.

Inference is withheld (`NA`) when either fit did not converge, when the LRT
statistic is materially negative, or when either fit is on the order
boundary. Tiny negative statistics within a scale-aware numerical tolerance
are rounded to zero; a material negative value indicates that the fitted
full model failed to attain the nested null likelihood and is not converted
silently to `p = 1`.

**Sparse-peak behavior.** The unpenalized MLE can fail to exist or its
information can become singular when a peak has very few nonzero observations.
In review simulations the exact-path `NA` rate was 0% for dense peaks, 3% at
`alpha = (-2.5, -4)` with `n = 300`, 38% at `alpha = (-3.5, -5)` with
`n = 300`, and 46% at `alpha = (-2.5, -4)` with `n = 100`. These figures
describe those simulation regimes rather than a universal sparsity curve.
Withholding inference is conservative for the affected peak because a failed
fit cannot become a false positive, but the increasing `NA` rate reduces
power by removing peaks from analysis. Simulated type-I error among returned
p-values was not inflated in the reviewed regimes.

**Boundary caveat.** The χ² approximation assumes the MLE is interior to
the feasible region. When the order constraint `α_1 ≥ ... ≥ α_T` is active
at the optimum, the reference distribution generally involves a
problem-specific mixture. The direction of the resulting distortion is not
universal. Consequently, `compare_models_cumu` checks both full and null
fits, warns, and returns `NA` rather than an ordinary chi-square p-value.

**Boundary detection threshold.** `compare_models_cumu` flags peaks where
any `exp(ã_t) < 1e−3` for `t ≥ 2`. At that threshold,
`α_{t−1} − α_t < 0.001`, which is effectively zero on the logit scale —
adjacent cumulative probabilities differ by less than ~0.025 percentage
points near `p = 0.5`. This is a numerical guard, not a test of whether a gap
is statistically distinguishable from zero; that would require its sampling
uncertainty. Tightening it (say `1e−4`) flags fewer numerically collapsed
gaps. It is an internal argument to `compare_models_cumu` and is not currently
part of the public API.

## 10. Implementation status and remaining work

- Option B likelihood, score, expected information, ordered-threshold map,
  unpenalized Fisher scoring, and ordinary LRT are implemented.
- `max_T` is a top-code in the public exact API: category `T` means `T` or
  more. Internal likelihood helpers reject responses outside `0:T` clearly.
- Logistic category differences are computed in log space and scoring uses
  stable ratios. Fisher-scoring updates use step-halving, and convergence
  requires both a small step and a small score on the free parameters.
- The default remains `method = "stacked"` for backward compatibility.
- Option A and a statistically derived mean-bias-reduction method remain
  future work. Option A must use the expected-information construction in
  §5.2 whenever Fisher information is required.

Out of scope: changes to the legacy binary `pacs_test_logit` pipeline.
