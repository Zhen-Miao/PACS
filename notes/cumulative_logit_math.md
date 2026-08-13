# Cumulative-logit PACS — exact likelihood with capture-rate correction

This note replaces the current "stack-and-treat-as-binary" approximation in
`R/PACS_test_cumulative.R` with a proper cumulative-logit formulation that
handles the per-cell capture rate `q_i` consistently with the existing binary
model in `R/PACS_test_logit.R`.

The derivation here is the prerequisite for an implementation; no code is
changed by this document.

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

### 5.2 Observed information under Option A (Louis identity)

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
information `J^{obs}(θ)`.

For the Firth penalty in §6 we use `J^{obs}(θ̂)` evaluated at the current
iterate. This is consistent with Heinze & Schemper's recommendation and
with the existing binary code, which also evaluates `infor_mat` at the
current iterate rather than taking an outer expectation over `M`.

Computing the *expected* Fisher information under Option A would require
an additional outer expectation over `M_i` and is not needed for the
implementation.

## 6. Firth penalty

Carry over the existing penalisation. The penalised log-likelihood is

```
ℓ*(θ) = ℓ(θ) + (1/2) log det I(θ),
```

with `I(θ)` the expected Fisher information from §5.1 under Option B, or
the observed information `J^{obs}(θ)` from §5.2 under Option A. The
penalty score is

```
∂ℓ*/∂θ_r = ∂ℓ/∂θ_r + (1/2) tr( I^{-1} ∂I/∂θ_r ),
```

which generalises `loss_grad_pen`. Because `I(θ)` now mixes `α` and `β`
through a non-trivial cross-block, we cannot factor the working weights
through a single `sqrt(W) · X` matrix as in the binary case. The cleanest
implementation is to assemble `I(θ)` block-by-block from §5.1 and apply the
matrix-derivative identity directly; the Cholesky of `I` can be reused for
both the IRLS step and the Firth gradient.

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
- The Firth penalty in §6 reduces to `loss_grad_pen`.

This gives us a precise parity test to validate the implementation.

## 9. Likelihood-ratio test under the new model

`compare_models` in `R/differential_identification.R` currently computes the
penalised LRT using the binary working weights even when called from the
cumulative wrapper. After the change, the full and null models are fit with
the cumulative likelihood from §3, and the Firth-corrected LRT uses

```
2 [ ℓ*_full(θ̂_full) − ℓ*_null(θ̂_null) ]   ~   χ²_{df},
```

with `df` equal to the number of `β` components constrained to zero under
the null. Both models re-estimate all `α_t` from the same `T`-dimensional
threshold space, so the thresholds contribute zero to the df.

**Boundary caveat.** The χ² approximation assumes the MLE is interior to
the feasible region. When the order constraint `α_1 ≥ ... ≥ α_T` is
active at the optimum (i.e. some `exp(ã_t)` is at the lower clip), the
parameter is on the boundary and the LRT follows a mixture of χ²
distributions per Self–Liang. In practice we should flag fits that hit the
clip and either widen the clip or report a warning rather than a p-value.

**Boundary detection threshold.** `compare_models_cumu` flags peaks where
any `exp(ã_t) < 1e−3` for `t ≥ 2`. At that threshold,
`α_{t−1} − α_t < 0.001`, which is effectively zero on the logit scale —
adjacent cumulative probabilities differ by less than ~0.025 percentage
points near `p = 0.5`. The threshold is conservative; tightening it (say
`1e−4`) would suppress fewer peaks but risk false negatives on borderline
fits. We can expose this as an argument once usage patterns clarify
whether tuning is needed.

### 9.1 Saddlepoint adjustment for scalar tests (df = 1)

The first-order χ² approximation has O(n⁻¹) error in the tail
probabilities. For genome-wide multiple testing at stringent α (e.g.
5 × 10⁻⁶), this error can push the actual rejection rate above the
nominal level. The Barndorff-Nielsen r* statistic achieves O(n⁻³/²)
tail accuracy by combining the signed root of the LRT with a Wald-type
correction.

**Setup.** Partition θ = (λ, ψ) where ψ is the scalar test parameter
(the single β component held to zero under the null) and λ collects
the nuisance parameters (ã₁, ..., ã_T, remaining β's). The penalised
null MLE θ̂*_null satisfies ∂ℓ*/∂λ = 0 while ψ = 0; the penalised
full MLE θ̂*_full satisfies ∂ℓ*/∂θ = 0.

**Signed root.** Define

```
r = sign(ψ̂*_full) · √w,    w = 2[ℓ*(θ̂*_full) − ℓ*(θ̂*_null)].
```

Under H₀, r ~ N(0, 1) to first order; `pchisq(w, 1)` is the two-sided
version.

**Profile score at the null.** Because ∂ℓ*/∂λ = 0 at θ̂*_null, the
profile score for ψ simplifies:

```
S*_{ψ·λ}(θ̂*_null)  =  ∂ℓ*/∂ψ |_{θ̂*_null}.
```

This is the total penalised score (likelihood score + Firth penalty
score) for the test parameter at the constrained null MLE.

**Profile information at the full MLE.** Using the Fisher information
I(θ) as proxy for the negative Hessian (consistent with Fisher-scoring
IRLS), the profile observed information is the Schur complement:

```
J*_{ψψ·λ}  =  I_{ψψ} − I_{ψλ} I_{λλ}⁻¹ I_{λψ},
```

evaluated at θ̂*_full.

**r* formula.**

```
u = S*_{ψ·λ}(θ̂*_null) / √J*_{ψψ·λ}(θ̂*_full),
r* = r + (1/r) · log(u / r).
```

The two-sided p-value is `2 Φ(−|r*|)`.

**Fallback conditions.** The formula is undefined when:
- df > 1 (multivariate test; Skovgaard's 2001 extension would be needed),
- |r| < ε (test statistic ≈ 0; p ≈ 1 regardless),
- u/r ≤ 0 (sign mismatch; indicates numerical instability),
- J*_{ψψ·λ} ≤ 0 (non-positive profile information).

In all cases `compare_models_cumu` falls back to `pchisq(w, df)`.

**Empirical finding: saddlepoint does not improve on Firth-corrected
chi-squared.** Simulation at `n = 300` and `n = 80` (T = 2, df = 1)
shows the Firth-corrected LRT chi-squared is already very well
calibrated (KS p > 0.4 at n = 300, > 0.9 at n = 80). The r*
correction consistently *degrades* calibration (KS p < 0.01) because
Firth's penalty already absorbs the O(n⁻¹) bias that r* targets.

The root cause is that Firth regularization and the Barndorff-Nielsen
r* are **alternative approaches to the same asymptotic deficiency**.
Firth acts on the estimator (penalising the likelihood to reduce MLE
bias), r* acts on the reference distribution (correcting the
chi-squared tail). Both target the O(n⁻¹) error term, and combining
them produces over-correction: the penalised LRT statistic from `ℓ*`
already has Bartlett-type calibration improvements baked in, so the
additional r* adjustment double-counts the correction.

Additionally, the standard r* formula assumes score and information
from the same likelihood; our implementation must use the *unpenalized*
score at the Firth null MLE (where `∂ℓ*/∂λ = 0` but `∂ℓ/∂λ ≠ 0`),
requiring the full profile score formula rather than the simplified
`S_{ψ·λ} = S_ψ`. This adds computational cost without calibration
benefit.

The `pvalue_method = "saddlepoint"` option is retained for research
comparisons and for potential use in non-Firth settings, but
`"chisq"` remains the recommended default. See Kenne Pagui, Salvan &
Sartori (2017) for a comparison of bias-reduction approaches in
standard logistic regression.

## 10. Implementation plan

Phased rollout to keep the existing pipeline working:

1. **Phase 1 (Option B, parity).** Add `R/param_estimate_cumu.R` with
   `loss_fun_cumu`, `loss_gradient_cumu`, `infor_mat_cumu`,
   `loss_grad_pen_cumu`, `irls_iter_cumu`, `irls_iter_cumu_null`. At `T = 1`
   these must match the existing binary functions to machine precision; add a
   test in `tests/testthat` that exercises this parity.
2. **Phase 1 wiring.** Add `method = c("exact", "stacked")` to
   `pacs_test_cumu`. Default remains `"stacked"` for back-compat in this PR.
   Provide a vignette example comparing the two on simulated `T = 2` data.
3. **Phase 2 (LRT).** Update `compare_models` (or add
   `compare_models_cumu`) to compute the penalised LRT with the new
   information matrix, used when `method = "exact"`.
4. **Phase 3 (Option A).** Add a thinning-mode flag and the Louis-identity
   information matrix. Validate against Phase 1 by simulating from each
   model and confirming approximate parity at `q_i ≈ 1`.
5. **Phase 4 (default switch).** Once Phases 1–3 are validated, switch
   `pacs_test_cumu`'s default to `method = "exact"` and deprecate the
   stacking path.

Out of scope for this PR: any changes to the binary `pacs_test_logit`
pipeline.
