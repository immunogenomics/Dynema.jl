
# ---------------------------------------------------------------------------- #
#                    Statistic types and small-sample factors                  #
# ---------------------------------------------------------------------------- #
#
# Dynema reports one of four statistic types, named after the column that
# holds the statistic in the summary table:
#
#                      1 df (one tested term)     > 1 df (joint test)
#   small = false      "z"  (signed, N(0,1))      "χ²"  (χ²_q)
#   small = true       "t"  (signed, t_{G-1})     "F"   (F_{q, G-1})
#
# "z"/"t" are signed, two-sided statistics; "χ²"/"F" are non-negative,
# upper-tailed. `_signed` distinguishes the two families everywhere a
# bootstrap tail count or a reference distribution is needed.
_signed(stattype) = stattype == "z" || stattype == "t"

"""
    smallsample_factors(small, G, nobs) -> (clusteradj, smallsample, dof_r)

Small-sample factors of the cluster-robust score test, replicating
`WildBootTests.jl`'s ML (score) mode with one-way clustering and `small = true`
(`init.jl`: `clusteradj = small`, `smallsample = _dof / (Nobs - robust)`,
`dof_r = G - 1`), which is how Dynema invokes it:

  - `clusteradj = G / (G - 1)` multiplies the cluster-robust denominator
    (the "meat"), where `G` is the number of clusters (donors);
  - `smallsample = nobs / (nobs - 1)` rescales the statistic (its square root
    for 1 df). WildBootTests derives `_dof = nobs - k`, but in ML mode Dynema
    supplies no `predexog`, so its `k` is 0. For single-cell data
    (`nobs` = number of cells) the factor is essentially 1; it is kept so the
    direct implementation reproduces the library to numerical precision;
  - `dof_r = G - 1` is the residual degrees of freedom of the reference
    t (1 df) or F (`q` df) distribution.

With `small = false` all three are neutral: `(1, 1, 0)`, and the reference
distributions are N(0,1) (1 df) and χ² (`q` df).
"""
function smallsample_factors(small::Bool, G::Integer, nobs::Integer)
    small || return (1.0, 1.0, 0.0)
    G >= 2 || throw(ArgumentError("small = true needs at least 2 clusters (donors); got $G"))
    nobs >= 2 || throw(ArgumentError("small = true needs at least 2 observations; got $nobs"))
    return (G / (G - 1), nobs / (nobs - 1), Float64(G - 1))
end

"""
    pass_counts(stat, boot, stattype) -> (c1, c2)

Allocation-free tail counts of `stat` against a bootstrap distribution chunk:
for the signed types ("z", "t"), the (greater, less) counts (whose running
sums stay additive across chunks, unlike their min); for the upper-tailed
types ("χ²", "F"), `(count(boot .>= stat), 0)`.
"""
function pass_counts(stat, boot, stattype)
    if _signed(stattype)
        (count(b -> stat > b, boot), count(b -> stat < b, boot))
    else
        (count(b -> b >= stat, boot), 0)
    end
end

# Final tail count from accumulated (c1, c2) totals
combine_pass(c1, c2, stattype) = _signed(stattype) ? min(c1, c2) : c1

function pass(stat, boot, stattype)
    c1, c2 = pass_counts(stat, boot, stattype)
    combine_pass(c1, c2, stattype)
end

# From a precomputed tail count and total replication count
function compute_pvalue(pass_obs::Integer, nboot::Integer, stattype)

    if _signed(stattype)

        pass_obs == 0 ? 2 / (nboot + 1) : 2 * pass_obs / nboot

    elseif stattype == "χ²" || stattype == "F"

        pass_obs == 0 ? 1 / (nboot + 1) : pass_obs / nboot

    else

        throw(ArgumentError("unknown statistic type \"$stattype\" (expected z, t, χ² or F)"))

    end

end

compute_pvalue(stat, boot, stattype) =
    compute_pvalue(pass(stat, boot, stattype), length(boot), stattype)


"""
    beta_approx_pvalue(stat, boot, stattype, q; dof_r = 0) -> Float64

FastQTL-style beta approximation of the bootstrap p-value, lifting the
empirical `p_boot`'s resolution floor of ~1/B. Both the observed statistic
and every bootstrap draw are transformed to *nominal* p-values under the
analytical reference distribution: two-sided normal for `"z"`, two-sided
Student's t with `dof_r` degrees of freedom for `"t"`, chi-square with `q`
degrees of freedom for `"χ²"`, and F with `(q, dof_r)` degrees of freedom for
`"F"` (`dof_r = G - 1`, see [`smallsample_factors`](@ref)). If that reference were exactly
correct, the bootstrap nominal p-values would be Uniform(0,1) = Beta(1,1);
instead, a `Beta(a, b)` is fit to them by maximum likelihood (with a
method-of-moments fallback), capturing any miscalibration, and the observed
nominal p-value is evaluated under the fitted Beta's CDF. The result is a
smooth p-value that extrapolates below 1/B through the parametric form --
exactly how FastQTL extrapolates beyond its permutation count, applied here
to the score-bootstrap distribution. Returns `NaN` if the fit is degenerate.
"""
function beta_approx_pvalue(stat, boot, stattype, q::Int; dof_r::Real = 0)

    tonominal = if stattype == "z"
        x -> 2 * ccdf(Normal(), abs(x))
    elseif stattype == "t"
        x -> 2 * ccdf(TDist(dof_r), abs(x))
    elseif stattype == "F"
        x -> ccdf(FDist(q, dof_r), x)
    elseif stattype == "χ²"
        x -> ccdf(Chisq(q), x)
    else
        throw(ArgumentError("unknown statistic type \"$stattype\" (expected z, t, χ² or F)"))
    end

    ps = clamp.(tonominal.(boot), 1e-12, 1 - 1e-12)

    bfit = try
        fit_mle(Beta, ps)
    catch
        # Method-of-moments fallback (also FastQTL's starting values)
        m, v = mean(ps), var(ps)
        (v <= 0 || v >= m * (1 - m)) && return NaN
        k = m * (1 - m) / v - 1
        Beta(max(m * k, 1e-8), max((1 - m) * k, 1e-8))
    end

    return cdf(bfit, min(tonominal(stat), 1.0))

end


"""
    crvetest(R, r; resp, scores, betas, A, clustid, small = false)

Analytical CRVE score (Lagrange multiplier) test via
`WildBootTests.scoretest`. WildBootTests is an *optional* dependency: this
generic implementation is only needed as a fallback (multiway clustering,
formulas outside the fast path) and for the fast path's one-time
self-check; standard one-way-clustered mapping uses Dynema's own
[`crvetest_direct`](@ref) instead. The real method lives in the
DynemaWildBootTestsExt extension; this stub errors with install
instructions when the extension isn't loaded.
"""
function crvetest(R, r; kwargs...)
    error("This code path needs the optional WildBootTests package (used as a fallback " *
          "for multiway clustering or non-standard formulas). Install and load it with:\n" *
          "    using Pkg; Pkg.add(\"WildBootTests\"); using WildBootTests\n" *
          "alongside `using Dynema` (the dynema-map CLI does this automatically when --boot is used).")
end


# One-way CRVE score-test ingredients shared by every analytical routine:
# the numerator `numer = (A R')' S` (a one-step Newton estimate of the tested
# coefficient(s) from the null fit) and the recentered CRVE "meat" transformed
# by `A R'`, `denom = K'K` (no small-sample factor applied).
function _score_numer_denom(R::AbstractMatrix, A::AbstractMatrix,
                            Sg::AbstractMatrix, clustshare::AbstractVector)

    ARt = A * R'                          # k × q
    S = vec(sum(Sg, dims = 1))            # total score sum (k)
    numer = ARt' * S                      # q

    K = Sg * ARt                          # G × q
    K .-= clustshare * (S' * ARt)         # recenter, as WildBootTests does for scorebs
    denom = K' * K                        # q × q CRVE denominator

    return numer, denom

end

# p-value of a signed 1-df statistic: N(0,1) when dof_r == 0, t(dof_r) otherwise
_p1(stat, dof_r) = dof_r == 0 ? ccdf(Chisq(1), stat^2) : 2 * ccdf(TDist(dof_r), abs(stat))

"""
    crvetest_direct(R, A, Sg, clustshare; small = false, nobs = 0)

Direct implementation of the analytical CRVE score test for the shared-null
fast path, replicating `WildBootTests.scoretest`'s ml-mode statistic
(`reps = 0`, `imposenull = true`, `scorebs = true`, one-way clustering)
without the library's per-call setup:

  - `numer = (A R')' S`, with `S` the total score sum;
  - cluster-level `K = Sg * (A R')`, recentered by each cluster's share of
    observations, so `denom = K'K` is the CRVE meat transformed by `A R'`;
  - `small = false` (default): q = 1: `z = numer / sqrt(denom)`,
    `p = ccdf(Chisq(1), z²)` (two-tailed); q > 1: `χ² = numer' denom⁻¹ numer`,
    `p = ccdf(Chisq(q), χ²)`;
  - `small = true`: with `(cm, ss, dof_r) = smallsample_factors(true, G, nobs)`
    (`G = size(Sg, 1)` clusters), the denominator becomes `cm * denom`,
    q = 1: `t = sqrt(ss) * numer / sqrt(cm * denom)`,
    `p = 2 * ccdf(TDist(G - 1), |t|)`; q > 1:
    `F = (ss / q) * numer' (cm * denom)⁻¹ numer`,
    `p = ccdf(FDist(q, G - 1), F)`.

`Sg` is the G × k matrix of per-cluster score sums; `clustshare[g]` is
cluster g's share of observations; `nobs` (number of observations/cells) is
required when `small = true`.

Returns `(stat, p, stattype, beta, se, dof_r)` in the same shape `crvetest`
results are consumed (plus `pk` for q > 1). `beta = numer` is the one-step
(score-based) estimate of each tested coefficient from the null fit and
`se = sqrt(cm * diag(denom) / ss)` its cluster-robust standard error, so that
`beta ./ se` reproduces the reported z/t statistic of each tested term (for
q = 1, exactly `stat`). `stattype` is `"z"`/`"t"` for q = 1 and `"χ²"`/`"F"`
for q > 1.
"""
function crvetest_direct(R::AbstractMatrix, A::AbstractMatrix,
                         Sg::AbstractMatrix, clustshare::AbstractVector;
                         small::Bool = false, nobs::Integer = 0)

    numer, denom = _score_numer_denom(R, A, Sg, clustshare)
    cm, ss, dof_r = smallsample_factors(small, size(Sg, 1), nobs)

    q = size(R, 1)
    se = sqrt.(cm .* diag(denom) ./ ss)   # SE of each tested coefficient

    if q == 1
        stat = numer[1] / se[1]
        return (stat = stat, p = _p1(stat, dof_r), stattype = small ? "t" : "z",
                beta = numer, se = se, dof_r = dof_r)
    else
        χ² = dot(numer, Symmetric(cm .* denom) \ numer)
        # Per-component 1-df tests of each tested coefficient, from the same
        # restricted fit: (numer_k / se_k) against the 1-df reference (see
        # crve_percomponent). Free by-products of the joint statistic;
        # consumed by the per-context interaction columns.
        pk = [_p1(numer[k] / se[k], dof_r) for k in 1:q]
        stat = small ? (ss / q) * χ² : χ²
        p = small ? ccdf(FDist(q, dof_r), stat) : ccdf(Chisq(q), stat)
        return (stat = stat, p = p, stattype = small ? "F" : "χ²", pk = pk,
                beta = numer, se = se, dof_r = dof_r)
    end

end


# Per-cluster score sums (G × k) and each cluster's share of observations,
# from row-level scores and one-way integer cluster ids 1..G.
function _cluster_score_sums(scores::AbstractMatrix, clustid::AbstractVector{<:Integer})

    G, n = maximum(clustid), length(clustid)
    Sg = zeros(G, size(scores, 2))
    clustshare = zeros(G)
    @inbounds for i in 1:n
        g = clustid[i]
        clustshare[g] += 1.0
        @views Sg[g, :] .+= scores[i, :]
    end
    clustshare ./= n
    return Sg, clustshare

end


"""
    crve_percomponent(R, A, Sg, clustshare; small = false, nobs = 0) -> Vector{Float64}
    crve_percomponent(R, A; scores, clustid, small = false) -> Vector{Float64}

Per-component 1-df CRVE score p-values for each tested row of `R`, computed
from the *same* restricted (null) fit as the joint test: with `numer` and
`denom = K'K` exactly as in [`crvetest_direct`](@ref), and
`se_k = sqrt(cm * denom[k,k] / ss)` (factors from
[`smallsample_factors`](@ref); both 1 when `small = false`),

    p_k = 2 * ccdf(Normal(), |numer[k] / se_k|)     (small = false)
    p_k = 2 * ccdf(TDist(G - 1), |numer[k] / se_k|)  (small = true)

These are tests of each coefficient under the joint null (all tested
coefficients zero) -- exact size under that null, but under partial
alternatives with correlated contexts, a true effect in one context can
project onto another's score. They decompose the joint statistic; the joint
test remains the primary inference.

The keyword method builds per-cluster score sums from row-level `scores`
and one-way integer cluster ids (used by the generic, non-workspace path).
"""
function crve_percomponent(R::AbstractMatrix, A::AbstractMatrix,
                           Sg::AbstractMatrix, clustshare::AbstractVector;
                           small::Bool = false, nobs::Integer = 0)

    numer, denom = _score_numer_denom(R, A, Sg, clustshare)
    cm, ss, dof_r = smallsample_factors(small, size(Sg, 1), nobs)
    return [_p1(numer[k] / sqrt(cm * denom[k, k] / ss), dof_r) for k in 1:size(R, 1)]

end

"""
    crve_se(X, y, β, groups; small = false) -> Vector{Float64}

One-way cluster-robust (CRVE) standard errors of the *unrestricted*
Poisson/log-link coefficients `β` (the `se_wald_<term>` output columns, reported
alongside the unrestricted estimates selected with `betas`).

Sandwich estimator `V = (cm / ss) A (Σ_g s_g s_g') A`, where

  - `A = (X'diag(μ̂)X)⁻¹` is the model-based bread,
  - `s_g` is the sum of the scores `(y - μ̂) .* X` over the cells of cluster
    (donor) `g`,
  - everything is evaluated at the unrestricted fit `μ̂ = exp.(Xβ)`,
  - `cm = ss = 1` when `small = false`; with `small = true` they are the same
    small-sample factors as the score test
    ([`smallsample_factors`](@ref): `cm = G/(G-1)`, `ss = nobs/(nobs-1)`),
    so that SEs are on the same footing as the reported t statistics.

The SE of coefficient `j` is `sqrt(V[j, j])`, on the same (natural-log) scale
as `β[j]`. `β[j] / se[j]` is therefore a Wald z-statistic (`small = false`)
or a Wald t-statistic on `G - 1` degrees of freedom (`small = true`).

!!! note
    These SEs are *not* the quantities behind Dynema's p-values. The `z`/`t`
    and `p` columns come from a score test at the restricted (null) fit, with
    recentered scores (see [`crvetest_direct`](@ref)); this SE comes from the
    unrestricted fit, so `|β| / se` will generally differ from `|z|` (they
    agree asymptotically, and closely for small effects). With
    `small = false` no small-sample correction is applied, matching the score
    test, so SEs are optimistic with few donors. The score-test `p` (or
    `p_boot`) remains the primary inference; use `se_wald_<term>` for
    effect-size uncertainty and meta-analysis.

Returns `NaN`s under multiway clustering (`size(groups, 2) != 1`), which is
not implemented here.
"""
function crve_se(X::AbstractMatrix, y::AbstractVector, β::AbstractVector, groups::AbstractMatrix;
                 small::Bool = false)

    k = size(X, 2)
    size(groups, 2) == 1 || return fill(NaN, k)

    clustid = Int.(vec(groups[:, 1]))
    μ = exp.(X * β)
    A = inv(Symmetric(X' * (μ .* X)))
    scores = (y .- μ) .* X
    Sg = zeros(maximum(clustid), k)
    @inbounds for i in eachindex(clustid), j in 1:k
        Sg[clustid[i], j] += scores[i, j]
    end
    cm, ss, _ = smallsample_factors(small, size(Sg, 1), length(clustid))
    V = (cm / ss) .* (A * (Sg' * Sg) * A)
    return sqrt.(max.(diag(V), 0.0))

end

# Which unrestricted coefficients get an `se_wald_<term>` column: the
# genotype terms only (`G` and every `G & <context>`); nuisance covariates
# and context main effects are reported as estimates without SEs.
is_genotype_term(name::AbstractString) = name == "G" || startswith(name, "G & ")

# Wide one-row DataFrame of SE columns named `se_wald_<coefname>` (see `crve_se`)
# for the genotype terms, column-aligned with the coefficient columns.
function se_frame(se::AbstractVector, names)
    keep = findall(is_genotype_term, String.(names))
    DataFrame(transpose(se[keep]), "se_wald_" .* String.(names)[keep])
end

function crve_percomponent(R::AbstractMatrix, A::AbstractMatrix;
                           scores::AbstractMatrix, clustid::AbstractVector{<:Integer},
                           small::Bool = false)

    Sg, clustshare = _cluster_score_sums(scores, clustid)
    return crve_percomponent(R, A, Sg, clustshare; small = small, nobs = length(clustid))

end
