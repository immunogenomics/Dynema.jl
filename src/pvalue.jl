
"""
    pass_counts(stat, boot, stattype) -> (c1, c2)

Allocation-free tail counts of `stat` against a bootstrap distribution chunk:
for "z", the (greater, less) counts (whose running sums stay additive across
chunks, unlike their min); for "χ²", `(count(boot .>= stat), 0)`.
"""
function pass_counts(stat, boot, stattype)
    if stattype == "z"
        (count(b -> stat > b, boot), count(b -> stat < b, boot))
    else
        (count(b -> b >= stat, boot), 0)
    end
end

# Final tail count from accumulated (c1, c2) totals
combine_pass(c1, c2, stattype) = stattype == "z" ? min(c1, c2) : c1

function pass(stat, boot, stattype)
    c1, c2 = pass_counts(stat, boot, stattype)
    combine_pass(c1, c2, stattype)
end

# From a precomputed tail count and total replication count
function compute_pvalue(pass_obs::Integer, nboot::Integer, stattype)

    if stattype == "z"

        pass_obs == 0 ? 2 / (nboot + 1) : 2 * pass_obs / nboot

    elseif stattype == "χ²"

        pass_obs == 0 ? 1 / (nboot + 1) : pass_obs / nboot

    end

end

compute_pvalue(stat, boot, stattype) =
    compute_pvalue(pass(stat, boot, stattype), length(boot), stattype)


"""
    beta_approx_pvalue(stat, boot, stattype, q) -> Float64

FastQTL-style beta approximation of the bootstrap p-value, lifting the
empirical `p_boot`'s resolution floor of ~1/B. Both the observed statistic
and every bootstrap draw are transformed to *nominal* p-values under the
analytical reference distribution (two-sided normal for `"z"`, chi-square
with `q` degrees of freedom for `"χ²"`). If that reference were exactly
correct, the bootstrap nominal p-values would be Uniform(0,1) = Beta(1,1);
instead, a `Beta(a, b)` is fit to them by maximum likelihood (with a
method-of-moments fallback), capturing any miscalibration, and the observed
nominal p-value is evaluated under the fitted Beta's CDF. The result is a
smooth p-value that extrapolates below 1/B through the parametric form --
exactly how FastQTL extrapolates beyond its permutation count, applied here
to the score-bootstrap distribution. Returns `NaN` if the fit is degenerate.
"""
function beta_approx_pvalue(stat, boot, stattype, q::Int)

    tonominal(x) = stattype == "z" ? 2 * ccdf(Normal(), abs(x)) : ccdf(Chisq(q), x)

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


"""
    crvetest_direct(R, A, Sg, clustshare)

Direct implementation of the analytical CRVE score test for the shared-null
fast path, replicating `WildBootTests.scoretest`'s ml-mode statistic
(`reps = 0`, `imposenull = true`, `scorebs = true`, `small = false`,
one-way clustering) without the library's per-call setup:

  - `numer = (A R')' S`, with `S` the total score sum;
  - cluster-level `K = Sg * (A R')`, recentered by each cluster's share of
    observations, so `denom = K'K` is the CRVE meat transformed by `A R'`;
  - q = 1: `z = numer / sqrt(denom)`, `p = ccdf(Chisq(1), z²)` (two-tailed);
    q > 1: `χ² = numer' denom⁻¹ numer`, `p = ccdf(Chisq(q), χ²)`.

`Sg` is the G × k matrix of per-cluster score sums; `clustshare[g]` is
cluster g's share of observations. Returns `(stat, p, stattype)` in the same
shape `crvetest` results are consumed.
"""
function crvetest_direct(R::AbstractMatrix, A::AbstractMatrix,
                         Sg::AbstractMatrix, clustshare::AbstractVector)

    ARt = A * R'                          # k × q
    S = vec(sum(Sg, dims = 1))            # total score sum (k)
    numer = ARt' * S                      # q

    K = Sg * ARt                          # G × q
    K .-= clustshare * (S' * ARt)         # recenter, as WildBootTests does for scorebs
    denom = K' * K                        # q × q CRVE denominator

    q = size(R, 1)
    if q == 1
        z = numer[1] / sqrt(denom[1, 1])
        return (stat = z, p = ccdf(Chisq(1), z^2), stattype = "z")
    else
        χ² = dot(numer, Symmetric(denom) \ numer)
        # Per-component 1-df tests of each tested coefficient, from the same
        # restricted fit: numer_k² against the k-th diagonal of the CRVE
        # denominator (see crve_percomponent). Free by-products of the joint
        # statistic; consumed by the per-context interaction columns.
        pk = [ccdf(Chisq(1), numer[k]^2 / denom[k, k]) for k in 1:q]
        return (stat = χ², p = ccdf(Chisq(q), χ²), stattype = "χ²", pk = pk)
    end

end


"""
    crve_percomponent(R, A, Sg, clustshare) -> Vector{Float64}
    crve_percomponent(R, A; scores, clustid) -> Vector{Float64}

Per-component 1-df CRVE score p-values for each tested row of `R`, computed
from the *same* restricted (null) fit as the joint test: with `numer` and
`K` exactly as in [`crvetest_direct`](@ref),

    p_k = ccdf(Chisq(1), numer[k]² / (K'K)[k,k])

These are tests of each coefficient under the joint null (all tested
coefficients zero) -- exact size under that null, but under partial
alternatives with correlated contexts, a true effect in one context can
project onto another's score. They decompose the joint statistic; the joint
test remains the primary inference.

The keyword method builds per-cluster score sums from row-level `scores`
and one-way integer cluster ids (used by the generic, non-workspace path).
"""
function crve_percomponent(R::AbstractMatrix, A::AbstractMatrix,
                           Sg::AbstractMatrix, clustshare::AbstractVector)

    ARt = A * R'
    S = vec(sum(Sg, dims = 1))
    numer = ARt' * S
    K = Sg * ARt
    K .-= clustshare * (S' * ARt)
    return [ccdf(Chisq(1), numer[k]^2 / dot(view(K, :, k), view(K, :, k)))
            for k in 1:size(R, 1)]

end

function crve_percomponent(R::AbstractMatrix, A::AbstractMatrix;
                           scores::AbstractMatrix, clustid::AbstractVector{<:Integer})

    G, n = maximum(clustid), length(clustid)
    Sg = zeros(G, size(scores, 2))
    clustshare = zeros(G)
    @inbounds for i in 1:n
        g = clustid[i]
        clustshare[g] += 1.0
        @views Sg[g, :] .+= scores[i, :]
    end
    clustshare ./= n
    return crve_percomponent(R, A, Sg, clustshare)

end
