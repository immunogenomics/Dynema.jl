using Dynema
using Test
using LinearAlgebra
using Random
using Distributions
using DataFrames
using StatsModels

# Internal helpers under test (not exported)
using Dynema: smallsample_factors, crvetest_direct, crve_percomponent,
              crve_se, pass_counts, combine_pass, compute_pvalue, beta_approx_pvalue,
              check_reserved_columns

# Synthetic score-test ingredients: G clusters, k parameters, tested columns `tested`
function synth(; G = 12, k = 4, nobs = 480, seed = 1)
    rng = MersenneTwister(seed)
    Sg = randn(rng, G, k)
    M = randn(rng, k, k)
    A = inv(M' * M + k * I)
    counts = rand(rng, 1:5, G)
    clustshare = counts ./ sum(counts)
    return A, Sg, clustshare, nobs
end

@testset "Dynema.jl" begin

    @testset "reserved output column names" begin
        # ordinary terms pass
        @test check_reserved_columns(["(Intercept)", "G", "tbru_age", "gPC1", "G & CV1"]) === nothing
        # fixed reserved names, including the t / F statistic types
        for bad in ("p", "z", "t", "χ²", "F", "variant", "p_boot", "pos")
            @test_throws ArgumentError check_reserved_columns(["(Intercept)", "G", bad])
        end
        # per-context p-value columns
        @test check_reserved_columns(["G", "p_CV1x"], ["p_CV1", "p_CV2"]) === nothing
        @test_throws ArgumentError check_reserved_columns(["G", "p_CV1"], ["p_CV1", "p_CV2"])
        # se_wald_<term> of another term
        @test_throws ArgumentError check_reserved_columns(["G", "age", "se_wald_age"])
        @test check_reserved_columns(["G", "se_wald_age"]) === nothing   # no term `age` -> no clash
        # the message names the offenders
        err = nothing
        try
            check_reserved_columns(["G", "p", "pos"])
        catch e
            err = e
        end
        @test err isa ArgumentError
        @test occursin("\"p\"", err.msg) && occursin("\"pos\"", err.msg)

        # end to end: a covariate called `pos` is rejected up front
        rng = MersenneTwister(2)
        G, per = 6, 10; n = G * per
        didx = repeat(1:G, inner = per)
        geno = Float64.(rand(rng, 0:2, G, 2)[didx, :])
        y = Float64[rand(rng, Poisson(2.0)) for _ in 1:n]
        meta = DataFrame(pos = randn(rng, n))
        @test_throws ArgumentError map_locus(@formula(0 ~ 1 + G + pos); pheno = y, geno = geno,
                                             meta = meta, groups = string.("d", didx), termtest = "G")
    end

    @testset "smallsample_factors" begin
        @test smallsample_factors(false, 12, 480) == (1.0, 1.0, 0.0)
        cm, ss, dofr = smallsample_factors(true, 12, 480)
        @test cm ≈ 12 / 11
        @test ss ≈ 480 / 479
        @test dofr == 11.0
        @test_throws ArgumentError smallsample_factors(true, 1, 480)
    end

    @testset "crvetest_direct: z/t, χ²/F and SE consistency" begin
        A, Sg, share, nobs = synth()
        G = size(Sg, 1)

        # 1 df ------------------------------------------------------------
        R1 = BitMatrix([0 0 0 1])
        d0 = crvetest_direct(R1, A, Sg, share)                       # small = false
        d1 = crvetest_direct(R1, A, Sg, share; small = true, nobs = nobs)

        @test d0.stattype == "z" && d1.stattype == "t"
        # Statistic is exactly beta / se in both modes
        @test d0.stat ≈ d0.beta[1] / d0.se[1]
        @test d1.stat ≈ d1.beta[1] / d1.se[1]
        # beta (one-step estimate) does not depend on `small`; SE does
        @test d0.beta ≈ d1.beta
        cm, ss, _ = smallsample_factors(true, G, nobs)
        @test d1.se[1] ≈ d0.se[1] * sqrt(cm / ss)
        # Reference distributions
        @test d0.p ≈ 2 * ccdf(Normal(), abs(d0.stat))
        @test d1.p ≈ 2 * ccdf(TDist(G - 1), abs(d1.stat))
        @test d1.p >= d0.p   # t(G-1) + inflated variance is never more significant

        # > 1 df ------------------------------------------------------------
        R2 = BitMatrix([0 0 1 0; 0 0 0 1])
        j0 = crvetest_direct(R2, A, Sg, share)
        j1 = crvetest_direct(R2, A, Sg, share; small = true, nobs = nobs)
        @test j0.stattype == "χ²" && j1.stattype == "F"
        @test j0.p ≈ ccdf(Chisq(2), j0.stat)
        @test j1.p ≈ ccdf(FDist(2, G - 1), j1.stat)
        # per-component p-values use the same 1-df reference as the 1-df test
        @test j1.pk ≈ crve_percomponent(R2, A, Sg, share; small = true, nobs = nobs)
        @test j0.pk ≈ crve_percomponent(R2, A, Sg, share)
        @test all(j1.pk .>= j0.pk)
        @test j1.pk[2] ≈ 2 * ccdf(TDist(G - 1), abs(j1.beta[2] / j1.se[2]))
    end

    @testset "crve_percomponent keyword method (row-level scores)" begin
        rng = MersenneTwister(7)
        G, per, k = 9, 25, 3
        clustid = repeat(1:G, inner = per)
        n = length(clustid)
        scores = randn(rng, n, k)
        M = randn(rng, k, k); A = inv(M' * M + k * I)
        R = BitMatrix([0 0 1])

        Sg = zeros(G, k); share = zeros(G)
        for i in 1:n
            Sg[clustid[i], :] .+= scores[i, :]; share[clustid[i]] += 1
        end
        share ./= n

        for small in (false, true)
            @test crve_percomponent(R, A; scores = scores, clustid = clustid, small = small) ≈
                  crve_percomponent(R, A, Sg, share; small = small, nobs = n)
        end
    end

    @testset "crve_se (unrestricted Wald-type CRVE SE)" begin
        rng = MersenneTwister(3)
        G, per = 10, 30
        clustid = repeat(1:G, inner = per); n = length(clustid)
        X = hcat(ones(n), randn(rng, n), randn(rng, n))
        β = [0.4, 0.1, -0.2]
        y = Float64[rand(rng, Poisson(m)) for m in exp.(X * β)]
        groups = reshape(clustid, :, 1)

        # Naive per-cluster sandwich
        μ = exp.(X * β)
        bread = inv(X' * (μ .* X))
        meat = zeros(3, 3)
        for g in 1:G
            idx = clustid .== g
            s = X[idx, :]' * (y[idx] .- μ[idx])
            meat .+= s * s'
        end
        V = bread * meat * bread

        @test crve_se(X, y, β, groups) ≈ sqrt.(diag(V))
        cm, ss, _ = smallsample_factors(true, G, n)
        @test crve_se(X, y, β, groups; small = true) ≈ sqrt.(diag(V) .* (cm / ss))
        # multiway clustering is not implemented -> NaN
        @test all(isnan, crve_se(X, y, β, hcat(clustid, clustid)))
    end

    @testset "bootstrap helpers: t/F behave like z/χ²" begin
        boot = [-2.0, -1.0, -0.5, 0.2, 0.9, 1.4, 2.5]
        # signed (z, t): two-sided, min of tail counts
        for st in ("z", "t")
            c1, c2 = pass_counts(1.0, boot, st)
            @test (c1, c2) == (count(<(1.0), boot), count(>(1.0), boot))
            @test combine_pass(c1, c2, st) == min(c1, c2)
            @test compute_pvalue(3, 100, st) ≈ 2 * 3 / 100
            @test compute_pvalue(0, 100, st) ≈ 2 / 101
        end
        # upper-tailed (χ², F)
        for st in ("χ²", "F")
            c1, c2 = pass_counts(1.0, boot, st)
            @test (c1, c2) == (count(>=(1.0), boot), 0)
            @test combine_pass(c1, c2, st) == c1
            @test compute_pvalue(3, 100, st) ≈ 3 / 100
            @test compute_pvalue(0, 100, st) ≈ 1 / 101
        end
        @test_throws ArgumentError compute_pvalue(3, 100, "bogus")

        # beta approximation: reference distribution changes with the type
        rng = MersenneTwister(5)
        zboot = randn(rng, 4000)
        pz = beta_approx_pvalue(2.0, zboot, "z", 1)
        pt = beta_approx_pvalue(2.0, zboot, "t", 1; dof_r = 5)
        @test 0 < pz < 0.2 && 0 < pt < 0.2
        @test pz != pt
        fboot = rand(rng, FDist(2, 8), 4000)
        pF = beta_approx_pvalue(6.0, fboot, "F", 2; dof_r = 8)
        @test 0 < pF < 0.2
        @test_throws ArgumentError beta_approx_pvalue(1.0, zboot, "bogus", 1)
    end

    @testset "map_locus end to end: small switches z/t and adds SEs" begin
        rng = MersenneTwister(11)
        G, per, nv = 14, 30, 3
        didx = repeat(1:G, inner = per); n = length(didx)
        donor = string.("d", didx)
        geno = Float64.(rand(rng, 0:2, G, nv)[didx, :])
        ctx = randn(rng, n)
        μ = exp.(0.3 .+ 0.15 .* geno[:, 1] .+ 0.1 .* ctx .* geno[:, 1] .+ (0.2 .* randn(rng, G))[didx])
        y = Float64[rand(rng, Poisson(m)) for m in μ]
        meta = DataFrame(C = ctx)
        cmss = (G / (G - 1)) / (n / (n - 1))

        # ---- main effect (shared null), 1 df ------------------------------
        run(small) = map_locus(@formula(0 ~ 1 + G); pheno = y, geno = geno, meta = meta,
                               groups = donor, termtest = "G", betas = :all, small = small)
        s0 = get_summary(run(false)); s1 = get_summary(run(true))
        @test get_stattype(run(false)) == "z" && get_stattype(run(true)) == "t"
        @test "z" in names(s0) && "t" in names(s1) && !("t" in names(s0))
        @test !("beta" in names(s0)) && !("se" in names(s0))   # no one-step outputs
        @test s1.p ≈ 2 .* ccdf.(TDist(G - 1), abs.(s1.t))
        @test all(s1.p .>= s0.p)
        @test all(isfinite, s0.se_wald_G) && all(s0.se_wald_G .> 0)
        @test !("se_wald_(Intercept)" in names(s0))   # SEs for genotype terms only
        @test s1.se_wald_G ≈ s0.se_wald_G .* sqrt(cmss)
        @test s0.G ≈ s1.G                       # unrestricted estimates unaffected

        # ---- single-context interaction (per-variant null), 1 df -----------
        runi(small) = map_locus(@formula(0 ~ 1 + G + C + G & C); pheno = y, geno = geno,
                                meta = meta, groups = donor, termtest = "G & C",
                                betas = :all, small = small)
        i0 = get_summary(runi(false)); i1 = get_summary(runi(true))
        @test !("beta" in names(i0)) && !("se" in names(i1))
        @test all(isfinite, i1[!, "se_wald_G & C"])

        # ---- joint (2 df) test: χ² vs F, no single beta / se ----
        runj(small) = map_locus(@formula(0 ~ 1 + G + C + G & C); pheno = y, geno = geno,
                                meta = meta, groups = donor, termtest = ["G", "G & C"],
                                betas = :none, small = small)
        j0 = get_summary(runj(false)); j1 = get_summary(runj(true))
        @test "χ²" in names(j0) && "F" in names(j1)
        @test !("beta" in names(j0))
        @test j1.p ≈ ccdf.(FDist(2, G - 1), j1.F)
    end

    # Private TBRU regression test: local-only
    let f = joinpath(@__DIR__, "tbru_mini_tests.jl")
        isfile(f) && include(f)
    end

end
