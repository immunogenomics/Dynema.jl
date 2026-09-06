#!/usr/bin/env julia
# ---------------------------------------------------------------------------- #
#                             dynema_prepare_vcf.jl                            #
# ---------------------------------------------------------------------------- #
#
# One-time preparation of a VCF for `dynema-map --vcf`: BGZF-compresses it
# with bgzip if needed and builds the tabix index -- using the htslib
# binaries bundled with Dynema (htslib_jll), so no system bgzip/tabix
# install is required. The input file is NEVER modified: any compression
# writes a new, Dynema-namespaced file alongside it, so nothing of the
# user's can be overwritten. The three states a VCF arrives in:
#
#   - plain .vcf                       -> compressed to a new <name>_dynema.vcf.bgz, indexed
#   - .vcf.gz that is plain gzip       -> recompressed to a new <name>_dynema.vcf.bgz, indexed
#     (tabix rejects non-BGZF gzip)
#   - .vcf.gz already BGZF             -> indexed (no-op if a .tbi/.csi exists)
#
# The VCF must be coordinate-sorted: Dynema deliberately never modifies or
# reorders a VCF's records, so an unsorted file errors with instructions to
# sort it first with external tooling (e.g. bcftools sort). Thin CLI
# wrapper around `Dynema.prepare_vcf` (a core library function, callable
# from any Julia session without this wrapper).
#
# Usage:
#   ./dynema-prepare-vcf --vcf genotypes.vcf.gz
#
# Run with --help for the full list of options.
#
# ---------------------------------------------------------------------------- #
#                        Self-bootstrap on first run                          #
# ---------------------------------------------------------------------------- #
# Shares bin/Project.toml with dynema_map.jl -- if that script has already
# been run once, this one needs no extra setup.

import Pkg

const CLI_DIR = @__DIR__
Pkg.activate(CLI_DIR)

first_run = !isfile(joinpath(CLI_DIR, "Manifest.toml"))
if first_run
    println("First run detected: setting up the shared bin/ environment " *
            "(installing dependencies and, if needed, updating Julia's package " *
            "registry -- this can take several minutes; please wait)...")
    flush(stdout)
end
Pkg.develop(Pkg.PackageSpec(path = joinpath(CLI_DIR, "..")); io = first_run ? stdout : devnull)
Pkg.instantiate()

using ArgParse
using Dynema # prepare_vcf

# ---------------------------------------------------------------------------- #
#                              Argument parsing                                #
# ---------------------------------------------------------------------------- #

function parse_commandline()

    s = ArgParseSettings(
        prog = "dynema_prepare_vcf.jl",
        description = "BGZF-compress (if needed) and tabix-index a coordinate-sorted VCF for dynema-map, using Dynema's bundled htslib -- no system bgzip/tabix required. Run once per VCF.",
    )

    @add_arg_table! s begin
        "--vcf"
            help = "Path to the VCF: plain .vcf (compressed to .vcf.gz alongside; original kept), or .vcf.gz (recompressed in place with bgzip if it is plain gzip rather than BGZF)."
            arg_type = String
            required = true
        "--force"
            help = "Rebuild the tabix index (and redo compression) even if up-to-date outputs already exist. The input file itself is never modified."
            action = :store_true
    end

    return parse_args(s)

end

function main()

    args = parse_commandline()

    t0 = time()
    path = prepare_vcf(vcf = args["vcf"], force = args["force"])
    # When a compression step ran, the ready file is a NEW file (the input is
    # never modified) -- point at it, since that's what --vcf should get.
    println("Done in $(round(time() - t0, digits = 1))s." *
            (path == args["vcf"] ? "" : " Pass this file to dynema-map: --vcf $path"))

end

main()
