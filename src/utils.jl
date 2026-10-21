function printlnv(args...; verbose = true)
    if verbose
        println(args...)
    end
end


function verify_nobs_map_locus(geno::AbstractMatrix, pheno::AbstractVector, 
                                meta::AbstractDataFrame, groups::AbstractDataFrame)

    geno_size = size(geno, 1)
    pheno_size = length(pheno)
    meta_size = nrow(meta)
    groups_size = nrow(groups)



   n_obs = Dict(:variant => geno_size, :expr => pheno_size, 
                 :meta => meta_size, :groups => groups_size)

   n_obs_vals  = collect(values(n_obs))

   if(!all(x -> x == first(n_obs_vals), n_obs_vals))
       
       println("Verify number of observations of input files")
       println("Number of observations per input data\n")
       for (key, value) in pairs(n_obs)
           println(key, " = ", value)
           
       end
       throw("Number of observations in input data do not match between input files")
       
   end


end


# Column names Dynema itself writes to the per-variant summary table (the
# statistic column is named after its type, see `smallsample_factors`).
# A model term (covariate/context) with one of these names would produce a
# coefficient column of the same name, which either errors (duplicate column
# names) or silently overwrites the Dynema column.
const RESERVED_OUTPUT_COLUMNS =
    ("variant", "chr", "pos", "z", "t", "χ²", "F", "p", "p_boot", "p_boot_approx")

"""
    check_reserved_columns(terms, pcnames = nothing)

Stops with an `ArgumentError` naming the offending terms when a model term
(covariate or context column of the formula) would collide with a column
Dynema writes to the summary table: the fixed names in
`RESERVED_OUTPUT_COLUMNS`, the per-context `p_<context>` columns (`pcnames`),
or the `se_wald_<term>` column of another term (a term literally called
`se_wald_x` next to a term `x`). The fix is to rename the column in the
metadata.
"""
function check_reserved_columns(terms::AbstractVector, pcnames::Union{Nothing, AbstractVector} = nothing)
    terms = String.(terms)
    reserved = Set{String}(RESERVED_OUTPUT_COLUMNS)
    pcnames === nothing || union!(reserved, String.(pcnames))
    clash = sort(unique(filter(terms) do t
        t in reserved ||
            (startswith(t, "se_wald_") && t[length("se_wald_")+1:end] in terms)
    end))
    isempty(clash) && return nothing
    throw(ArgumentError(
        "model term(s) " * join(("\"$c\"" for c in clash), ", ") * " have the same name as output " *
        "column(s) Dynema writes (reserved names: " * join(RESERVED_OUTPUT_COLUMNS, ", ") *
        ", p_<context>, se_wald_<term>). Rename the column(s) in the metadata (--meta) and rerun."))
end


# Assume:
# R is a BitVector (or any logical vector)
# betas0 is a vector with length equal to count(!, R) (i.e., number of `false` in R)
# We want a new vector betas_filled of same size as R, where:
#   - betas_filled[i] = betas0[j] if !R[i]
#   - betas_filled[i] = 0       if R[i]

# Ensembl-style ids carry annotation-version suffixes (ENSG00000163131.10) in
# some sources (e.g. GTF-derived annotations) but not others (e.g. 10x
# features files); gene-identifier comparisons ignore them.
stripver(x::AbstractString) = startswith(x, "ENS") ? String(first(split(x, '.'))) : String(x)


function insert_zeros(R::AbstractVector{Bool}, betas0::AbstractVector)
    betas_filled = similar(R, eltype(betas0))  # create output vector of correct type
    idx_betas0 = 1
    for i in eachindex(R)
        if R[i]
            betas_filled[i] = 0
        else
            betas_filled[i] = betas0[idx_betas0]
            idx_betas0 += 1
        end
    end
    return betas_filled
end
