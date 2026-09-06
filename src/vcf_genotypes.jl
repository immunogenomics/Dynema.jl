# ---------------------------------------------------------------------------- #
#                               vcf_genotypes.jl                               #
# ---------------------------------------------------------------------------- #
#
# Core VCF -> donor x variant dosage matrix extraction, usable directly from
# Julia code (not just via bin/dynema-map) -- see
# `extract_geno_dataframe` below, the public entry point.
#
# Assumes the VCF already carries per-sample genotype dosages (FORMAT field
# `DS`) and/or genotype probabilities (FORMAT field `GP`) -- as produced by
# standard imputation pipelines (Minimac, IMPUTE2/5, Beagle, etc.). Hard-call
# genotypes (`GT`) are not read; a variant with neither `DS` nor `GP` is
# skipped. Region seeking is delegated to the `tabix` executable bundled by
# `htslib_jll` (a Julia package artifact -- no system htslib install needed);
# all GP/DS-to-dosage parsing happens here in Julia.

const VCF_SAMPLE_START = 10 # VCF fixed columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT, then samples

"""
    read_gene_bed(bed) -> Vector of (gene, chr, tss, strand)

Reads the bed-like file that tells Dynema *what* to map and *where*: a
plain-text (possibly gzipped) whitespace/tab-separated table with columns
`chr`, `start`, `end`, `gene`, `strand` -- standard 6-column BED (with a
score column before the strand) is also accepted, and `#` comments/header
rows are skipped. Each data row is one gene to map (so a multi-row file
defines a batch): its gene column (a name/symbol or a gene id) names the
gene, and the TSS is derived from its positions the way FastQTL does -- the
`start` position on the plus strand, the `end` position on the minus
strand. Positions are taken as given; note classic BED starts are 0-based
(one bp off a GTF-style start), which is immaterial for cis-window queries.
Duplicated gene names (which would collide on output filenames) and empty
files are errors. Internal helper for `resolve_region` and the CLIs; not
exported.
"""
function read_gene_bed(bed::AbstractString)

    genes = NamedTuple{(:gene, :chr, :tss, :strand), Tuple{String, String, Int, String}}[]

    io = open_maybe_gzip(bed)
    try
        for line in eachline(io)
            l = strip(line)
            (isempty(l) || startswith(l, "#")) && continue
            f = split(l)
            length(f) >= 5 ||
                error("Malformed line in $bed (need at least chr, start, end, gene, strand): '$l'")
            tryparse(Int, f[2]) === nothing && continue  # header row
            strand = length(f) >= 6 && f[6] in ("+", "-", ".") ? String(f[6]) : String(f[5])
            start = parse(Int, f[2])
            stop = parse(Int, f[3])
            tss = strand == "-" ? stop : start
            push!(genes, (gene = String(f[4]), chr = String(f[1]), tss = tss, strand = strand))
        end
    finally
        close(io)
    end

    isempty(genes) &&
        error("No genes found in $bed; each data row must be one gene (chr, start, end, gene, strand)")

    dupes = [g for (g, c) in countmap([g.gene for g in genes]) if c > 1]
    isempty(dupes) ||
        error("Duplicated gene name(s) in $bed ($(join(first(dupes, 5), ", "))); " *
              "each gene must appear once (its name determines the output filename)")

    return genes

end

"""
    resolve_region(; chr, tss, bed, window) -> (chr, tss, start_pos, end_pos)

Resolves the chromosome/TSS -- either from a single-gene bed-like `bed` file
(FastQTL-style: gene start on the + strand, gene end on the - strand; see
[`read_gene_bed`](@ref)), or directly from `chr`/`tss` -- then applies
`window` to get the [start_pos, end_pos] region to query. Internal helper
for `extract_geno_dataframe`; not exported.
"""
function resolve_region(; chr = nothing, tss = nothing, bed = nothing, window::Int)

    if bed !== nothing
        genes = read_gene_bed(bed)
        length(genes) == 1 ||
            error("$bed contains $(length(genes)) genes; this extraction handles one gene at a time -- " *
                  "pass a single-gene file (or chr/tss directly), or use dynema-map for multi-gene batches")
        b = genes[1]
        chr === nothing || string(chr) == b.chr ||
            @warn "chr ($chr) differs from the chromosome in the annotation ($(b.chr)); using $(b.chr) for the VCF query"
        chr = b.chr
        tss = b.tss
    elseif chr === nothing || tss === nothing
        error("Either a single-gene bed-like file, or both chr and tss, must be provided to locate the cis-window")
    end

    start_pos = max(1, tss - window)
    end_pos = tss + window

    return chr, tss, start_pos, end_pos

end

"""
    list_vcf_contigs(vcf) -> Vector{String}

Lists the chromosome/contig names present in `vcf`'s tabix index (`tabix -l`),
without querying any actual variant data. Internal helper for `verify_chr`;
not exported.
"""
function list_vcf_contigs(vcf::AbstractString)

    out = try
        htslib_jll.tabix() do tabix_exe
            read(`$tabix_exe -l $vcf`, String)
        end
    catch e
        error("tabix failed to list chromosomes/contigs in $vcf: $e")
    end

    return split(out, '\n'; keepempty = false)

end

"""
    verify_chr(vcf, chr)

Errors with a clear, actionable message if `chr` isn't one of `vcf`'s
indexed chromosome/contig names. This is the single most common reason a
cis-window query silently returns zero variants -- a chromosome-naming
mismatch (e.g. requesting `chr1` against a VCF indexed with `1`, or vice
versa) -- so it's checked directly rather than left to surface as an
unexplained "no variants found" downstream. Suggests the flipped
`chr`/no-`chr` spelling when it matches, since that's by far the most common
case. Internal helper for `run_tabix`; not exported.
"""
function verify_chr(vcf::AbstractString, chr::AbstractString)

    contigs = list_vcf_contigs(vcf)
    chr in contigs && return nothing

    alt = startswith(chr, "chr") ? chr[4:end] : "chr" * chr
    suggestion = alt in contigs ? " Did you mean '$alt'?" : ""

    shown = join(first(contigs, 10), ", ")
    more = length(contigs) > 10 ? ", ..." : ""

    error("Chromosome '$chr' not found in $vcf's index.$suggestion " *
          "Available chromosome(s) in this VCF: $shown$more. Check the annotation's " *
          "(or chr's) chromosome naming convention against the VCF.")

end

"""
    is_bgzf(path) -> Bool

Whether `path` is BGZF-compressed (the blocked gzip variant htslib's tabix
requires) as opposed to plain gzip or uncompressed: checks the gzip magic
bytes, the FEXTRA flag, and the leading 'BC' extra-subfield identifier that
BGZF blocks carry. Internal helper for `prepare_vcf`; not exported.
"""
function is_bgzf(path::AbstractString)
    open(path) do io
        magic = read(io, 14)
        length(magic) == 14 &&
            magic[1] == 0x1f && magic[2] == 0x8b &&
            (magic[4] & 0x04) != 0 &&
            magic[13] == UInt8('B') && magic[14] == UInt8('C')
    end
end

"""
    prepare_vcf(; vcf, force = false, verbose = true) -> String

One-time preparation of a VCF for Dynema's indexed region queries, using the
htslib binaries bundled with the package (`htslib_jll`) -- no system
`bgzip`/`tabix` install needed. **The input file is never modified**: any
needed compression writes a *new*, Dynema-namespaced file alongside it,
`<name>_dynema.vcf.bgz` (`.bgz` is htslib's extension for BGZF files) --
a name no user file plausibly occupies, so nothing of the user's can be
overwritten. The three states a VCF arrives in:

  - plain `.vcf`: BGZF-compressed to a new `<name>_dynema.vcf.bgz`;
  - `.vcf.gz` that is plain gzip rather than BGZF (tabix rejects these):
    recompressed to a new `<name>_dynema.vcf.bgz`;
  - already BGZF: used as is (no new file).

On a rerun, an existing `<name>_dynema.vcf.bgz` is reused when its
decompressed contents still match the input, and regenerated when they
don't (a stale copy of an input that has since changed). The resulting
file is then tabix-indexed (skipped if an up-to-date `.tbi`/`.csi` already
exists, unless `force = true`; an index older than the file it indexes is
rebuilt). Returns the path of the ready, indexed file -- which is what to
pass to `dynema-map --vcf` (it differs from the input when a compression
step ran).

The VCF must be coordinate-sorted: Dynema deliberately never modifies or
reorders a VCF's records, so an unsorted file errors with instructions to
sort it with external tooling (e.g. `bcftools sort`) first.
"""
function prepare_vcf(; vcf::AbstractString, force::Bool = false, verbose::Bool = true)

    isfile(vcf) || error("VCF file not found: $vcf")

    path = vcf
    rewrote = false  # (re)compressed in this call -> any existing index is stale
    if endswith(vcf, ".gz") && is_bgzf(vcf)
        # already BGZF: used as is, nothing to write
    else
        # Both compression cases write to one predictable, Dynema-namespaced
        # output name -- <name>_dynema.vcf.bgz -- so a Dynema-made file can
        # never be confused with (or collide with) the user's own files.
        stem = endswith(vcf, ".vcf.gz") ? chop(vcf; tail = 7) :
               endswith(vcf, ".vcf")    ? chop(vcf; tail = 4) :
               endswith(vcf, ".gz")     ? chop(vcf; tail = 3) : vcf
        bgz = stem * "_dynema.vcf.bgz"
        endswith(vcf, ".gz") &&
            printlnv("$vcf is gzip- but not BGZF-compressed (tabix requires BGZF)"; verbose)
        rewrote = bgzf_compress(vcf, bgz; force, verbose)
        path = bgz
    end

    idx = isfile(path * ".tbi") ? path * ".tbi" :
          isfile(path * ".csi") ? path * ".csi" : nothing
    stale = idx !== nothing && (rewrote || mtime(idx) < mtime(path))
    if !force && idx !== nothing && !stale
        printlnv("$path is already indexed; nothing to do (use force = true to rebuild)"; verbose)
    else
        stale && printlnv("$idx is out of date (older than the VCF); rebuilding it..."; verbose)
        printlnv("Indexing $path with tabix..."; verbose)
        tabix_index(path)
    end

    return String(path)

end

# BGZF-compress `src` to a NEW file `dest`, streaming through the bundled
# bgzip and swapping in atomically via a temp file. `src` is never touched.
# `dest` carries Dynema's own `_dynema` suffix (see `prepare_vcf`), so a
# file already sitting there is a product of a previous run, never the
# user's primary data: it is reused when its decompressed contents still
# match `src` (see `same_decompressed_prefix`) and regenerated when they
# don't (a stale copy of an input that has since changed). Returns `true`
# when `dest` was (re)written, `false` when an existing copy was reused.
function bgzf_compress(src::AbstractString, dest::AbstractString; force::Bool, verbose::Bool)

    if isfile(dest) && !force
        if is_bgzf(dest) && same_decompressed_prefix(src, dest)
            printlnv("$dest already exists and matches the contents of $src; using it"; verbose)
            return false
        end
        printlnv("$dest exists but no longer matches the contents of $src (stale copy of " *
                 "an earlier input?); regenerating it..."; verbose)
    end

    printlnv("BGZF-compressing $src to $dest (the input file is never modified)..."; verbose)
    tmp = dest * ".tmp"
    try
        htslib_jll.bgzip() do bgzip_exe
            open(tmp, "w") do out_io
                proc = open(pipeline(`$bgzip_exe -c`; stdout = out_io), "w")
                src_io = open_maybe_gzip(src)
                try
                    while !eof(src_io)
                        write(proc, read(src_io, 1 << 20))
                    end
                finally
                    close(src_io)
                    close(proc)
                end
            end
        end
        mv(tmp, dest; force = true)
    catch
        rm(tmp; force = true)
        rethrow()
    end

    return true

end

# Read exactly `n` decompressed bytes from `io` (fewer only at EOF).
function _readn(io::IO, n::Int)
    buf = read(io, n)
    while length(buf) < n && !eof(io)
        append!(buf, read(io, n - length(buf)))
    end
    return buf
end

# Do two (possibly gzip/BGZF-compressed) files agree over the first `n`
# bytes of their *decompressed* contents -- including ending together if
# shorter? Cheap identity check used to recognize a previous run's
# recompression output without decompressing multi-GB files in full: for a
# VCF, the first megabytes cover the header, sample columns, and the first
# variant records, so an unrelated file colliding on all of that is not a
# realistic risk.
function same_decompressed_prefix(a::AbstractString, b::AbstractString; n::Int = 1 << 22)
    ia = open_maybe_gzip(a)
    ib = open_maybe_gzip(b)
    try
        remaining = n
        while remaining > 0
            want = min(remaining, 1 << 16)
            ca = _readn(ia, want)
            cb = _readn(ib, want)
            ca == cb || return false          # includes differing lengths
            length(ca) < want && return true  # both ended together
            remaining -= want
        end
        return true
    finally
        close(ia)
        close(ib)
    end
end

# Index `path` with tabix, translating the common "not coordinate-sorted"
# rejection into an actionable error. Dynema deliberately never reorders or
# rewrites a VCF's records -- sorting is the user's call, with their own
# tooling, so there is no risk of this package damaging primary data.
function tabix_index(path::AbstractString)
    htslib_jll.tabix() do tabix_exe
        errbuf = Pipe()
        proc = run(pipeline(ignorestatus(`$tabix_exe -f -p vcf $path`); stderr = errbuf); wait = false)
        close(errbuf.in)
        msg_task = @async read(errbuf, String)
        wait(proc)
        msg = fetch(msg_task)
        success(proc) && return
        if occursin(r"sorted|out of order|unsorted"i, msg)
            error("$path is not coordinate-sorted, so tabix cannot index it. " *
                  "Dynema never modifies or reorders your VCF, so please sort it yourself first -- e.g. " *
                  "`bcftools sort -Oz -o sorted.vcf.gz $path` -- and run dynema-prepare-vcf on the sorted file.")
        end
        error("tabix failed to index $path: $(strip(msg))")
    end
end

"""
    vcf_samples(vcf) -> Vector{String}

Sample names from `vcf`'s header (the #CHROM line), read via `tabix -H`
without touching any variant data. Internal helper (used by dynema-map's
--check mode); not exported.
"""
function vcf_samples(vcf::AbstractString)
    isfile(vcf) || error("VCF file not found: $vcf")
    out = try
        htslib_jll.tabix() do tabix_exe
            read(`$tabix_exe -H $vcf`, String)
        end
    catch e
        error("tabix failed to read the header of $vcf: $e")
    end
    lines = split(out, '\n'; keepempty = false)
    i = findlast(l -> startswith(l, "#CHROM"), lines)
    i === nothing && error("No #CHROM header line found in $vcf")
    return String.(split(lines[i], '\t'))[VCF_SAMPLE_START:end]
end

"""
    run_tabix(vcf, chr, start_pos, end_pos; verbose=true) -> (header_fields, data_lines)

Uses the `tabix` executable bundled by `htslib_jll` to fetch the VCF header
(for sample names) and the data lines overlapping the region in a single
indexed query. Internal helper for `extract_geno_dataframe`; not exported.
"""
function run_tabix(vcf::AbstractString, chr::AbstractString, start_pos::Int, end_pos::Int; verbose::Bool = true)

    isfile(vcf) || error("VCF file not found: $vcf")
    (isfile(vcf * ".tbi") || isfile(vcf * ".csi")) ||
        error("No tabix index found for $vcf. Run `./bin/dynema-prepare-vcf --vcf $vcf` once (uses the bundled htslib; no system tabix install needed), or `tabix -p vcf $vcf` if you have tabix.")

    verify_chr(vcf, chr)

    region = "$(chr):$(start_pos)-$(end_pos)"
    printlnv("Querying $vcf for region $region with tabix (via htslib_jll)..."; verbose)

    out = try
        htslib_jll.tabix() do tabix_exe
            read(`$tabix_exe -h $vcf $region`, String)
        end
    catch e
        error("tabix failed to query $vcf for region $region: $e")
    end

    lines = split(out, '\n'; keepempty = false)
    header_idx = findfirst(l -> startswith(l, "#CHROM"), lines)
    header_idx === nothing && error("Could not find a #CHROM header line in tabix output for $vcf")
    header_fields = split(lines[header_idx], '\t')
    data_lines = filter(l -> !startswith(l, "#"), lines)

    return header_fields, data_lines

end

"""
    variant_dosage(fields, sample_start, field_pref) -> (variant_id, dosage_or_nothing, reason_or_nothing)

Parses one VCF data line into a variant id and a per-sample dosage vector,
preferring GP or DS per `field_pref` ("auto", "GP", "DS"). Returns
`dosage = nothing` with a skip `reason` if the variant should be dropped
(multiallelic, or missing the requested field). Internal helper for
`extract_geno_dataframe`; not exported.
"""
function variant_dosage(fields::Vector{<:AbstractString}, sample_start::Int, field_pref::AbstractString)

    chrom, pos, id, ref, alt = fields[1], fields[2], fields[3], fields[4], fields[5]
    variant_id = id == "." ? "$(chrom)_$(pos)_$(ref)_$(alt)" : id

    occursin(",", alt) && return variant_id, nothing, "multiallelic"

    format_keys = split(fields[9], ':')
    idx_gp = findfirst(==("GP"), format_keys)
    idx_ds = findfirst(==("DS"), format_keys)

    use_field = if field_pref == "GP"
        idx_gp === nothing ? nothing : :GP
    elseif field_pref == "DS"
        idx_ds === nothing ? nothing : :DS
    else # auto: prefer GP, fall back to DS
        idx_gp !== nothing ? :GP : (idx_ds !== nothing ? :DS : nothing)
    end

    use_field === nothing && return variant_id, nothing, "neither GP nor DS present (FORMAT=$(fields[9]))"

    idx = use_field == :GP ? idx_gp : idx_ds
    sample_fields = @view fields[sample_start:end]
    n = length(sample_fields)
    dosage = Vector{Union{Missing,Float64}}(undef, n)

    for (i, sf) in enumerate(sample_fields)
        subfields = split(sf, ':')
        raw = idx <= length(subfields) ? subfields[idx] : "."
        dosage[i] = if raw == "." || isempty(raw)
            missing
        elseif use_field == :DS
            parse(Float64, raw)
        else # GP: three genotype probabilities (0/0, 0/1, 1/1) -> expected dosage
            p = parse.(Float64, split(raw, ','))
            length(p) == 3 ? (p[2] + 2 * p[3]) : missing
        end
    end

    return variant_id, dosage, nothing

end

"""
    extract_geno_dataframe(; vcf, chr=nothing, tss=nothing, bed=nothing,
                            window=500_000, field="auto", samples_file=nothing,
                            donor_col="donor_id", maf=0.0, max_missing=0.1, verbose=true)

Core VCF cis-window extraction: resolves the region around a TSS (from a
single-gene bed-like `bed` file -- columns chr, start, end, gene, strand --
FastQTL-style: gene start on the + strand, gene end on the - strand; or
given directly via `chr`/`tss`),
queries a bgzipped, tabix-indexed VCF for that region, converts each
variant's `GP`/`DS` FORMAT field to an expected allele dosage (preferring
`GP`, per `field`), optionally remaps VCF sample names to donor ids via
`samples_file` (columns `vcf_id`, `donor_id`), and drops variants that are
multiallelic, lack a usable dosage field, or exceed `max_missing`, or fall
below `maf`. Remaining missing dosages are mean-imputed.

Returns a `NamedTuple` with:
- `geno::DataFrame`: one `donor_col` column plus one column per retained
  variant, ready to hand to `map_locus` (via `expand_genotypes`) directly.
- `positions::Vector{Int}`: each retained variant's genomic position (the VCF
  POS field), aligned with `geno`'s variant columns.
- `chr`, `tss`, `start_pos`, `end_pos`: the resolved cis-window.
- `n_samples_vcf`, `n_samples_matched`: sample counts before/after `samples_file` matching.
- `n_variants_total`, `n_multiallelic`, `n_no_field`, `n_high_missing`, `n_retained`:
  how many variants were found in the region and why any were skipped.

Set `verbose=false` to suppress the built-in progress `println`s -- e.g. for
a caller (like the `dynema-map` CLI) that wants to
render its own progress messages from the returned counts instead.
"""
function extract_geno_dataframe(; vcf::AbstractString,
                                 chr::Union{Nothing,AbstractString} = nothing,
                                 tss::Union{Nothing,Int} = nothing,
                                 bed::Union{Nothing,AbstractString} = nothing,
                                 window::Int = 500_000,
                                 field::AbstractString = "auto",
                                 samples_file::Union{Nothing,AbstractString} = nothing,
                                 donor_col::AbstractString = "donor_id",
                                 maf::Float64 = 0.0,
                                 max_missing::Float64 = 0.1,
                                 verbose::Bool = true)

    chr, tss, start_pos, end_pos = resolve_region(chr = chr, tss = tss, bed = bed, window = window)
    printlnv("Gene TSS = $(chr):$(tss); cis-window = $(start_pos)-$(end_pos) (+/- $(window) bp)"; verbose)

    header_fields, data_lines = run_tabix(vcf, chr, start_pos, end_pos; verbose)
    vcf_samples = String.(header_fields[VCF_SAMPLE_START:end])

    # ------------------------------ Sample mapping ----------------------------- #

    donor_ids = vcf_samples
    keep_idx = collect(1:length(vcf_samples))
    if samples_file !== nothing
        smap = CSV.read(samples_file, DataFrame)
        for c in ("vcf_id", "donor_id")
            c in names(smap) || error("--samples is missing required column '$c'")
        end
        id_to_donor = Dict(zip(smap.vcf_id, smap.donor_id))
        keep_idx = findall(s -> haskey(id_to_donor, s), vcf_samples)
        isempty(keep_idx) && error("None of the VCF sample names matched 'vcf_id' in --samples")
        n_unmatched = length(vcf_samples) - length(keep_idx)
        n_unmatched > 0 && @warn "$n_unmatched VCF sample(s) not found in --samples were dropped"
        donor_ids = [id_to_donor[vcf_samples[i]] for i in keep_idx]
    end
    n_samples = length(keep_idx)
    printlnv("N samples in VCF = $(length(vcf_samples)); N samples retained after sample matching = $n_samples"; verbose)

    out_df = DataFrame()
    out_df[!, donor_col] = donor_ids

    if isempty(data_lines)
        @warn "No variants found in $(chr):$(start_pos)-$(end_pos); returning donor-id-only table"
        return (; geno = out_df, positions = Int[], chr, tss, start_pos, end_pos,
                  n_samples_vcf = length(vcf_samples), n_samples_matched = n_samples,
                  n_variants_total = 0, n_multiallelic = 0, n_no_field = 0,
                  n_high_missing = 0, n_retained = 0)
    end

    # ------------------------------ Parse variants ------------------------------ #

    variant_ids = String[]
    variant_pos = Int[]
    columns = Vector{Float64}[]
    n_multiallelic = 0
    n_no_field = 0
    n_high_missing = 0

    for line in data_lines
        fields = split(line, '\t')
        vid, dosage, reason = variant_dosage(fields, VCF_SAMPLE_START, field)

        if dosage === nothing
            reason == "multiallelic" ? (n_multiallelic += 1) : (n_no_field += 1)
            continue
        end

        dosage = dosage[keep_idx]
        n_missing = count(ismissing, dosage)
        if n_missing / n_samples > max_missing
            n_high_missing += 1
            continue
        end
        if n_missing > 0
            m = mean(skipmissing(dosage))
            dosage = coalesce.(dosage, m)
        end
        dosage = Float64.(dosage)

        maf_i = mean(dosage) / 2
        maf_i = min(maf_i, 1 - maf_i)
        maf_i < maf && continue

        push!(variant_ids, vid)
        push!(variant_pos, parse(Int, fields[2]))
        push!(columns, dosage)
    end

    if verbose
        println("Variants in region: $(length(data_lines))")
        println("  skipped (multiallelic):             $n_multiallelic")
        println("  skipped (missing GP/DS field):       $n_no_field")
        println("  skipped (missingness > $(max_missing)):     $n_high_missing")
        println("  retained after MAF >= $(maf) filter:  $(length(variant_ids))")
    end

    for (vid, col) in zip(variant_ids, columns)
        out_df[!, vid] = col
    end

    return (; geno = out_df, positions = variant_pos, chr, tss, start_pos, end_pos,
              n_samples_vcf = length(vcf_samples), n_samples_matched = n_samples,
              n_variants_total = length(data_lines), n_multiallelic, n_no_field,
              n_high_missing, n_retained = length(variant_ids))

end
