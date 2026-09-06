# Command-line interface (CLI)

**You don't need to write any Julia to run *Dynema*!**. The `bin/` folder of the
[Dynema.jl repository](https://github.com/joseah/Dynema.jl) ships a
self-contained command-line tool:

- **`dynema-map`**: This tool maps one or more genes (a batch, defined by
  the rows of `--bed`) against their cis genetic variants -- reads gene
  expression, single-cell metadata, and genotypes, runs [`map_locus`](@ref)
  per gene, and writes one summary statistics table per gene.
- **`dynema-prepare-expr`** (optional, run once per expression matrix):
  builds a gene-major index (`.dgx`) next to a Matrix Market export, after
  which `dynema-map` loads any single gene in milliseconds instead of
  scanning the whole multi-GB matrix per run.
- **`dynema-prepare-vcf`** (optional, run once per VCF): BGZF-compresses
  and tabix-indexes a coordinate-sorted VCF using Dynema's bundled htslib --
  for anyone whose VCF doesn't come with a `.tbi` index (or is plain gzip
  rather than BGZF). No system bgzip/tabix install needed.
- **`dynema-prepare-bed`** (optional, run once per study): builds
  `dynema-map`'s bed-like gene file(s) from a GTF annotation -- keeping only
  genes present in the expression data's features file, picking an
  unambiguous identifier per gene, and optionally splitting into fixed-size
  chunks for HPC batching.

These scripts bootstrap their own Julia environment automatically the first
time one is run (installing `Dynema` plus a handful of CLI-only
dependencies), so all you need installed is Julia (>= 1.9)
itself. 


```bash
curl -L https://github.com/joseah/Dynema.jl/archive/refs/heads/main.tar.gz | tar -xz
mv Dynema.jl-main Dynema.jl
cd Dynema.jl
./bin/dynema-map --help
```

If `julia` isn't on your `PATH`, invoke the underlying script directly
instead: `julia --project=bin bin/dynema_map.jl [options]`.

## What `dynema-map` needs

`dynema-map` inputs:

- **Gene expression** (Matrix Market, `--expr-prefix`): the standard sparse
  export (`PREFIX.mtx`, `PREFIX.features`, `PREFIX.barcodes` -- as produced by
  Seurat/scanpy/10x pipelines), containing raw counts. Run `dynema-prepare-expr`
  on it once to build a `.dgx` index, and every later `dynema-map` run loads its
  gene in milliseconds. The matrix is never loaded in full, whatever its size.
- **Genotypes** (VCF, `--vcf`): a bgzipped, tabix-indexed VCF with genotype
  dosages, from which the tested gene's cis-window variants are extracted in
  about a second. Not compressed/indexed yet? Run `dynema-prepare-vcf` on it
  once.
- **The gene(s) to map**: a bed-like plain text file (`--bed`) with columns
  chr, start, end, gene, strand -- one data row per gene. It serves two
  purposes at once: its gene column (a name/symbol or gene id) specifies
  *what* to map, and its positions/strand give each TSS, derived
  FastQTL-style -- the gene's start position on the plus strand, its end
  position on the minus strand. A multi-row file defines a batch: the genes
  are mapped one after another in the same run, each writing its own output
  file. (A small curated table, not a full GTF -- reading it stays instant.)
- **Metadata**: plain text file with cell id, donor id, cell-state
  contexts, and any donor or single-cell covariates, one row per cell.
- **Type of eQTL effect**: either main, interaction, or total.

These two indexed, random-access formats are the only expression/genotype
inputs Dynema accepts -- plain-text expression or genotype tables don't
scale (they must be parsed in full, on every run, using memory proportional
to the whole file). With Matrix Market + VCF, per-gene runs pay seconds of
I/O, not minutes, at any study size.

## Running Dynema

Single-cell eQTL studies start from a single genome-wide,
bgzipped and tabix-indexed VCF (`.vcf.gz` + `.vcf.gz.tbi`/`.csi`) and from
a Matrix Market
count export (`.mtx[.gz]` + `.features[.gz]` + `.barcodes[.gz]`). `dynema-map` 
reads both directly. 

Here's a simple example just to illustrate the parameter usage (do not run):

```bash
./bin/dynema-map \
  --expr-prefix "$input/expr" \
  --meta "$input/meta.tsv" \
  --vcf "$input/genotypes.vcf.gz" \
  --bed "$input/genes.bed" \
  --window 500000 \
  --covariates scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1,gPC2,gPC3,gPC4,gPC5,ePC1,ePC2,ePC3,ePC4,ePC5 \
  --interaction-with cytotoxicity,treg_activation,central_memory \
  --donor-col donor_id \
  --cell-id-col cell_id \
  --effect interaction \
  --out interaction
```

Arguments:

- `--expr-prefix` specifies the basename of the matrix market files (e.g. [expr].mtx.gz, 
  [expr].barcodes.gz, [expr].features.gz). If a [expr].dgx index built by
  `dynema-prepare-expr` sits next to them, the gene loads from it in
  milliseconds; otherwise the matrix is scanned once per run (minutes for
  multi-GB files -- index it!).
- `--vcf`: VCF file (*.vcf or *.vcf.gx) with **genotype dosages, not hard calls.** 
  The VCF must already carry per-sample genotype dosages (`DS`) and/or genotype 
  probabilities (`GP`) -- as produced by standard imputation pipelines (Minimac, IMPUTE2/5,
  Beagle, etc.). `--field auto` (the default) prefers `GP` over `DS` when a
  variant has both; hard-call genotypes (`GT`) are never read. 
  Must be BGZF-compressed and accompanied by a tabix index file (*.tbi) --
  run `dynema-prepare-vcf` once if yours isn't.
- `--bed`: bed-like file (plain or gzipped) specifying the gene(s) to map:
  one data row per gene with columns chr, start, end, gene, strand (a
  standard 6-column BED with a score column also works; header/`#` lines are
  skipped). Each row's gene column -- a gene name/symbol (`CTSS`) or a gene
  id (`ENSG00000163131`) -- names that gene: with a single-column features file
  (e.g. a Seurat export) it must match that column exactly; with a 10x
  features file (gene_id, gene_name, modality) it is searched against
  gene_name (column 2) first and, failing that, against gene_id (column 1) --
  no match is an error (Ensembl id version suffixes are ignored). The TSS is
  derived the same way FastQTL does: the gene's start position on the plus
  strand, its end position on the minus strand. The annotation's chromosome
  naming must match the VCF's (e.g. `chr1` vs `1`). A cis-window is built
  around the derived TSS.
- `--window`: *cis* region half-width in bp around the TSS. By default 500000 (0.5 Mb).
- `--meta`: plain text file with cell id, donor id, cell-state contexts, and any donor or 
  single-cell covariates, one row per cell.
- `covariates`: list of covariates. Must be comma-separated and should match column names in 
  metadata file.
- `donor-col`: Column name in metadata file including the donor ids. This must match the VCF donor ids too.
- `cell-col`: Column name in metadata file including the cell ids/barcodes. This must match the barcodes 
  provided in the expression data.
- `--effect`: Type(s) of single-cell eQTL to test, comma-separated: any of
  `main`, `interaction`, and `total` (e.g. `--effect main,interaction`).
  Multiple effects share each gene's data extraction and write separate
  files, `<out>_<effect>_<gene>.tsv`. In a multi-effect run, `main` uses the
  classic model without G × context terms.
- `--interaction-with`: comma-separated context column(s) in the metadata file to
  test G × context interactions for; their main effects are added to the model
  automatically. Required when `--effect` includes interaction/total. Contexts
  you only want to adjust for (without an interaction) belong in `--covariates`.
- `--out`: output *prefix*. Each gene writes its own summary statistics
  table named `<out>_<gene>.tsv` (or `<out><gene>.tsv` if the prefix ends in
  `/`; directories are created as needed; with no `--out`, just
  `<gene>.tsv`). This makes it easy to distinguish analyses of the same
  genes: `--out main` and `--out interaction` give `main_CTSS.tsv` and
  `interaction_CTSS.tsv`. With `--vcf`, every output automatically includes
  `chr`/`pos` columns (taken from the VCF), and each invocation additionally
  writes `<out>_summary.tsv` with one lead-variant row per gene × effect.
- `--betas`: which variants get effect-size estimates (unrestricted-model
  coefficients) attached as extra output columns -- the test itself never
  needs them. `lead` (the default) fits only the lead variant(s) (smallest
  p-value, including exact ties; other rows are left empty), `all` fits every
  variant (substantially slower), `none` skips them.
- `--boot`: also compute bootstrap p-values via adaptive score bootstrapping
  (recommended for small or imbalanced cohorts). Adds `p_boot` (empirical,
  floored at ~2/B) and `p_boot_approx` (a FastQTL-style beta approximation of
  the bootstrap distribution that extrapolates smoothly below that floor).
  Uses the optional WildBootTests package, installed automatically into the
  `bin/` environment on first use.
- `--workers N`: start N local worker processes and distribute each gene's
  variants across them (implies `--parallel`). Worth it for long runs
  (interaction effects, bootstrapping); short runs are dominated by fixed
  startup costs.


## Mapping many genes: batches

Dynema maps genome-wide studies the same way FastQTL does: split the gene
list into chunks and let your HPC scheduler parallelize over chunks. Each
`dynema-map` invocation takes one bed-like chunk file and maps its genes
sequentially -- metadata, the model formula, and (via
the `.dgx` index) the expression matrix are all loaded once per invocation,
so per-gene overhead stays low. A gene that fails (e.g. absent from the
features file) is reported and skipped without sinking the rest of the
batch, and each gene writes its own `<out>_<gene>.tsv`. For example, with a
genome-wide study:

```bash
# once: build 100-gene chunk beds from your GTF, matched to your features file
./bin/dynema-prepare-bed --gtf gencode.gtf.gz --features expr.features.gz \
  --chunk-size 100 --out beds/chunk

# then submit one job per chunk, e.g. (SLURM):
#   ./bin/dynema-map --bed beds/chunk_001.bed --expr-prefix expr --vcf genotypes.vcf.gz \
#     ... --out main --log main_chunk_001.log
```

See [the prepare-bed extra](prepare_bed.md) for how identifiers are chosen
and which GTF genes are kept.

Each invocation also writes `<out>_summary.tsv` -- one row per gene × effect
with its lead variant and statistics -- so concatenating the chunk summaries
gives the study-wide top-associations table without parsing the per-gene
files.

Two flags make large batches safer: run once with `--check` before
submitting -- it validates everything in seconds (files parse, metadata
columns exist, every gene is found in the expression data, chromosomes match
the VCF index, all metadata donors have genotypes) and exits without
mapping; and submit jobs with `--skip-existing` so a killed or partially
completed job can simply be resubmitted -- finished gene × effect outputs
are skipped (their lead statistics re-read into the summary) and only the
remaining work runs.

## Learn more

Run `./bin/dynema-map --help` for the complete list of options (variant
filtering with `--variants`, VCF dosage-field control with `--field`, sample-id
remapping with `--samples`, MAF/missingness filters, bootstrap schedules with
`--B`, and more).


## Examples: main, interaction, and total effects

The [Quick start (end to end)](quick_start.md) walks through all three of
`--effect`'s modes on a small simulated demo dataset, one command each.

