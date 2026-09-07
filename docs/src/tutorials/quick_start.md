# Quick start (end to end)

This walkthrough goes from raw files to mapped single-cell eQTL effects (main, interaction, or total), using only the command line. It runs on a small simulated demo dataset (download below) -- 300 donors, 4,050 cells, and 3 genes (*CTSS*, *ACTB*, and *TSPAN32*) with 500 cis-variants each.

For your own study, the same four inputs are needed -- just swap the file
names:

1. **Genotypes**: a coordinate-sorted VCF with imputed dosages (`DS`) and/or
   genotype probabilities (`GP`) -- here `genotypes.vcf.gz`.
2. **Expression**: a Market Exchange Format (MEX) count-matrix (as exported by
   Seurat/scanpy/10x pipelines) -- here `expr.mtx.gz` + `expr.features.gz` +
   `expr.barcodes.gz`.
3. **Metadata**: one row per cell with the cell id, donor id,
   contexts (e.g. cell states), and covariates -- here `meta.tsv`, with donor column
   `donor_id`, cell-id column `cell_id`, contexts
   `cytotoxicity,treg_activation,central_memory`, and normalized/scaled covariates
   `scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1-5,ePC1-5`.
4. **Gene coordinates**: chromosome, start, end, and strand of the variants
   to test -- here `genes.bed` (normally derived from GTF annotation).

## Step 0: Get Julia, Dynema and demo data

- Install [Julia >= 1.9](https://github.com/julialang/juliaup). Super easy!

!!! note
    To use *Dynema*, no Julia coding is needed! Everything is command-line and simple!

- Download *Dynema*:

```bash
curl -L https://github.com/immunogenomics/Dynema.jl/archive/refs/heads/main.tar.gz | tar -xz
mv Dynema.jl-main Dynema
cd Dynema
```

- Download the demo dataset (~1 MB) from [Dynema_datasets](https://github.com/immunogenomics/Dynema_datasets):

```bash
input=quick_start_data
mkdir -p "$input"

curl -L https://github.com/immunogenomics/Dynema_datasets/archive/refs/heads/main.tar.gz | \
  tar -xz -C "$input" --strip-components=3 "Dynema_datasets-main/data/quick_start"
```

!!! note
    No need to worry about installing *Dynema* either! Once you use it on the command line, *Dynema* and all dependencies are automatically installed. This step happens only once and takes just a few minutes. 


## Step 1: Index the VCF (once per VCF)

`dynema-map` extracts each gene's variants to test (e.g. cis-window) from the VCF from a VCF file (BGZF-compressed and tabix-indexed). We provide the tool `dynema-prepare-vcf` in case your input file is not already BGZF-compressed and/or tabix-indexed. 


`dynema-prepare-vcf` automatically BGZF-compresses (if required) and indexes the VCF file:

```bash
./bin/dynema-prepare-vcf --vcf "$input/genotypes.vcf.gz"
```

Your input file is **never modified**: any needed compression writes a new
file alongside it, and the command reports the ready file to use. The demo
VCF is plain gzip (not BGZF), so here that ready file is
`$input/genotypes_dynema.vcf.bgz` -- the file all the commands below pass to
`--vcf`.

!!! note
    This uses Dynema's bundled `htslib` (no system `bgzip`/`tabix` needed!) and
    handles plain `.vcf`, plain-gzip `.vcf.gz`, and already-BGZF files alike --
    see [the prepare-vcf extra](prepare_vcf.md). If your VCF already ships with
    a `.tbi` index, skip this step and pass it to `--vcf` directly.

## Step 2: Index the expression matrix (once per matrix)

Extracting the expression of a single gene normally requires scanning the whole
Matrix Market file. Building a gene index `.dgx` once turns that
into a milliseconds-scale read for every later run:

```bash
./bin/dynema-prepare-expr --expr-prefix "$input/expr"
```

!!! note 
    On a real study-scale matrix (thousands of genes, hundreds of
    thousands of cells, a multi-GB `.mtx.gz`) the benefit of indexing becomes obvious (each gene loading in ~0.5 s instead of minutes).

## Step 3: Specify variants to test

A small bed-like file tells `dynema-map` *what* to map and *where*: one row
per gene with `chr`, `start`, `end`, `gene`, `strand`. The TSS is derived
FastQTL-style: the start position for `+` strand, and the end position for the`-` strand. The demo ships three genes as examples:

```bash
cat "$input/genes.bed"
```

```
chr1	150730706	150765778	CTSS	-
chr7	5527151	5530601	ACTB	-
chr11	2323216	2339430	TSPAN32	+
```

For your own study, write the file for a batch of genes you want to analyze. Preparing bed files as batches facilitates parallelization in HPC infrastructures.

!!! note
    We also provide a tool to build bed files straight from your GTF annotation with `dynema-prepare-bed` -- see [the prepare-bed extra](prepare_bed.md), facilitating parallelization.


## Step 4: Map single-cell eQTL effects

**Main effect**: A convential eQTL test (context-independent). We are trying to answer the question: *is there a constant/average effect of a variant on regulating gene expression*?

The cell-state contexts enter as ordinary covariates here:

```bash
./bin/dynema-map \
  --expr-prefix "$input/expr" \
  --meta "$input/meta.tsv" \
  --vcf "$input/genotypes_dynema.vcf.bgz" \
  --bed "$input/genes.bed" \
  --covariates scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1,gPC2,gPC3,gPC4,gPC5,ePC1,ePC2,ePC3,ePC4,ePC5,cytotoxicity,treg_activation,central_memory \
  --interaction-with cytotoxicity,treg_activation,central_memory \
  --donor-col donor_id --cell-id-col cell_id \
  --effect main \
  --out "$input/main"
```

---

**Multi-context Interaction effect**: This test leverages the heterogeneity in single-cell data to infer dynamic effect of variants on gene expression. We are trying to answer the question: *does the variant's effect on gene expression **change** across the contexts* (e.g. `cytotoxicity,treg_activation,central_memory`)?

The context's main effects are added to the model automatically via `--interaction-with`:

```bash
./bin/dynema-map \
  --expr-prefix "$input/expr" \
  --meta "$input/meta.tsv" \
  --vcf "$input/genotypes_dynema.vcf.bgz" \
  --bed "$input/genes.bed" \
  --covariates scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1,gPC2,gPC3,gPC4,gPC5,ePC1,ePC2,ePC3,ePC4,ePC5 \
  --interaction-with cytotoxicity,treg_activation,central_memory \
  --donor-col donor_id --cell-id-col cell_id \
  --effect interaction \
  --out interaction
```


!!! note
    Alongside the joint (multi-context) p-value `p`, Dynema automatically
    returns the **single-context interaction p-values** as one extra column per
    context (here `p_cytotoxicity`, `p_treg_activation`, `p_central_memory`) --
    computed from the same model fit, at no extra cost. These show *which*
    context(s) drive a joint signal; the joint `p` remains the primary test
    (pass `--per-context false` to omit them).

---

**Total effect**: this test captures both main effect and its interactions jointly. We are trying to answer the question: *is there **any** genetic effect, constant or context-dependent regulating gene expression?*:

```bash
./bin/dynema-map \
  --expr-prefix "$input/expr" \
  --meta "$input/meta.tsv" \
  --vcf "$input/genotypes_dynema.vcf.bgz" \
  --bed "$input/genes.bed" \
  --covariates scaled_age,sex,scaled_log_nUMI,percent_mito,gPC1,gPC2,gPC3,gPC4,gPC5,ePC1,ePC2,ePC3,ePC4,ePC5 \
  --interaction-with cytotoxicity,treg_activation,central_memory \
  --donor-col donor_id --cell-id-col cell_id \
  --effect total \
  --out "$input/total"
```


!!! note
    You can also run multiple tests with a single command (e.g.
    `--effect main,interaction`). Each gene's data is then loaded and
    extracted only once, shared across the tests. The effect name is
    appended to each output file -- `<out>_<effect>_<gene>.tsv`, e.g.
    `results_main_CTSS.tsv` and `results_interaction_CTSS.tsv` with
    `--out results` -- and the summary (one row per gene × effect) stays in
    a single `<out>_summary.tsv`. Note that the `main` test in this mode
    omits the untested interaction terms from its model, exactly matching
    the standalone command above.

## Step 6: Read the results

Each command-line run writes one file per gene plus a summary:

- **`<out>_<gene>.tsv`** (e.g. `main_CTSS.tsv`) -- one row per cis variant:
  `variant`, `chr`, `pos` (from the VCF), the score statistic (`z`, or `χ²`
  for joint tests), the analytical `p`, and effect-size estimates for the
  lead variant(s) (the `--betas lead` default; `--betas all` fits every
  variant).
- **`<out>_summary.tsv`** -- A summary of the lead variants for all tests across all genes.
- **`<out>.log`** -- the full log messages.


Because the data are simulated with known effects, the results should tell a clean story: *TSPAN32* has a hit near its TSS in all three runs (a main effect plus a `G × cytotoxicity` interaction); *CTSS*  has no **main** effect but has considerable **interaction** and **total** effects (its simulated effect exists only through `treg_activation`); *ACTB* stays null everywhere.

## Where to next

- Map thousands of genes by chunking the bed file and submitting one job per
  chunk -- see [batching](command_line.md#Mapping-many-genes:-batches),
  including `--check` and `--skip-existing` for safe HPC runs.
- Add bootstrap *p*-values for small cohorts with `--boot`.
- Every option: `./bin/dynema-map --help` and the
  [CLI overview](command_line.md).
