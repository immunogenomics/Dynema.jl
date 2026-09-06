# Extra: Building gene bed files from a GTF: `dynema-prepare-bed`

`dynema-map --bed` is driven by a small bed-like gene file -- one row per
gene with `chr`, `start`, `end`, `gene`, `strand` (tab-separated), from
which the TSS is derived FastQTL-style (start on the `+` strand, end on the
`-` strand). For a handful of genes it is easiest to write by hand:

```
chr1	150730706	150765778	CTSS	-
chr11	2323216	2339430	TSPAN32	+
```

Gene identifiers are matched against the expression features file by gene
name first, then gene id -- so `CTSS` and `ENSG00000163131` both work, and
Ensembl `.N` version suffixes are ignored.

For a genome-wide run, build the file(s) straight from your GTF annotation
(GENCODE, Ensembl, or the Cell Ranger reference's `genes.gtf`; plain or
gzipped), matched against the expression data's features file:

```bash
./bin/dynema-prepare-bed --gtf genes.gtf.gz --features expr.features.gz \
  --chunk-size 100 --out beds/chunk
```

With `--chunk-size 100` this writes `beds/chunk_001.bed`,
`beds/chunk_002.bed`, ... of at most 100 genes each, in GTF order (so
chunks are genomically contiguous) -- the intended genome-wide pattern is
then **one `dynema-map` job per chunk file** on your HPC scheduler, with
`--check` and `--skip-existing` for safe submission and resubmission (see
[batching](command_line.md#Mapping-many-genes:-batches)). Omit
`--chunk-size` (or pass `0`) to write everything to a single `<out>.bed`.

For each GTF gene record the tool picks the identifier that `dynema-map`
will resolve unambiguously at map time:

- the **gene name/symbol**, when it is unique both within the GTF and in
  the features file (dynema-map searches feature gene names first);
- the **gene id** otherwise (symbol absent, or duplicated in the features
  file or annotated at multiple loci in the GTF);
- genes **not present in the features file at all are dropped** -- they
  could never be mapped anyway.

!!! note
    Passing `--features` is strongly recommended: without it every GTF gene
    record is written using its gene name, with no guarantee it resolves
    against your expression data. `--feature-type` (default `gene`) selects
    which GTF records are read -- the default skips transcript, exon, and
    CDS lines.

!!! warning "Use a gene-level annotation"
    UCSC-style GTFs (e.g. `hg38.knownGene.gtf.gz`) do not work: they have no
    gene-level records and no `gene_name` attribute -- their `gene_id` is
    actually a *transcript* id (`ENST...`) -- so nothing can be matched to
    an expression features file. Use a GENCODE or Ensembl annotation
    (`gencode.vXX.annotation.gtf.gz`) or the Cell Ranger reference's
    `genes.gtf`, which is what your expression data was quantified against
    anyway.
