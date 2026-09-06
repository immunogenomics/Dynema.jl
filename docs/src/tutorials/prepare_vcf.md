# Extra: Compressing and indexing a VCF: `dynema-prepare-vcf`

`dynema-map --vcf` needs its VCF to be **BGZF-compressed** (bgzip) and
**tabix-indexed** (`.vcf.gz` + `.vcf.gz.tbi`), so that each gene's
cis-window can be extracted in about a second instead of scanning the whole
file. Imputation pipelines usually hand you exactly that -- but if yours
didn't, `dynema-prepare-vcf` gets you there in one command, using the htslib
binaries bundled with Dynema (**no system `bgzip`/`tabix` install needed**):

```bash
./bin/dynema-prepare-vcf --vcf genotypes.vcf.gz
```

**Your input file is never modified**: any needed compression writes a
*new* file alongside it, named `<name>_dynema.vcf.bgz` (`.bgz` is htslib's
extension for BGZF files). The `_dynema` suffix marks the file as
Dynema-generated, so nothing of yours can be overwritten -- and on a rerun
an existing copy is reused if it still matches the input, or regenerated if
the input has changed since. The three states a VCF arrives in:

- a **plain `.vcf`** is BGZF-compressed to a new `<name>_dynema.vcf.bgz`
  alongside it, then indexed;
- a **`.vcf.gz` that is plain gzip** rather than BGZF -- a common trap: it
  looks identical, decompresses identically, but tabix rejects it -- is
  detected and recompressed to a new `<name>_dynema.vcf.bgz` alongside it,
  then indexed;
- an **already-BGZF `.vcf.gz`** is simply indexed (a no-op if a
  `.tbi`/`.csi` index already exists; pass `--force` to rebuild).

The command prints the path of the ready file at the end -- that is what
you pass to `dynema-map --vcf` (it differs from the input whenever a
compression step ran). Run it once per VCF.

!!! note
    tabix requires a **coordinate-sorted** VCF. Dynema deliberately never
    modifies or reorders a VCF's records -- your primary data stays exactly
    as it is -- so if yours turns out not to be sorted, the command errors
    with instructions to sort it yourself first (e.g.
    `bcftools sort -Oz -o sorted.vcf.gz genotypes.vcf.gz`) and prepare the
    sorted file. Also remember that Dynema reads genotype *dosages* (`DS`)
    or genotype probabilities (`GP`), not hard calls: preparing the file
    only writes a recompressed copy and an index.
