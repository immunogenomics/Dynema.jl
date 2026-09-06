```@meta
CurrentModule = Dynema
```

# Dynema.jl

See our pre-print using this [link](https://www.biorxiv.org/content/10.64898/2026.08.25.747138v1)

*Dynema* (Dynamic eQTL mapping in single cells) is a fast and calibrated method to map eQTL effects genome-wide at true single-cell resolution. We argue that an eQTL effect is different for each individual cell. Dynema can decompose such effects into: i) context-independent (main), ii) context-dependent (interaction), and iii) total (main and interaction).

**No Julia coding required:** Dynema ships a command-line interface that runs the full mapping workflow -- straight from a VCF and a Matrix Market count matrix to per-gene summary statistics tables, in batches of genes ready to parallelize on any HPC scheduler -- entirely from the terminal. See the [Command-line overview](tutorials/command_line.md) to get started. A Julia API is also available for anyone who wants to call Dynema directly from their own scripts or pipelines -- see the [API Reference](functions.md).

## Getting started

- [Quick start (end to end)](tutorials/quick_start.md) from input files
  (VCF + Matrix Market expression) through data indexing, ad single-cell  eQTL tests on a small simulated dataset.
- [Command-line overview](tutorials/command_line.md) familiarizes with 
  running *Dynema* and its parameters from a VCF and a Matrix Market count 
  matrix, without writing any Julia code.
- [API Reference](functions.md) documents every exported function, for
  anyone calling Dynema directly from Julia.


## Getting help

Please open an issue on [GitHub](https://github.com/immunogenomics/Dynema.jl/issues)
for bug reports or feature requests.
