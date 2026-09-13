<img src="docs/src/assets/logo.svg" alt="Dynema.jl logo" width="180" align="right">

# Dynema (Dynamic eQTL mapping for single cells)

[![Build Status](https://github.com/immunogenomics/Dynema.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/immunogenomics/Dynema.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/immunogenomics/Dynema.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://immunogenomics.github.io/Dynema.jl/dev/)
[![Docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://immunogenomics.github.io/Dynema.jl/stable/)

*Dynema* is a method to map single-cell eQTL effects at real cellular resolution.

*Dynema*'s generalized framework enables testing complex regulatory effects including:

- **Main effects**: A standard eQTL effect, independent of any context
- **Interaction effects**: eQTL effects that change depending on one (**single-context**) or multiple (**multi-context**) contexts
- **Total effects**: Joint effect of main and interaction eQTL components. This effect captures any genetic signal driven by either main or interaction eQTL effects


*Dynema* scales to genome-wide analysis, accounts for repeated measurements (multiple cells per donor), and provides calibrated *p*-values by using cluster-robust variance estimators (CRVEs). Additionally, it provides robust inferences in extreme scenarios such as small number of donors via optional adaptive score bootstrapping (built on [WildBootTests.jl](https://github.com/droodman/WildBootTests.jl)).


# Installation

*Dynema* requires Julia >= 1.9 (easily installed with [juliaup](https://github.com/julialang/juliaup)).

**You don't need to write any Julia code to run *Dynema*!** The `bin/` folder ships a self-contained command-line interface (`dynema-map`) that installs *Dynema* and all dependencies for you.

```bash
curl -L https://github.com/immunogenomics/Dynema.jl/archive/refs/heads/main.tar.gz | tar -xz
mv Dynema.jl-main Dynema.jl
cd Dynema.jl
./bin/dynema-map --help
```

For more details and tutorials, see the [documentation website](https://immunogenomics.github.io/Dynema.jl/dev/).
