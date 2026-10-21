```@meta
CurrentModule = Dynema
```

# Internals

These types are not exported, but are documented here for contributors and
for users who want to work with `expand_geno`'s return value directly.

```@docs
Dynema.ExpandedGeno
Dynema.ExpandedGenoView
```

The CRVE score test itself, and the per-context 1-df decomposition reported
by multi-context interaction tests (the `p_<context>` columns; see
`map_locus`'s `percontext` keyword):

```@docs
Dynema.smallsample_factors
Dynema.crvetest_direct
Dynema.crve_percomponent
```

The cluster-robust standard errors of the unrestricted coefficients (the
`se_wald_<term>` columns that accompany the effect-size estimates; see `map_locus`'s
`betas` keyword). These are Wald-type SEs at the unrestricted fit, distinct
from the null-fit score test above:

```@docs
Dynema.crve_se
```

!!! note
    `Dynema.DynemaModel`, the struct returned by [`map_locus`](@ref), is
    intentionally accessed only through its `get_*`/`set_*` accessors (see
    the [API Reference](functions.md)) rather than by touching its fields
    directly, so that the internal layout can change without breaking user
    code.
