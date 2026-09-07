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
Dynema.crvetest_direct
Dynema.crve_percomponent
```

!!! note
    `Dynema.DynemaModel`, the struct returned by [`map_locus`](@ref), is
    intentionally accessed only through its `get_*`/`set_*` accessors (see
    the [API Reference](functions.md)) rather than by touching its fields
    directly, so that the internal layout can change without breaking user
    code.
