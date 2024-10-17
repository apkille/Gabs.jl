# [Tutorials](@id Tutorials)

## Visualizations

Quasiprobability distributions of Gaussian states can be visualized with [Makie.jl](https://github.com/MakieOrg/Makie.jl). Gabs.jl
currently has support for the following distributions, which can be called with the keyword argument `dist`:
- [`wigner`](@ref)
- [`wignerchar`](@ref)

```@example
using Gabs, CairoMakie
state = vacuumstate()
q, p = collect(-5.0:0.5:5.0), collect(-5.0:0.5:5.0) # phase space coordinates
heatmap(q, p, state, dist = :wigner)
```

## Using Custom Arrays

## GPU Acceleration

## Multithreading

## Benchmarking and Profiling

