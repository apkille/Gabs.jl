# [Tutorials](@id Tutorials)

## Visualizations

Quasiprobability distributions of Gaussian states can be visualized with [Makie.jl](https://github.com/MakieOrg/Makie.jl). Gabs.jl
currently has support for the following distributions, which can be called with the keyword argument `dist`:
- [`wigner`](@ref)
- [`wignerchar`](@ref)

Below is a code example that plots the wigner distribution of a vacuum state:
```@example
using Gabs, CairoMakie
state = vacuumstate()
q, p = collect(-5.0:0.5:5.0), collect(-5.0:0.5:5.0) # phase space coordinates
heatmap(q, p, state, dist = :wigner)
```
Of course, more plotting sugar can be added to this example with internal Makie attributes:
```@example
using Gabs, CairoMakie
state = vacuumstate()
q, p = collect(-5.0:0.5:5.0), collect(-5.0:0.5:5.0) # phase space coordinates
fig = Figure(fontsize=20, size = (375, 300), fonts = (; regular="CMU Serif"))
ax = Axis(fig[1,1], xlabel = L"q", ylabel = L"p")
hm = heatmap!(ax, q, p, state, dist = :wigner, colormap = :heat)
Colorbar(fig[1,2], hm)
fig
```
## Using Custom Arrays

## GPU Acceleration

## Multithreading

## Benchmarking and Profiling

