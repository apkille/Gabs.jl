# [Tutorials](@id Tutorials)

```@meta
DocTestSetup = quote
    using Gabs
end
```

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

The types defined in Gabs.jl, such as [`GaussianState`](@ref), [`GaussianUnitary`](@ref), and
[`GaussianChannel`](@ref) are agnostic array wrappers, therefore custom array types are supported 
in this package. Let us explore this feature in detail, using [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) as a case study.

Say we wanted to create a coherent state that wraps around pure Julia arrays. Then we would call
`coherentstate` with a single complex argument:
```jldoctest
julia> coherentstate(1.0-im)
GaussianState for 1 mode.
mean: 2-element Vector{Float64}:
  1.4142135623730951
 -1.4142135623730951
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
Here, we initialized a [`GaussianState`](@ref) object that wraps around a `Vector` of size two (mean vector)
and a `Matrix` of size two by two (covariance matrix). If we wanted to optimize and boost the performance of our code,
then we would wrap our [`GaussianState`](@ref) around different array types, for instance, an [`SVector`](https://juliaarrays.github.io/StaticArrays.jl/stable/api/#StaticArraysCore.SVector) from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) or a [`SparseVector`](https://docs.julialang.org/en/v1/stdlib/SparseArrays/#SparseArrays.SparseVector) from [SparseArrays.jl](https://github.com/JuliaSparse/SparseArrays.jl). Any defined method in Gabs.jl that (i) creates a custom type or (ii) transforms a custom type
can specify an array type in its first (and second) arguments. Let's see an example with [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl):
```jldoctest
julia> using StaticArrays

julia> state = coherentstate(SVector{2}, SMatrix{2,2}, 1.0-im)
GaussianState for 1 mode.
mean: 2-element SVector{2, Float64} with indices SOneTo(2):
  1.4142135623730951
 -1.4142135623730951
covariance: 2×2 SMatrix{2, 2, Float64, 4} with indices SOneTo(2)×SOneTo(2):
 1.0  0.0
 0.0  1.0

julia> state ⊗ state
GaussianState for 2 modes.
mean: 4-element SVector{4, Float64} with indices SOneTo(4):
  1.4142135623730951
 -1.4142135623730951
  1.4142135623730951
 -1.4142135623730951
covariance: 4×4 SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0
```
In the above example, we created a [`GaussianState`](@ref) containing custom array types from
[StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl), then took the tensor product
of the state with itself, which preserved the custom array types.

!!! note
    If we have an array wrapper that initializes both vectors
    and matrices, then we can specifiy the array type with a single argument. For instance, if we want
    to initialize a state that contains `Array`s holding numbers of type `Float32` rather
    than `Float64`, then we could simply pass `Array{Float32}` to any relevant Gabs.jl method:
    ```jldoctest
    julia> state = displace(Array{Float32}, 1.0-im)
    GaussianUnitary for 1 mode.
    displacement: 2-element Vector{Float32}:
      1.4142135
     -1.4142135
    symplectic: 2×2 Matrix{Float32}:
     1.0  0.0
     0.0  1.0
    ```

## GPU Acceleration

## Multithreading

## Benchmarking and Profiling

