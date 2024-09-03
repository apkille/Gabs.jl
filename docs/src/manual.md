# [Manual](@id Manual)

```@meta
DocTestSetup = quote
    using Gabs
end
```

Simply put, `Gabs.jl` is a package for creating and transforming Gaussian bosonic systems. This section discusses the "lower level" tools for simulating such phenomena.

## Gaussian States

The star of this package is the [`GaussianState`](@ref) type, which allows us to initialize
and manipulate a phase space description of an arbitrary Gaussian state.

```@docs; canonical = false
GaussianState
```

Functions to create instances of elementary Gaussian states are provided as part of the package API. 
Listed below are supported predefined Gaussian states:

- [`vacuumstate`](@ref)
- [`thermalstate`](@ref)
- [`coherentstate`](@ref)
- [`squeezedstate`](@ref)
- [`eprstate`](@ref)

Detailed discussions and mathematical descriptions for each of these states are given in the
[Gaussian Zoos](@ref) page.

If we were operating in the state (Fock) space, and wanted to describe multi-mode Gaussian states,
we would take the tensor product of multiple density operators. That method, however,
is quite computationally expensive and requires a finite truncation of the Fock basis. To create
such state vector simulations, we recommend using the [`QuantumOptics.jl`](https://github.com/qojulia/QuantumOptics.jl) library. For our purposes in the phase space, we can manually create multi-mode Gaussian systems with a direct sum, which can be called with either [`directsum`](@ref) or `⊕`, the direct sum symbol
which can be typed in the Julia REPL as `\oplus<TAB>`. Take the following example, where we
produce a 3-mode Gaussian state that consists of a coherent state, vacuumstate, and squeezed state:

```jldoctest
julia> coherentstate(-1.0) ⊕ vacuumstate() ⊕ squeezedstate(0.25, pi/4)
GaussianState
mean: 6-element Vector{Float64}:
 -1.4142135623730951
  0.0
  0.0
  0.0
  0.0
  0.0
covariance: 6×6 Matrix{Float64}:
 1.0  0.0  0.0  0.0  0.0       0.0
 0.0  1.0  0.0  0.0  0.0       0.0
 0.0  0.0  1.0  0.0  0.0       0.0
 0.0  0.0  0.0  1.0  0.0       0.0
 0.0  0.0  0.0  0.0  0.379578  0.184235
 0.0  0.0  0.0  0.0  0.184235  0.748048
```

## Gaussian Unitary Operators

```@docs; canonical = false
GaussianUnitary
```

## Gaussian Channels

```@docs; canonical = false
GaussianChannel
```