# Gabs.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://apkille.github.io/Gabs.jl/stable)
[![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://apkille.github.io/Gabs.jl/dev)
[![Build status (Github Actions)](https://github.com/apkille/Gabs.jl/workflows/CI/badge.svg)](https://github.com/apkille/Gabs.jl/actions)
[![codecov](https://codecov.io/github/apkille/Gabs.jl/graph/badge.svg?token=JWMOD4FY6P)](https://codecov.io/github/apkille/Gabs.jl)
[![JET static analysis](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Gabs.jl is a numerical tooling package for simulating **Gaussian quantum information**.

Gaussian states and operators have the convenient property that they can be
characterized by low-dimensional matrices in the phase space representation.
Thus, a large class of continuous variable quantum information can be efficiently
simulated on a classical computer, lending to applications in quantum cryptography, quantum machine learning, integrated quantum photonics, and more. Gabs.jl provides a high-level [Julia](https://julialang.org) interface for straightforward and high performance implementations of Gaussian quantum systems.

See the detailed [suggested readings & references page](https://apkille.github.io/Gabs.jl/dev/bibliography/) for a background on quantum information with Gaussian states.

## Installation

To install Gabs.jl, start Julia and run the following command:

```julia
using Pkg
Pkg.add("Gabs")
```
To use the package, run the command

```julia
using Gabs
```

Now, the entire library is loaded into the current workspace, with access to its
high-level interface and predefined objects.

## Examples
### Gaussian States

<details>
 <summaryClick me! ></summary>
<p>

A wide variety of predefined methods to create a specific instance of the [`GaussianState`](https://apkille.github.io/Gabs.jl/dev/API/#Gabs.GaussianState) type are available. For a full description of the API for Gaussian states, see the [State Zoo section](https://apkille.github.io/Gabs.jl/dev/zoos/#State-Zoo) of the documentation. Let's examine a few well-known examples with the Julia REPL:

```julia
julia> using Gabs

julia> sympbasis = QuadPairBasis(1)
QuadPairBasis(1)

julia> s1 = vacuumstate(sympbasis)
GaussianState for 1 mode.
  symplectic basis: QuadPairBasis
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0

julia> s2 = coherentstate(sympbasis, rand(ComplexF64))
GaussianState for 1 mode.
  symplectic basis: QuadPairBasis
mean: 2-element Vector{Float64}:
 1.1000447533324929
 0.38900397266196973
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0

julia> s3 = squeezedstate(sympbasis, rand(Float64), rand(Float64))
GaussianState for 1 mode.
  symplectic basis: QuadPairBasis
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 0.414711   0.0537585
 0.0537585  0.609798
```

Tensor products of Gaussian states are called with `⊗` or `tensor`:

```julia
julia> tp = s1 ⊗ s2 ⊗ s3
GaussianState for 3 modes.
  symplectic basis: QuadPairBasis
mean: 6-element Vector{Float64}:
 0.0
 0.0
 1.1000447533324929
 0.38900397266196973
 0.0
 0.0
covariance: 6×6 Matrix{Float64}:
 1.0  0.0  0.0  0.0  0.0        0.0
 0.0  1.0  0.0  0.0  0.0        0.0
 0.0  0.0  1.0  0.0  0.0        0.0
 0.0  0.0  0.0  1.0  0.0        0.0
 0.0  0.0  0.0  0.0  0.414711   0.0537585
 0.0  0.0  0.0  0.0  0.0537585  0.609798
```

Partial traces of Gaussian states are called with `ptrace`:

```julia
julia> ptrace(tp, [1, 3])
GaussianState for 2 modes.
  symplectic basis: QuadPairBasis
mean: 4-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
covariance: 4×4 Matrix{Float64}:
 1.0  0.0  0.0        0.0
 0.0  1.0  0.0        0.0
 0.0  0.0  0.414711   0.0537585
 0.0  0.0  0.0537585  0.609798
```
</p>
</details>

### Gaussian Operators

<details>
 <summaryClick me! ></summary>
<p>

Gabs.jl contains many predefined methods to create instances of [`GaussianUnitary`](https://apkille.github.io/Gabs.jl/dev/API/#Gabs.GaussianUnitary) and [`GaussianChannel`](https://apkille.github.io/Gabs.jl/dev/API/#Gabs.GaussianChannel) types. For a full description of the API for Gaussian operators, see the [Operator Zoo section](https://apkille.github.io/Gabs.jl/dev/zoos/#Operator-Zoo) of the documentation. Let's see a few well-known examples:

```julia
julia> using Gabs

julia> sympbasis = QuadBlockBasis(1)
QuadBlockBasis(1)

julia> un = displace(sympbasis, rand(ComplexF64))
GaussianUnitary for 1 mode.
  symplectic basis: QuadBlockBasis
displacement: 2-element Vector{Float64}:
 1.0534226538602067
 0.16968892481522801
symplectic: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0

julia> ch = amplifier(sympbasis, rand(Float64), rand(Int64))
GaussianChannel for 1 mode.
  symplectic basis: QuadBlockBasis
displacement: 2-element Vector{Float64}:
 0.0
 0.0
transform: 2×2 Matrix{Float64}:
 1.05403  0.0
 0.0      1.05403
noise: 2×2 Matrix{Float64}:
 5.93927e17  0.0
 0.0         5.93927e17
```

Applications of these operators on states can be called in-place with `apply!` and out-of-place with `*`:

```julia
julia> s = vacuumstate(sympbasis);

julia> apply!(s, un); s
GaussianState for 1 mode.
  symplectic basis: QuadBlockBasis
mean: 2-element Vector{Float64}:
 1.0534226538602067
 0.16968892481522801
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0

julia> tp = (ch ⊗ ch ⊗ ch) * (s ⊗ s ⊗ s)
GaussianState for 3 modes.
  symplectic basis: QuadBlockBasis
mean: 6-element Vector{Float64}:
 1.1103377371945864
 1.1103377371945864
 1.1103377371945864
 0.17885700114374567
 0.17885700114374567
 0.17885700114374567
covariance: 6×6 Matrix{Float64}:
 5.93927e17  0.0         0.0         0.0         0.0         0.0
 0.0         5.93927e17  0.0         0.0         0.0         0.0
 0.0         0.0         5.93927e17  0.0         0.0         0.0
 0.0         0.0         0.0         5.93927e17  0.0         0.0
 0.0         0.0         0.0         0.0         5.93927e17  0.0
 0.0         0.0         0.0         0.0         0.0         5.93927e17
```
</p>
</details>