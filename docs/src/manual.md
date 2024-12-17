# [Manual](@id Manual)

```@meta
DocTestSetup = quote
    using Gabs
end
```

Simply put, Gabs is a package for creating and transforming Gaussian bosonic systems. This section discusses the "lower level" tools for simulating such phenomena, with
mathematical explanations when appropriate. For comprehensive reviews of Gaussian
quantum information, see the [suggested readings page](@ref References).

## The Symplectic Formalism

The underlying geometry of Gaussian informatics in the phase space is *symplectic*. From the basic canonical commutation
relations (CCRs) of quantized continuous variable systems, manifestations of the symplectic group $\text{Sp}(2N, \mathbb{R})$
appear everywhere. In Gabs, symplectic basis types must be defined from the beginning. Here's how they are laid out in this library:

| canonical ordering | symplectic form | basis type |
| :---: | :---: | :---: |
| $(\hat{x}_1, \hat{p}_1, \cdots, \hat{x}_N, \hat{p}_N)$ | $\begin{pmatrix} 0 & 1 \\ - 1 & 0 \end{pmatrix} \otimes \mathbf{I}_N$ | [`QuadPairBasis`](@ref) |
| $(\hat{x}_1, \cdots, \hat{x}_N, \hat{p}_1, \cdots, \hat{p}_N)$ | $\begin{pmatrix} 0 & \mathbf{I}_N \\ -\mathbf{I}_N & 0 \end{pmatrix}$ | [`QuadBlockBasis`](@ref) |

Each symplectic basis type is wrapped around the number of bosonic modes $N$. We can compose a larger symplectic basis
with [`directsum`](@ref) or `⊕`, the direct sum symbol which can be typed in the Julia REPL as `\oplus<TAB>`:

```jldoctest
julia> b = QuadPairBasis(2)
QuadPairBasis(2)

julia> b ⊕ b
QuadPairBasis(4)
```
Of course, this type of behavior will occur implicitly when we take tensor products of Gaussian states and operators, as discussed
in the following sections.

!!! note
    A matrix $\mathbf{S}$ of size $2N\times 2N$ is symplectic when it satisfies the relation $\mathbf{S} \mathbf{\Omega} \mathbf{S}^{\text{T}} = \mathbf{\Omega},$ where $\mathbf{\Omega}$ is an invertible skew-symmetric matrix known as the *symplectic form*.

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
such state vector simulations, we recommend using the [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) library. For our purposes in the phase space, we can manually create multi-mode Gaussian systems with a tensor product, which can be called with either [`tensor`](@ref) or `⊗`, the Kronecker product symbol
which can be typed in the Julia REPL as `\otimes<TAB>`. Take the following example, where we produce a 3-mode Gaussian state that consists of a coherent state, vacuumstate, and squeezed state:

```jldoctest
julia> basis = QuadPairBasis(1);

julia> coherentstate(basis, -1.0+im) ⊗ vacuumstate(basis) ⊗ squeezedstate(basis, 0.25, pi/4)
GaussianState for 3 modes.
  symplectic basis: QuadPairBasis
mean: 6-element Vector{Float64}:
 -1.4142135623730951
  1.4142135623730951
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

Note that in the above example, we defined the symplectic basis to be of type [`QuadPairBasis`](@ref). If we wanted the canonical field operators to be ordered blockwise, then we would call [`QuadBlockBasis`](@ref) instead:

```jldoctest
julia> basis = QuadBlockBasis(1);

julia> coherentstate(basis, -1.0+im) ⊗ vacuumstate(basis) ⊗ squeezedstate(basis, 0.25, pi/4)
GaussianState for 3 modes.
  symplectic basis: QuadBlockBasis
mean: 6-element Vector{Float64}:
 -1.4142135623730951
  0.0
  0.0
  1.4142135623730951
  0.0
  0.0
covariance: 6×6 Matrix{Float64}:
 1.0  0.0  0.0       0.0  0.0  0.0
 0.0  1.0  0.0       0.0  0.0  0.0
 0.0  0.0  0.379578  0.0  0.0  0.184235
 0.0  0.0  0.0       1.0  0.0  0.0
 0.0  0.0  0.0       0.0  1.0  0.0
 0.0  0.0  0.184235  0.0  0.0  0.748048
```

## Gaussian Operators

To transform Gaussian states into Gaussian states, we need Gaussian maps. There are various ways to construct Gaussian transformations, which we will discuss in this section.

### Gaussian Unitaries

Let's begin with the simplest Gaussian transformation, a unitary transformation, which can be created with the [`GaussianUnitary`](@ref) type:

```@docs; canonical = false
GaussianUnitary
```

This is a rather clean way to characterize a large group of Gaussian transformations on
an `N`-mode Gaussian bosonic system. As long as we have a displacement vector of size `2N` and symplectic matrix of size `2N x 2N`, we can create a Gaussian transformation. 

This library has a number of predefined Gaussian unitaries, which are listed below:

- [`displace`](@ref)
- [`squeeze`](@ref)
- [`twosqueeze`](@ref)
- [`phaseshift`](@ref)
- [`beamsplitter`](@ref)
  
Detailed discussions and mathematical descriptions for each of these unitaries are given in the [Gaussian Zoos](@ref) page.

### Gaussian Channels

Noisy bosonic channels are an important model for describing the interaction between a Gaussian state and its environment. Similar to Gaussian unitaries, Gaussian channels are linear bosonic channels that map Gaussian states to Gaussian states. Such objects can be created with the [`GaussianChannel`](@ref) type:

```@docs; canonical = false
GaussianChannel
```

Listed below are a list of predefined Gaussian channels supported by Gabs:

- [`attenuator`](@ref)
- [`amplifier`](@ref)
  
!!! note
    When its noise matrix $\mathbf{N} = \mathbf{0}$ and transform operator $\mathbf{T}$ is a symplectic matrix, a Gaussian channel is a unitary operator. Any predefined Gaussian unitary
    method can be called with an additional noise matrix to create a [`GaussianChannel`](@ref) object. For instance, a noisy displacement operator can be called with [`displace`](@ref) as follows:

    ```jldoctest
    julia> basis = QuadPairBasis(1);

    julia> noise = [1.0 -2.0; 4.0 -3.0];

    julia> displace(basis, 1.0-im, noise)
    GaussianChannel for 1 mode.
      symplectic basis: QuadPairBasis
    displacement: 2-element Vector{Float64}:
      1.4142135623730951
     -1.4142135623730951
    transform: 2×2 Matrix{Float64}:
     1.0  0.0
     0.0  1.0
    noise: 2×2 Matrix{Float64}:
     1.0  -2.0
     4.0  -3.0
    ```