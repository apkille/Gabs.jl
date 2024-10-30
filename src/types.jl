"""
Defines a Gaussian state for an N-mode bosonic system over a 2N-dimensional phase space.

## Fields

- `mean`: The mean vector of length 2N.
- `covar`: The covariance matrix of size 2N x 2N.
- `nmodes`: The number of modes N of the Gaussian state.

## Mathematical description of a Gaussian state

An ``N``-mode Gaussian state, ``\\hat{\\rho}(\\mathbf{\\bar{x}}, \\mathbf{V})``, is a density
operator characterized by two statistical moments: a mean vector ``\\mathbf{\\bar{x}}`` of
length ``2N`` and covariance matrix ``\\mathbf{V}`` of size ``2N\\times 2N``. By definition,
the Wigner representation of a Gaussian state is a Gaussian function.

## Example

```jldoctest
julia> coherentstate(1.0+im)
GaussianState for 1 mode.
mean: 2-element Vector{Float64}:
 1.4142135623730951
 1.4142135623730951
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
"""
struct GaussianState{M,V,N<:Int} <: StateVector{M,V}
    mean::M
    covar::V
    nmodes::N
    function GaussianState(m::M, v::V, n::N) where {M,V,N}
        all(size(v) .== length(m) .== 2*n) || throw(DimensionMismatch(STATE_ERROR))
        return new{M,V,N}(m, v, n)
    end
end
GaussianState(mean::M, covar::V) where {M,V} = GaussianState(mean, covar, Int(length(mean)/2))

Base.:(==)(x::GaussianState, y::GaussianState) = x.mean == y.mean && x.covar == y.covar
Base.isapprox(x::GaussianState, y::GaussianState) = isapprox(x.mean,y.mean) && isapprox(x.covar,y.covar)
function Base.show(io::IO, mime::MIME"text/plain", x::GaussianState)
    Base.summary(io, x)
    print(io, "\nmean: ")
    Base.show(io, mime, x.mean)
    print(io, "\ncovariance: ")
    Base.show(io, mime, x.covar)
end


"""
Defines a Gaussian unitary for an N-mode bosonic system over a 2N-dimensional phase space.

## Fields

- `disp`: The displacement vector of length 2N.
- `symplectic`: The symplectic matrix of size 2N x 2N.
- `nmodes`: The number of modes N for the Gaussian unitary.

## Mathematical description of a Gaussian unitary

An ``N``-mode Gaussian unitary, ``U(\\mathbf{d}, \\mathbf{S})``, is a unitary
operator characterized by a displacement vector ``\\mathbf{d}`` of length ``2N`` and symplectic
matrix ``\\mathbf{S}`` of size ``2N\\times 2N``, such that its action on a Gaussian state
results in a Gaussian state. More specifically, a Gaussian unitary transformation on a
Gaussian state ``\\hat{\\rho}(\\mathbf{\\bar{x}}, \\mathbf{V})`` is described by its maps on
the statistical moments of the Gaussian state:

```math
\\mathbf{\\bar{x}} \\to \\mathbf{S} \\mathbf{\\bar{x}} + \\mathbf{d}, \\quad
\\mathbf{V} \\to \\mathbf{S} \\mathbf{V} \\mathbf{S}^{\\text{T}}.
```

## Example

```jldoctest
julia> displace(1.0+im)
GaussianUnitary for 1 mode.
displacement: 2-element Vector{Float64}:
 1.4142135623730951
 1.4142135623730951
symplectic: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
"""
struct GaussianUnitary{D,S,N<:Int} <: AbstractOperator{D,S}
    disp::D
    symplectic::S
    nmodes::N
    function GaussianUnitary(d::D, s::S, n::N) where {D,S,N}
        all(size(s) .== length(d) .== 2*n) || throw(DimensionMismatch(UNITARY_ERROR))
        return new{D,S,N}(d, s, n)
    end
end
GaussianUnitary(disp::D, symplectic::S) where {D,S} = GaussianUnitary(disp, symplectic, Int(length(disp)/2))

Base.:(==)(x::GaussianUnitary, y::GaussianUnitary) = x.disp == y.disp && x.symplectic == y.symplectic
function Base.show(io::IO, mime::MIME"text/plain", x::GaussianUnitary)
    Base.summary(io, x)
    print(io, "\ndisplacement: ")
    Base.show(io, mime, x.disp)
    print(io, "\nsymplectic: ")
    Base.show(io, mime, x.symplectic)
end

function Base.:(*)(op::GaussianUnitary, state::GaussianState)
    d, S, = op.disp, op.symplectic
    op.nmodes == state.nmodes || throw(DimensionMismatch(ACTION_ERROR))
    mean′ = S * state.mean .+ d
    covar′ = S * state.covar * transpose(S)
    return GaussianState(mean′, covar′, state.nmodes)
end
function apply!(state::GaussianState, op::GaussianUnitary)
    d, S = op.disp, op.symplectic
    state.nmodes == op.nmodes || throw(DimensionMismatch(ACTION_ERROR))
    state.mean .= S * state.mean .+ d
    state.covar .= S * state.covar * transpose(S)
    return state
end

"""
Defines a Gaussian channel for an N-mode bosonic system over a 2N-dimensional phase space.

## Fields

- `disp`: The displacement vector of length 2N.
- `transform`: The transformation matrix of size 2N x 2N.
- `noise`: The noise matrix of size 2N x 2N.
- `nmodes`: The number of modes N for the Gaussian channel.

## Mathematical description of a Gaussian channel

An ``N``-mode Gaussian channel, ``G(\\mathbf{d}, \\mathbf{T}, \\mathbf{N})``, is an
operator characterized by a displacement vector ``\\mathbf{d}`` of length ``2N``, as well as
a transformation matrix ``\\mathbf{T}`` and noise matrix ``\\mathbf{N}`` of size ``2N\\times 2N``,
such that its action on a Gaussian state results in a Gaussian state. More specifically, a Gaussian
channel action on a Gaussian state ``\\hat{\\rho}(\\mathbf{\\bar{x}}, \\mathbf{V})`` is
described by its maps on the statistical moments of the Gaussian state:

```math
\\mathbf{\\bar{x}} \\to \\mathbf{S} \\mathbf{\\bar{x}} + \\mathbf{d}, \\quad
\\mathbf{V} \\to \\mathbf{T} \\mathbf{V} \\mathbf{T}^{\\text{T}} + \\mathbf{N}.
```

## Example

```jldoctest
julia> noise = [1.0 -3.0; 4.0 2.0];

julia> displace(1.0+im, noise)
GaussianChannel for 1 mode.
displacement: 2-element Vector{Float64}:
 1.4142135623730951
 1.4142135623730951
transform: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
noise: 2×2 Matrix{Float64}:
 1.0  -3.0
 4.0   2.0
```
"""
struct GaussianChannel{D,T,N<:Int} <: AbstractOperator{D,T}
    disp::D
    transform::T
    noise::T
    nmodes::N
    function GaussianChannel(d::D, t::T, n::T, m::N) where {D,T,N}
        all(length(d) .== size(t) .== size(n) .== 2*m) || throw(DimensionMismatch(CHANNEL_ERROR))
        return new{D,T,N}(d, t, n, m)
    end
end
GaussianChannel(disp::D, transform::T, noise::T) where {D,T} = GaussianChannel(disp, transform, noise, Int(length(disp)/2))

Base.:(==)(x::GaussianChannel, y::GaussianChannel) = x.disp == y.disp && x.transform == y.transform && x.noise == y.noise
function Base.show(io::IO, mime::MIME"text/plain", x::GaussianChannel)
    Base.summary(io, x)
    print(io, "\ndisplacement: ")
    Base.show(io, mime, x.disp)
    print(io, "\ntransform: ")
    Base.show(io, mime, x.transform)
    print(io, "\nnoise: ")
    Base.show(io, mime, x.noise)
end
function Base.summary(io::IO, x::Union{GaussianState,GaussianUnitary,GaussianChannel})
    printstyled(io, nameof(typeof(x)); color=:blue)
    modenum = x.nmodes
    if isone(modenum)
        print(io, " for $(modenum) mode.")
    else
        print(io, " for $(modenum) modes.")
    end
end

function Base.:(*)(op::GaussianChannel, state::GaussianState)
    d, T, N = op.disp, op.transform, op.noise
    op.nmodes == state.nmodes || throw(DimensionMismatch(ACTION_ERROR))
    mean′ = T * state.mean .+ d
    covar′ = T * state.covar * transpose(T) .+ N
    return GaussianState(mean′, covar′, state.nmodes)
end
function apply!(state::GaussianState, op::GaussianChannel)
    d, T, N = op.disp, op.transform, op.noise
    state.nmodes == op.nmodes || throw(DimensionMismatch(ACTION_ERROR))
    state.mean .= T * state.mean .+ d
    state.covar .= T * state.covar * transpose(T) .+ N
    return state
end