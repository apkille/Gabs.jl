##
# Predefined Gaussian channels
##

function displace(::Type{Td}, ::Type{Tt}, repr::SymplecticRepr{N}, alpha::A, noise::M) where {Td,Tt,N<:Int,A,M}
    disp, transform = _displace(repr, alpha)
    return GaussianChannel(repr, Td(disp), Tt(transform), Tt(noise))
end
displace(::Type{T}, repr::SymplecticRepr{N}, alpha::A, noise::M) where {T,N<:Int,A,M} = displace(T, T, repr, alpha, noise)
function displace(repr::SymplecticRepr{N}, alpha::A, noise::M) where {N<:Int,A,M}
    disp, transform = _displace(repr, alpha)
    return GaussianChannel(repr, disp, transform, noise)
end

function squeeze(::Type{Td}, ::Type{Tt}, repr::SymplecticRepr{N}, r::R, theta::R, noise::M) where {Td,Tt,N<:Int,R,M}
    disp, transform = _squeeze(repr, r, theta)
    return GaussianChannel(repr, Td(disp), Tt(transform), Tt(noise))
end
squeeze(::Type{T}, repr::SymplecticRepr{N}, r::R, theta::R, noise::M) where {T,N<:Int,R,M} = squeeze(T, T, repr, r, theta, noise)
function squeeze(repr::SymplecticRepr{N}, r::R, theta::R, noise::M) where {N<:Int,R,M}
    disp, transform = _squeeze(repr, r, theta)
    return GaussianChannel(repr, disp, transform, noise)
end

function twosqueeze(::Type{Td}, ::Type{Tt}, repr::SymplecticRepr{N}, r::R, theta::R, noise::M) where {Td,Tt,N<:Int,R,M}
    disp, transform = _twosqueeze(repr, r, theta)
    return GaussianChannel(repr, Td(disp), Tt(transform), Tt(noise))
end
twosqueeze(::Type{T}, repr::SymplecticRepr{N}, r::R, theta::R, noise::M) where {T,N<:Int,R,M} = twosqueeze(T, T, repr, r, theta, noise)
function twosqueeze(repr::SymplecticRepr{N}, r::R, theta::R, noise::M) where {N<:Int,R,M}
    disp, transform = _twosqueeze(repr, r, theta)
    return GaussianChannel(repr, disp, transform, noise)
end

function phaseshift(::Type{Td}, ::Type{Tt}, repr::SymplecticRepr{N}, theta::R, noise::M) where {Td,Tt,N<:Int,R,M}
    disp, transform = _phaseshift(repr, theta)
    return GaussianChannel(repr, Td(disp), Tt(transform), Tt(noise))
end
phaseshift(::Type{T}, repr::SymplecticRepr{N}, theta::R, noise::M) where {T,N<:Int,R,M} = phaseshift(T, T, repr, theta, noise)
function phaseshift(repr::SymplecticRepr{N}, theta::R, noise::M) where {N<:Int,R,M}
    disp, transform = _phaseshift(repr, theta)
    return GaussianChannel(repr, disp, transform, noise)
end

function beamsplitter(::Type{Td}, ::Type{Tt}, repr::SymplecticRepr{N}, transmit::R, noise::M) where {Td,Tt,N<:Int,R,M}
    disp, transform = _beamsplitter(repr, transmit)
    return GaussianChannel(repr, Td(disp), Tt(transform), Tt(noise))
end
beamsplitter(::Type{T}, repr::SymplecticRepr{N}, transmit::R, noise::M) where {T,N<:Int,R,M} = beamsplitter(T, T, repr, transmit, noise)
function beamsplitter(repr::SymplecticRepr{N}, transmit::R, noise::M) where {N<:Int,R,M}
    disp, transform = _beamsplitter(repr, transmit)
    return GaussianChannel(repr, disp, transform, noise)
end

"""
    attenuator([Td=Vector{Float64}, Tt=Matrix{Float64},] theta<:Real, n<:Int)

Gaussian channel describing the coupling of an input
single mode Gaussian state and its environment via a beam splitter operation. The channel is paramatrized
by beam splitter rotation angle `theta` and thermal noise `n`.

## Mathematical description of an attenuator channel

An attenuator channel, ``\\mathcal{E}_{\\theta}^{n_{\\text{th}}}``, where ``\\theta`` is
the beam splitter rotation parameter and ``n_{\\text{th}} \\geq 1`` is the thermal noise parameter,
is characterized by the displacement vector ``\\mathbf{d}``, transformation matrix ``\\mathbf{T}``,
and noise matrix ``\\mathbf{N}``, expressed respectively as follows:

```math
\\mathbf{d} = \\mathbf{0},
\\quad \\mathbf{T} = \\cos\\theta\\mathbf{I},
\\qquad \\mathbf{N} = (\\sin\\theta)^2 n_{\\text{th}} \\mathbf{I}.
```

## Example

```jldoctest
julia> attenuator(CanonicalForm(1), pi/6, 3)
GaussianChannel for 1 mode in CanonicalForm representation.
displacement: 2-element Vector{Float64}:
 0.0
 0.0
transform: 2×2 Matrix{Float64}:
 0.866025  0.0
 0.0       0.866025
noise: 2×2 Matrix{Float64}:
 0.75  0.0
 0.0   0.75
```
"""
function attenuator(::Type{Td}, ::Type{Tt}, repr::SymplecticRepr{N}, theta::R, n::M) where {Td,Tt,N<:Int,R,M}
    disp, transform, noise = _attenuator(repr, theta, n)
    return GaussianChannel(repr, Td(disp), Tt(transform), Tt(noise))
end
attenuator(::Type{T}, repr::SymplecticRepr{N}, theta::R, n::M) where {T,N<:Int,R,M} = attenuator(T, T, repr, theta, n)
function attenuator(repr::SymplecticRepr{N}, theta::R, n::M) where {N<:Int,R,M}
    disp, transform, noise = _attenuator(repr, theta, n)
    return GaussianChannel(repr, disp, transform, noise)
end
function _attenuator(repr::CanonicalForm{N}, theta::R, n::M) where {N<:Int,R<:Real,M<:Int}
    nmodes = repr.nmodes
    disp = zeros(2*nmodes) 
    transform = Matrix{Float64}(cos(theta) * I, 2*nmodes, 2*nmodes)
    noise = Matrix{Float64}((sin(theta))^2 * n * I, 2*nmodes, 2*nmodes)
    return disp, transform, noise
end
function _attenuator(repr::CanonicalForm{N}, theta::R, n::M) where {N<:Int,R<:Vector,M<:Vector}
    nmodes = repr.nmodes
    disp = zeros(2*nmodes) 
    transform = zeros(2*nmodes, 2*nmodes)
    noise = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        ct, st = cos(theta[i]), sin(theta[i])
        ni = n[i]

        transform[2*i-1, 2*i-1] = ct
        transform[2*i, 2*i] = ct

        noise[2*i-1, 2*i-1] = st^2 * ni
        noise[2*i, 2*i] = st^2 * ni
    end
    return disp, transform, noise
end

"""
    amplifier([Td=Vector{Float64}, Tt=Matrix{Float64},] r<:Real, n<:Int)

Gaussian channel describing the interaction of an input
single mode Gaussian state and its environment via a two-mode squeezing operation. The channel is paramatrized
by squeezing amplitude parameter `r` and thermal noise `n`.

## Mathematical description of an amplifier channel

An amplifier channel, ``\\mathcal{A}_{\\theta}^{n_{\\text{th}}}``, where ``r`` is
the squeezing amplitude parameter and ``n_{\\text{th}} \\geq 1`` is the thermal noise parameter,
is characterized by the displacement vector ``\\mathbf{d}``, transformation matrix ``\\mathbf{T}``,
and noise matrix ``\\mathbf{N}``, expressed respectively as follows:

```math
\\mathbf{d} = \\mathbf{0},
\\quad \\mathbf{T} = \\cosh r\\mathbf{I},
\\qquad \\mathbf{N} = (\\sinh r)^2 n_{\\text{th}} \\mathbf{I}.
```

## Example

```jldoctest
julia> amplifier(CanonicalForm(1), 2.0, 3)
GaussianChannel for 1 mode in CanonicalForm representation.
displacement: 2-element Vector{Float64}:
 0.0
 0.0
transform: 2×2 Matrix{Float64}:
 3.7622  0.0
 0.0     3.7622
noise: 2×2 Matrix{Float64}:
 39.4623   0.0
  0.0     39.4623
```
"""
function amplifier(::Type{Td}, ::Type{Tt}, repr::SymplecticRepr{N}, r::R, n::M) where {Td,Tt,N<:Int,R,M}
    disp, transform, noise = _amplifier(repr, r, n)
    return GaussianChannel(repr, Td(disp), Tt(transform), Tt(noise))
end
amplifier(::Type{T}, repr::SymplecticRepr{N}, r::R, n::M) where {T,N<:Int,R,M} = amplifier(T, T, repr, r, n)
function amplifier(repr::SymplecticRepr{N}, r::R, n::M) where {N<:Int,R,M}
    disp, transform, noise = _amplifier(repr, r, n)
    return GaussianChannel(repr, disp, transform, noise)
end
function _amplifier(repr::CanonicalForm{N}, r::R, n::M) where {N<:Int,R<:Real,M<:Int}
    nmodes = repr.nmodes
    disp = zeros(2*nmodes) 
    transform = Matrix{Float64}(cosh(r) * I, 2*nmodes, 2*nmodes)
    noise = Matrix{Float64}((sinh(r))^2 * n * I, 2*nmodes, 2*nmodes)
    return disp, transform, noise
end
function _amplifier(repr::CanonicalForm{N}, r::R, n::M) where {N<:Int,R<:Vector,M<:Vector}
    nmodes = repr.nmodes
    disp = zeros(2*nmodes) 
    transform = zeros(2*nmodes, 2*nmodes)
    noise = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        cr, sr = cos(r[i]), sin(r[i])
        ni = n[i]

        transform[2*i-1, 2*i-1] = cr
        transform[2*i, 2*i] = cr

        noise[2*i-1, 2*i-1] = sr^2 * ni
        noise[2*i, 2*i] = sr^2 * ni
    end
    return disp, transform, noise
end

##
# Predefined operations on Gaussian channels
##

function tensor(::Type{Td}, ::Type{Tt}, op1::GaussianChannel, op2::GaussianChannel) where {Td,Tt}
    disp, transform, noise = _tensor(op1, op2)
    return GaussianChannel(op1.repr + op2.repr, Td(disp), Tt(transform), Tt(noise))
end
tensor(::Type{T}, op1::GaussianChannel, op2::GaussianChannel) where {T} = tensor(T, T, op1, op2)
function tensor(op1::GaussianChannel, op2::GaussianChannel)
    disp, transform, noise = _tensor(op1, op2)
    return GaussianChannel(op1.repr + op2.repr, disp, transform, noise)
end
function _tensor(op1::GaussianChannel, op2::GaussianChannel)
    disp1, disp2 = op1.disp, op2.disp
    repr1, repr2 = op1.repr, op2.repr
    length1, length2 = 2*repr1.nmodes, 2*repr2.nmodes
    slengths = length1 + length2
    trans1, trans2 = op1.transform, op2.transform
    # initialize direct sum of displacement vectors
    disp′ = zeros(slengths)
    @inbounds for i in eachindex(disp1)
        disp′[i] = disp1[i]
    end
    @inbounds for i in eachindex(disp2)
        disp′[i+length1] = disp2[i]
    end
    # initialize direct sum of transform matrix
    transform′ = zeros(slengths, slengths)
    axes1 = (Base.OneTo(length1), Base.OneTo(length1))
    @inbounds for i in axes1[1], j in axes1[2]
        transform′[i,j] = trans1[i,j]
    end
    axes2 = (Base.OneTo(length2), Base.OneTo(length2))
    @inbounds for i in axes2[1], j in axes2[2]
        transform′[i+length1,j+length1] = trans2[i,j]
    end
    noise1, noise2 = op1.noise, op2.noise
    # initialize direct sum of noise matrix
    noise′ = zeros(slengths, slengths)
    @inbounds for i in axes1[1], j in axes1[2]
        noise′[i,j] = noise1[i,j]
    end
    @inbounds for i in axes2[1], j in axes2[2]
        noise′[i+length1,j+length1] = noise2[i,j]
    end
    # extract output array types
    disp′′ = _promote_output_vector(typeof(disp1), typeof(disp2), disp′)
    transform′′ = _promote_output_matrix(typeof(trans1), typeof(trans2), transform′)
    noise′′ = _promote_output_matrix(typeof(noise1), typeof(noise2), noise′)
    return disp′′, transform′′, noise′′
end