##
# Predefined Gaussian channels
##

function displace(::Type{Td}, ::Type{Tt}, alpha::A, noise::N) where {Td,Tt,A<:Number,N}
    disp = sqrt(2) * Td([real(alpha), imag(alpha)])
    transform = Tt(Matrix{Float64}(I, 2, 2))
    return GaussianChannel(disp, transform, Tt(noise), 1)
end
displace(::Type{T}, alpha::A, noise::N) where {T,A<:Number,N} = displace(T, T, alpha, noise)
function displace(alpha::A, noise::N) where {A<:Number,N}
    disp = sqrt(2) * [real(alpha), imag(alpha)]
    transform = Matrix{Float64}(I, 2, 2)
    return GaussianChannel(disp, transform, noise, 1)
end

function squeeze(::Type{Td}, ::Type{Tt}, r::R, theta::R, noise::N) where {Td,Tt,R<:Real,N}
    disp = Td(zeros(2))
    cr, sr = cosh(r), sinh(r)
    t = sinh(r) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    transform = Tt(t)
    return GaussianChannel(disp, transform, Tt(noise), 1)
end
squeeze(::Type{T}, r::R, theta::R, noise::N) where {T,R<:Real,N} = squeeze(T, T, r, theta, noise)
function squeeze(r::R, theta::R, noise::N) where {R<:Real,N}
    disp = zeros(2)
    cr, sr = cosh(r), sinh(r)
    transform = sinh(r) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    return GaussianChannel(disp, transform, noise, 1)
end

function twosqueeze(::Type{Td}, ::Type{Tt}, r::R, theta::R, noise::N) where {Td,Tt,R<:Real,N}
    disp = Td(zeros(4))
    v1 = cosh(r) * Matrix{Float64}(I, 2, 2)
    v2 = sinh(r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    transform = Tt([v1 -v2; -v2 v1])
    return GaussianChannel(disp, transform, Tt(noise), 2)
end
twosqueeze(::Type{T}, r::R, theta::R, noise::N) where {T,R<:Real,N} = twosqueeze(T, T, r, theta, noise)
function twosqueeze(r::R, theta::R, noise::N) where {R<:Real,N}
    disp = zeros(4)
    v1 = cosh(r) * Matrix{Float64}(I, 2, 2)
    v2 = sinh(r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    transform = [v1 -v2; -v2 v1]
    return GaussianChannel(disp, transform, noise, 2)
end

function phaseshift(::Type{Td}, ::Type{Tt}, theta::R, noise::N) where {Td,Tt,R<:Real,N}
    disp = Td(zeros(2))
    transform = Tt([cos(theta) sin(theta); -sin(theta) cos(theta)])
    return GaussianChannel(disp, transform, Tt(noise), 1)
end
phaseshift(::Type{T}, theta::R, noise::N) where {T,R<:Real,N} = phaseshift(T, T, theta, noise)
function phaseshift(theta::R, noise::N) where {R<:Real,N}
    disp = zeros(2)
    transform = [cos(theta) sin(theta); -sin(theta) cos(theta)]
    return GaussianChannel(disp, transform, noise, 1)
end

function beamsplitter(::Type{Td}, ::Type{Tt}, transmit::R, noise::N) where {Td,Tt,R<:Real,N}
    disp = Td(zeros(4))
    I2 = Matrix{Float64}(I, 2, 2)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    transform = Tt([a1*I2 a2*I2; -a2*I2 a1*I2])
    return GaussianChannel(disp, transform, Tt(noise), 2)
end
beamsplitter(::Type{T}, transmit::R, noise::N) where {T,R<:Real,N} = beamsplitter(T, T, transmit, noise)
function beamsplitter(transmit::R, noise::N) where {R<:Real,N}
    disp = zeros(4)
    I2 = Matrix{Float64}(I, 2, 2)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    transform = [a1*I2 a2*I2; -a2*I2 a1*I2]
    return GaussianChannel(disp, transform, noise, 2)
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
julia> attenuator(pi/6, 3)
GaussianChannel for 1 mode.
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
function attenuator(::Type{Td}, ::Type{Tt}, theta::R, n::N) where {Td,Tt,R<:Real,N<:Int}
    disp = zeros(2)
    transform = cos(theta) * Matrix{Float64}(I, 2, 2)
    noise = (sin(theta))^2 * n * Matrix{Float64}(I, 2, 2)
    return GaussianChannel(Td(disp), Tt(transform), Tt(noise), 1)
end
attenuator(::Type{T}, theta::R, n::N) where {T,R<:Real,N<:Int} = attenuator(T, T, theta, n)
function attenuator(theta::R, n::N) where {R<:Real,N<:Int}
    disp = zeros(2)
    transform = cos(theta) * Matrix{Float64}(I, 2, 2)
    noise = (sin(theta))^2 * n * Matrix{Float64}(I, 2, 2)
    return GaussianChannel(disp, transform, noise, 1)
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
julia> amplifier(2.0, 3)
GaussianChannel for 1 mode.
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
function amplifier(::Type{Td}, ::Type{Tt}, r::R, n::N) where {Td,Tt,R<:Real,N<:Int}
    disp = zeros(2)
    transform = cosh(r) * Matrix{Float64}(I, 2, 2)
    noise = (sinh(r))^2 * n * Matrix{Float64}(I, 2, 2)
    return GaussianChannel(Td(disp), Tt(transform), Tt(noise), 1)
end
amplifier(::Type{T}, r::R, n::N) where {T,R<:Real,N<:Int} = amplifier(T, T, theta, n)
function amplifier(r::R, n::N) where {R<:Real,N<:Int}
    disp = zeros(2)
    transform = cosh(r) * Matrix{Float64}(I, 2, 2)
    noise = (sinh(r))^2 * n * Matrix{Float64}(I, 2, 2)
    return GaussianChannel(disp, transform, noise, 1)
end

##
# Predefined operations on Gaussian channels
##

function tensor(::Type{Td}, ::Type{Tt}, op1::GaussianChannel, op2::GaussianChannel) where {Td,Tt}
    disp′, transform′, noise′ = _tensor_fields(op1, op2)
    return GaussianChannel(Td(disp′), Tt(transform′), Tt(noise′), op1.nmodes + op2.nmodes)
end
tensor(::Type{T}, op1::GaussianChannel, op2::GaussianChannel) where {T} = tensor(T, T, op1, op2)
function tensor(op1::GaussianChannel, op2::GaussianChannel)
    disp′, transform′, noise′ = _tensor_fields(op1, op2)
    return GaussianChannel(disp′, transform′, noise′, op1.nmodes + op2.nmodes)
end
function _tensor_fields(op1::GaussianChannel, op2::GaussianChannel)
    disp1, disp2 = op1.disp, op2.disp
    length1, length2 = 2*op1.nmodes, 2*op2.nmodes
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