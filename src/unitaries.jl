##
# Predefined Gaussian unitaries
##

"""
    displace([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, alpha<:Number)
    displace([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, alpha<:Number, noise::Ts)

Gaussian operator that displaces the vacuum state into a coherent state, known
as the displacement operator. The symplectic representation is given by `basis`. The complex amplitude is given by `alpha`. Noise can
be added to the operation with `noise`.

## Mathematical description of a displacement operator

A displacement operator ``D(\\alpha)`` is defined by the operation
``D(\\alpha)|0\\rangle = |\\alpha\\rangle``, where ``\\alpha`` is
a complex amplitude. The operator ``D(\\alpha)`` is characterized by 
the displacement vector ``\\mathbf{d}`` and symplectic
matrix ``\\mathbf{S}``, expressed respectively as follows:

```math
\\mathbf{d} = \\sqrt{2}\\left(\\text{Re}(\\alpha), \\text{Im}(\\alpha)\\right)^{\\text{T}},
\\quad \\mathbf{S} = \\mathbf{I}.
```

## Example

```jldoctest
julia> displace(QuadPairBasis(1), 1.0+im)
GaussianUnitary for 1 mode in QuadPairBasis representation.
displacement: 2-element Vector{Float64}:
 1.4142135623730951
 1.4142135623730951
symplectic: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
"""
function displace(::Type{Td}, ::Type{Ts}, basis::SymplecticBasis{N}, alpha::A) where {Td,Ts,N<:Int,A}
    disp, symplectic = _displace(basis, alpha)
    return GaussianUnitary(basis, Td(disp), Ts(symplectic))
end
displace(::Type{T}, basis::SymplecticBasis{N}, alpha::A) where {T,N<:Int,A} = displace(T, T, basis, alpha)
function displace(basis::SymplecticBasis{N}, alpha::A) where {N<:Int,A}
    disp, symplectic = _displace(basis, alpha)
    return GaussianUnitary(basis, disp, symplectic)
end
function _displace(basis::QuadPairBasis{N}, alpha::A) where {N<:Int,A<:Number}
    nmodes = basis.nmodes
    disp = repeat([sqrt(2)*real(alpha), sqrt(2)*imag(alpha)], nmodes)
    symplectic = Matrix{Float64}(I, 2*nmodes, 2*nmodes)
    return disp, symplectic
end
function _displace(basis::QuadPairBasis{N}, alpha::A) where {N<:Int,A<:Vector}
    nmodes = basis.nmodes
    disp = sqrt(2) * reinterpret(Float64, alpha)
    symplectic = Matrix{Float64}(I, 2*nmodes, 2*nmodes)
    return disp, symplectic
end
function _displace(basis::QuadBlockBasis{N}, alpha::A) where {N<:Int,A<:Number}
    nmodes = basis.nmodes
    disp = repeat([sqrt(2)*real(alpha), sqrt(2)*imag(alpha)], inner = nmodes)
    symplectic = Matrix{Float64}(I, 2*nmodes, 2*nmodes)
    return disp, symplectic
end
function _displace(basis::QuadBlockBasis{N}, alpha::A) where {N<:Int,A<:Vector}
    nmodes = basis.nmodes
    re = reinterpret(Float64, alpha)
    disp = vcat(@view(re[1:2:end]), @view(re[2:2:end]))
    disp .*= sqrt(2)
    symplectic = Matrix{Float64}(I, 2*nmodes, 2*nmodes)
    return disp, symplectic
end

"""
    squeeze([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, r<:Real, theta<:Real)
    squeeze([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, r<:Real, theta<:Real, noise::Ts)

Gaussian operator that squeezes the vacuum state into a squeezed state, known
as the squeezing operator. The symplectic representation is given by `basis`. The amplitude and phase squeezing parameters 
are given by `r` and `theta`, respectively. Noise can be added to the operation
with `noise`.

## Mathematical description of a squeezing operator

A squeeze operator ``S(r, \\theta)`` is defined by the operation
``S(r, \\theta)|0\\rangle = |r, \\theta\\rangle``, where ``r`` and ``\\theta``
are the real amplitude and phase parameters, respectively. The operator 
``S(r, \\theta)`` is characterized by 
the displacement vector ``\\mathbf{d}`` and symplectic
matrix ``\\mathbf{S}``, expressed respectively as follows:

```math
\\mathbf{d} = \\mathbf{0},
\\quad \\mathbf{S} = \\cosh(r)\\mathbf{I} - \\sinh(r)\\mathbf{R}(\\theta),
```

where ``\\mathbf{R}(\\theta)`` is the rotation matrix.

## Example

```jldoctest
julia> squeeze(QuadPairBasis(1), 0.25, pi/4)
GaussianUnitary for 1 mode in QuadPairBasis representation.
displacement: 2-element Vector{Float64}:
 0.0
 0.0
symplectic: 2×2 Matrix{Float64}:
  0.852789  -0.178624
 -0.178624   1.21004
```
"""
function squeeze(::Type{Td}, ::Type{Ts}, basis::SymplecticBasis{N}, r::R, theta::R) where {Td,Ts,N<:Int,R}
    disp, symplectic = _squeeze(basis, r, theta)
    return GaussianUnitary(basis, Td(disp), Ts(symplectic))
end
squeeze(::Type{T}, basis::SymplecticBasis{N}, r::R, theta::R) where {T,N<:Int,R} = squeeze(T, T, basis, r, theta)
function squeeze(basis::SymplecticBasis{N}, r::R, theta::R) where {N<:Int, R}
    disp, symplectic = _squeeze(basis, r, theta)
    return GaussianUnitary(basis, disp, symplectic)
end
function _squeeze(basis::QuadPairBasis{N}, r::R, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    cr, sr = cosh(r), sinh(r)
    ct, st = cos(theta), sin(theta)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        symplectic[2*i-1, 2*i-1] = cr - sr*ct
        symplectic[2*i-1, 2*i] = -sr * st
        symplectic[2*i, 2*i-1] = -sr * st
        symplectic[2*i, 2*i] = cr + sr*ct
    end
    return disp, symplectic
end
function _squeeze(basis::QuadPairBasis{N}, r::R, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        cr, sr = cosh(r[i]), sinh(r[i])
        ct, st = cos(theta[i]), sin(theta[i])
        symplectic[2*i-1, 2*i-1] = cr - sr*ct
        symplectic[2*i-1, 2*i] = -sr * st
        symplectic[2*i, 2*i-1] = -sr * st
        symplectic[2*i, 2*i] = cr + sr*ct
    end
    return disp, symplectic
end
function _squeeze(basis::QuadBlockBasis{N}, r::R, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    cr, sr = cosh(r), sinh(r)
    ct, st = cos(theta), sin(theta)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        symplectic[i, i] = cr - sr*ct
        symplectic[i, i+nmodes] = -sr * st
        symplectic[i+nmodes, i] = -sr * st
        symplectic[i+nmodes, i+nmodes] = cr + sr*ct
    end
    return disp, symplectic
end
function _squeeze(basis::QuadBlockBasis{N}, r::R, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        cr, sr = cosh(r[i]), sinh(r[i])
        ct, st = cos(theta[i]), sin(theta[i])
        symplectic[i, i] = cr - sr*ct
        symplectic[i, i+nmodes] = -sr * st
        symplectic[i+nmodes, i] = -sr * st
        symplectic[i+nmodes, i+nmodes] = cr + sr*ct
    end
    return disp, symplectic
end

"""
    twosqueeze([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, r<:Real, theta<:Real)
    twosqueeze([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, r<:Real, theta<:Real, noise::Ts)

Gaussian operator that squeezes a two-mode vacuum state into a two-mode squeezed state,
known as the two-mode squeezing operator. The symplectic representation is given by `basis`. The amplitude and phase squeezing parameters 
are given by `r` and `theta`, respectively. Noise can be added to the operation
with `noise`.

## Mathematical description of a two-mode squeezing operator

A two-mode squeeze operator ``S_2(r, \\theta)`` is defined by the operation
``S_2(r, \\theta)|0\\rangle = |r, \\theta\\rangle``, where ``r`` and ``\\theta``
are the real amplitude and phase parameters, respectively. The operator 
``S_2(r, \\theta)`` is characterized by 
the displacement vector ``\\mathbf{d}`` and symplectic
matrix ``\\mathbf{S}``, expressed respectively as follows:

```math
\\mathbf{d} = \\mathbf{0},
\\quad \\mathbf{S} = \\begin{pmatrix}
                    \\cosh(r)\\mathbf{I} & -\\sinh(r)\\mathbf{R}(\\theta) \\\\
                    -\\sinh(r)\\mathbf{R}(\\theta) & \\cosh(r)\\mathbf{I} \\\\
                    \\end{pmatrix},
```

where ``\\mathbf{R}(\\theta)`` is the rotation matrix.

## Example

```jldoctest
julia> twosqueeze(QuadPairBasis(2), 0.25, pi/4)
GaussianUnitary for 2 modes in QuadPairBasis representation.
displacement: 4-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
symplectic: 4×4 Matrix{Float64}:
  1.03141    0.0       -0.178624  -0.178624
  0.0        1.03141   -0.178624   0.178624
 -0.178624  -0.178624   1.03141    0.0
 -0.178624   0.178624   0.0        1.03141
```
"""
function twosqueeze(::Type{Td}, ::Type{Ts}, basis::SymplecticBasis{N}, r::R, theta::R) where {Td,Ts,N<:Int,R}
    disp, symplectic = _twosqueeze(basis, r, theta)
    return GaussianUnitary(basis, Td(disp), Ts(symplectic))
end
twosqueeze(::Type{T}, basis::SymplecticBasis{N}, r::R, theta::R) where {T,N<:Int,R} = twosqueeze(T, T, basis, r, theta)
function twosqueeze(basis::SymplecticBasis{N}, r::R, theta::R) where {N<:Int,R}
    disp, symplectic = _twosqueeze(basis, r, theta)
    return GaussianUnitary(basis, disp, symplectic)
end
function _twosqueeze(basis::QuadPairBasis{N}, r::R, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    cr, sr = cosh(r), sinh(r)
    ct, st = cos(theta), sin(theta)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        symplectic[4*i-3, 4*i-3] = cr
        symplectic[4*i-3, 4*i-1] = -sr * ct
        symplectic[4*i-3, 4*i] = -sr * st

        symplectic[4*i-2, 4*i-2] = cr
        symplectic[4*i-2, 4*i-1] = -sr * st
        symplectic[4*i-2, 4*i] = sr * ct

        symplectic[4*i-1, 4*i-3] = -sr * ct
        symplectic[4*i-1, 4*i-2] = -sr * st
        symplectic[4*i-1, 4*i-1] = cr

        symplectic[4*i, 4*i-3] = -sr * st
        symplectic[4*i, 4*i-2] = sr * ct
        symplectic[4*i, 4*i] = cr
    end
    return disp, symplectic
end
function _twosqueeze(basis::QuadPairBasis{N}, r::R, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        cr, sr = cosh(r[i]), sinh(r[i])
        ct, st = cos(theta[i]), sin(theta[i])

        symplectic[4*i-3, 4*i-3] = cr
        symplectic[4*i-3, 4*i-1] = -sr * ct
        symplectic[4*i-3, 4*i] = -sr * st

        symplectic[4*i-2, 4*i-2] = cr
        symplectic[4*i-2, 4*i-1] = -sr * st
        symplectic[4*i-2, 4*i] = sr * ct

        symplectic[4*i-1, 4*i-3] = -sr * ct
        symplectic[4*i-1, 4*i-2] = -sr * st
        symplectic[4*i-1, 4*i-1] = cr

        symplectic[4*i, 4*i-3] = -sr * st
        symplectic[4*i, 4*i-2] = sr * ct
        symplectic[4*i, 4*i] = cr
    end
    return disp, symplectic
end
function _twosqueeze(basis::QuadBlockBasis{N}, r::R, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    cr, sr = cosh(r), sinh(r)
    ct, st = cos(theta), sin(theta)
    symplectic = zeros(2*nmodes, 2*nmodes)
    for i in Base.OneTo(Int(nmodes/2))
        symplectic[2*i-1, 2*i-1] = cr
        symplectic[2*i-1, 2*i] = -sr * ct
        symplectic[2*i, 2*i-1] = -sr * ct
        symplectic[2*i, 2*i] = cr

        symplectic[2*i-1, 2*i+nmodes] = -sr * st
        symplectic[2*i, 2*i+nmodes-1] = -sr * st

        symplectic[2*i+nmodes-1, 2*i+nmodes-1] = cr
        symplectic[2*i+nmodes-1, 2*i+nmodes] = sr * ct
        symplectic[2*i+nmodes, 2*i+nmodes-1] = sr * ct
        symplectic[2*i+nmodes, 2*i+nmodes] = cr

        symplectic[2*i+nmodes-1, 2*i] = -sr * st
        symplectic[2*i+nmodes, 2*i-1] = -sr * st
    end
    return disp, symplectic
end
function _twosqueeze(basis::QuadBlockBasis{N}, r::R, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(Int(nmodes/2))
        cr, sr = cosh(r[2*i]), sinh(r[2*i])
        ct, st = cos(theta[2*i]), sin(theta[2*i])

        symplectic[2*i-1, 2*i-1] = cr
        symplectic[2*i-1, 2*i] = -sr * ct
        symplectic[2*i, 2*i-1] = -sr * ct
        symplectic[2*i, 2*i] = cr

        symplectic[2*i-1, 2*i+nmodes] = -sr * st
        symplectic[2*i, 2*i+nmodes-1] = -sr * st

        symplectic[2*i+nmodes-1, 2*i+nmodes-1] = cr
        symplectic[2*i+nmodes-1, 2*i+nmodes] = sr * ct
        symplectic[2*i+nmodes, 2*i+nmodes-1] = sr * ct
        symplectic[2*i+nmodes, 2*i+nmodes] = cr

        symplectic[2*i+nmodes-1, 2*i] = -sr * st
        symplectic[2*i+nmodes, 2*i-1] = -sr * st
    end
    return disp, symplectic
end

"""
    phaseshift([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, theta<:Real)
    phaseshift([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, theta<:Real, noise::Ts)

Gaussian operator that rotates the phase of a given Gaussian mode by `theta`,
as the phase shift operator. The symplectic representation is given by `basis`. Noise can be added to the operation
with `noise`.

## Mathematical description of a phase shift operator

A phase shift operator ``U(\\theta)`` is defined by the operation
``U(\\Theta) = \\exp(-i\\theta\\hat{a}^{\\dagger}\\hat{a})``, where ``\\theta`` is
the phase parameter, and ``\\hat{a}^{\\dagger}`` and ``\\hat{a}`` are the raising
and lowering operators, respectively. The operator 
``U(\\theta)`` is characterized by 
the displacement vector ``\\mathbf{d}`` and symplectic
matrix ``\\mathbf{S}``, expressed respectively as follows:

```math
\\mathbf{d} = \\mathbf{0},
\\quad \\mathbf{S} = \\begin{pmatrix}
                        \\cos(\\theta) & \\sin(\\theta) \\\\
                        -\\sin(\\theta) & \\cos(\\theta) \\\\
                     \\end{pmatrix}.
```

## Example

```jldoctest
julia> phaseshift(QuadPairBasis(1), 3pi/4)
GaussianUnitary for 1 mode in QuadPairBasis representation.
displacement: 2-element Vector{Float64}:
 0.0
 0.0
symplectic: 2×2 Matrix{Float64}:
 -0.707107   0.707107
 -0.707107  -0.707107
```
"""
function phaseshift(::Type{Td}, ::Type{Ts}, basis::SymplecticBasis{N}, theta::R) where {Td,Ts,N<:Int,R}
    disp, symplectic = _phaseshift(basis, theta)
    return GaussianUnitary(basis, Td(disp), Ts(symplectic))
end
phaseshift(::Type{T}, basis::SymplecticBasis{N}, theta::R) where {T,N<:Int,R} = phaseshift(T, T, basis, theta)
function phaseshift(basis::SymplecticBasis{N}, theta::R) where {N<:Int,R}
    disp, symplectic = _phaseshift(basis, theta)
    return GaussianUnitary(basis, disp, symplectic)
end
function _phaseshift(basis::QuadPairBasis{N}, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    ct, st = cos(theta), sin(theta)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        symplectic[2*i-1, 2*i-1] = ct
        symplectic[2*i-1, 2*i] = st
        symplectic[2*i, 2*i-1] = -st
        symplectic[2*i, 2*i] = ct
    end
    return disp, symplectic
end
function _phaseshift(basis::QuadPairBasis{N}, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        ct, st = cos(theta[i]), sin(theta[i])
        symplectic[2*i-1, 2*i-1] = ct
        symplectic[2*i-1, 2*i] = st
        symplectic[2*i, 2*i-1] = -st
        symplectic[2*i, 2*i] = ct
    end
    return disp, symplectic
end
function _phaseshift(basis::QuadBlockBasis{N}, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    ct, st = cos(theta), sin(theta)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        symplectic[i, i] = ct
        symplectic[i, i+nmodes] = st
        symplectic[i+nmodes, i] = -st
        symplectic[i+nmodes, i+nmodes] = ct
    end
    return disp, symplectic
end
function _phaseshift(basis::QuadBlockBasis{N}, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        ct, st = cos(theta[i]), sin(theta[i])
        symplectic[i, i] = ct
        symplectic[i, i+nmodes] = st
        symplectic[i+nmodes, i] = -st
        symplectic[i+nmodes, i+nmodes] = ct
    end
    return disp, symplectic
end

"""
    beamsplitter([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, transmit<:Real)
    beamsplitter([Tm=Vector{Float64}, Ts=Matrix{Float64}], basis::SymplecticBasis, transmit<:Real, noise::Ts)

Gaussian operator that serves as the beam splitter transformation of a
two-mode Gaussian state, known as the beam splitter operator. The symplectic representation is given
by `basis`. The transmittivity of the operator is given by `transmit`. Noise can be added to the operation
with `noise`.

## Mathematical description of a beam splitter operator

A beam splitter operator ``B(\\tau)`` is defined by the operation
``B(\\theta) = \\exp\\left[\\theta(\\hat{a}^{\\dagger}\\hat{b} - \\hat{a}\\hat{b}^{\\dagger})\\right]``,
where ``\\theta`` is defined by ``\\tau = \\cos^2\\theta``, and ``\\hat{a}`` and ``\\hat{b}``
are the annihilation operators of the two modes, respectively. The operator 
``B(\\tau)`` is characterized by 
the displacement vector ``\\mathbf{d}`` and symplectic
matrix ``\\mathbf{S}``, expressed respectively as follows:

```math
\\mathbf{d} = \\mathbf{0},
\\quad \\mathbf{S} = \\begin{pmatrix}
                        \\sqrt{\\tau}\\mathbf{I} & \\sqrt{1-\\tau}\\mathbf{I} \\\\
                        -\\sqrt{1-\\tau}\\mathbf{I} & \\sqrt{\\tau}\\mathbf{I} \\\\
                     \\end{pmatrix}.
```

## Example

```jldoctest
julia> beamsplitter(QuadPairBasis(2), 0.75)
GaussianUnitary for 2 modes in QuadPairBasis representation.
displacement: 4-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
symplectic: 4×4 Matrix{Float64}:
  0.866025   0.0       0.5       0.0
  0.0        0.866025  0.0       0.5
 -0.5        0.0       0.866025  0.0
  0.0       -0.5       0.0       0.866025
```
"""
function beamsplitter(::Type{Td}, ::Type{Ts}, basis::SymplecticBasis{N}, transmit::R) where {Td,Ts,N<:Int,R}
    disp, symplectic = _beamsplitter(basis, transmit)
    return GaussianUnitary(basis, Td(disp), Ts(symplectic))
end
beamsplitter(::Type{T}, basis::SymplecticBasis{N}, transmit::R) where {T,N<:Int,R} = beamsplitter(T, T, basis, transmit)
function beamsplitter(basis::SymplecticBasis{N}, transmit::R) where {N<:Int,R}
    disp, symplectic = _beamsplitter(basis, transmit)
    return GaussianUnitary(basis, disp, symplectic)
end
function _beamsplitter(basis::QuadPairBasis{N}, transmit::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        symplectic[4*i-3, 4*i-3] = a1
        symplectic[4*i-3, 4*i-1] = a2

        symplectic[4*i-2, 4*i-2] = a1
        symplectic[4*i-2, 4*i] = a2
    
        symplectic[4*i-1, 4*i-3] = -a2
        symplectic[4*i-1, 4*i-1] = a1

        symplectic[4*i, 4*i-2] = -a2
        symplectic[4*i, 4*i] = a1
    end
    return disp, symplectic
end
function _beamsplitter(basis::QuadPairBasis{N}, transmit::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        a1, a2 = sqrt(transmit[i]), sqrt(1 - transmit[i])

        symplectic[4*i-3, 4*i-3] = a1
        symplectic[4*i-3, 4*i-1] = a2

        symplectic[4*i-2, 4*i-2] = a1
        symplectic[4*i-2, 4*i] = a2
    
        symplectic[4*i-1, 4*i-3] = -a2
        symplectic[4*i-1, 4*i-1] = a1

        symplectic[4*i, 4*i-2] = -a2
        symplectic[4*i, 4*i] = a1
    end
    return disp, symplectic
end
function _beamsplitter(basis::QuadBlockBasis{N}, transmit::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        symplectic[2*i-1, 2*i-1] = a2
        symplectic[2*i-1, 2*i] = a1
        symplectic[2*i, 2*i-1] = -a1
        symplectic[2*i, 2*i] = a2
    end
    return disp, symplectic
end
function _beamsplitter(basis::QuadBlockBasis{N}, transmit::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    disp = zeros(2*nmodes)
    symplectic = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(Int(nmodes/2))
        a1, a2 = sqrt(transmit[i]), sqrt(1 - transmit[i])

        symplectic[2*i-1, 2*i-1] = a2
        symplectic[2*i-1, 2*i] = a1
        symplectic[2*i, 2*i-1] = -a1
        symplectic[2*i, 2*i] = a2

        symplectic[2*i+nmodes-1, 2*i+nmodes-1] = a2
        symplectic[2*i+nmodes-1, 2*i+nmodes] = a1
        symplectic[2*i+nmodes, 2*i+nmodes-1] = -a1
        symplectic[2*i+nmodes, 2*i+nmodes] = a2
    end
    return disp, symplectic
end

##
# Operations on Gaussian unitaries
##

function tensor(::Type{Td}, ::Type{Ts}, op1::GaussianUnitary, op2::GaussianUnitary) where {Td,Ts}
    typeof(op1.basis) == typeof(op2.basis) || throw(ArgumentError(SYMPLECTIC_ERROR))
    disp, symplectic = _tensor(op1, op2)
    return GaussianUnitary(op1.basis + op2.basis, Td(disp), Ts(symplectic))
end
tensor(::Type{T}, op1::GaussianUnitary, op2::GaussianUnitary) where {T} = tensor(T, T, op1, op2)
function tensor(op1::GaussianUnitary, op2::GaussianUnitary)
    disp, symplectic = _tensor(op1, op2)
    return GaussianUnitary(op1.basis + op2.basis, disp, symplectic)
end
function _tensor(op1::GaussianUnitary{B1,D1,S1}, op2::GaussianChannel{B2,D2,S2}) where {B1<:QuadPairBasis,B2<:QuadPairBasis,D1,D2,S1,S2}
    basis1, basis2 = op1.basis, op2.basis
    nmodes1, nmodes2 = basis1.nmodes, basis2.nmodes
    nmodes = nmodes1 + nmodes2
    block1, block2 = Base.OneTo(2*nmodes1), Base.OneTo(2*nmodes2)
    # initialize direct sum of displacement vectors
    disp1, disp2 = op1.disp, op2.disp
    disp′ = zeros(2*nmodes)
    @inbounds for i in block1
        disp′[i] = disp1[i]
    end
    @inbounds for i in block2
        disp′[i+2*nmodes1] = disp2[i]
    end
    # initialize direct sum of symplectic matrices
    symp1, symp2 = op1.symplectic, op2.symplectic
    symp′ = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in block1, j in block1
        symp′[i,j] = symp1[i,j]
    end
    @inbounds for i in block2, j in block2
        symp′[i+2*nmodes1,j+2*nmodes1] = symp2[i,j]
        symp′[i+2*nmodes1,j+2*nmodes1] = symp2[i,j]
    end
    # extract output array types
    disp′′ = _promote_output_vector(typeof(disp1), typeof(disp2), disp′)
    symp′′ = _promote_output_matrix(typeof(symp1), typeof(symp2), symp′)
    return disp′′, symp′′
end
function _tensor(op1::GaussianUnitary{B1,D1,S1}, op2::GaussianUnitary{B2,D2,S2}) where {B1<:QuadBlockBasis,B2<:QuadBlockBasis,D1,D2,S1,S2}
    basis1, basis2 = op1.basis, op2.basis
    nmodes1, nmodes2 = basis1.nmodes, basis2.nmodes
    nmodes = nmodes1 + nmodes2
    block1, block2 = Base.OneTo(nmodes1), Base.OneTo(nmodes2)
    # initialize direct sum of displacement vectors
    disp1, disp2 = op1.disp, op2.disp
    disp′ = zeros(2*nmodes)
    @inbounds for i in block1
        disp′[i] = disp1[i]
        disp′[i+nmodes] = disp1[i+nmodes1]
    end
    @inbounds for i in block2
        disp′[i+nmodes1] = disp2[i]
        disp′[i+nmodes+nmodes1] = disp2[i+nmodes2]
    end
    # initialize direct sum of symplectic matrices
    symp1, symp2 = op1.symplectic, op2.symplectic
    symp′ = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in block1, j in block1
        symp′[i,j] = symp1[i,j]
        symp′[i,j+nmodes] = symp1[i,j+nmodes1]
        symp′[i+nmodes,j] = symp1[i+nmodes1,j]
        symp′[i+nmodes,j+nmodes] = symp1[i+nmodes1,j+nmodes1]
    end
    @inbounds for i in block2, j in block2
        symp′[i+nmodes1,j+nmodes1] = symp2[i,j]
        symp′[i+nmodes1,j+nmodes+nmodes1] = symp2[i,j+nmodes2]
        symp′[i+nmodes+nmodes1,j+nmodes1] = symp2[i+nmodes2,j]
        symp′[i+nmodes+nmodes1,j+nmodes+nmodes1] = symp2[i+nmodes2,j+nmodes2]
   end
    # extract output array types
    disp′′ = _promote_output_vector(typeof(disp1), typeof(disp2), disp′)
    symp′′ = _promote_output_matrix(typeof(symp1), typeof(symp2), symp′)
    return disp′′, symp′′
end