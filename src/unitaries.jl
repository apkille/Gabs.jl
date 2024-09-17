##
# Predefined Gaussian unitaries
##

"""
    displace([Tm=Vector{Float64}, Ts=Matrix{Float64}], alpha<:Number)
    displace([Tm=Vector{Float64}, Ts=Matrix{Float64}], alpha<:Number, noise::Ts)

Gaussian operator that displaces the vacuum state into a coherent state, known
as the displacement operator. The complex amplitude is given by `alpha`. Noise can
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
julia> displace(1.0+im)
GaussianUnitary
displacement: 2-element Vector{Float64}:
 1.4142135623730951
 1.4142135623730951
symplectic: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
"""
function displace(::Type{Td}, ::Type{Ts}, alpha::N) where {Td,Ts,N<:Number}
    disp = sqrt(2) * Td([real(alpha), imag(alpha)])
    symplectic = Ts(Matrix{Float64}(I, 2, 2))
    return GaussianUnitary(disp, symplectic)
end
displace(::Type{T}, alpha::N) where {T,N<:Number} = displace(T, T, alpha)
function displace(alpha::N) where {N<:Number}
    disp = sqrt(2) * [real(alpha), imag(alpha)]
    symplectic = Matrix{Float64}(I, 2, 2)
    return GaussianUnitary(disp, symplectic)
end

"""
    squeeze([Tm=Vector{Float64}, Ts=Matrix{Float64}], r<:Real, theta<:Real)
    squeeze([Tm=Vector{Float64}, Ts=Matrix{Float64}], r<:Real, theta<:Real, noise::Ts)

Gaussian operator that squeezes the vacuum state into a squeezed state, known
as the squeezing operator. The amplitude and phase squeezing parameters 
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
julia> squeeze(0.25, pi/4)
GaussianUnitary
displacement: 2-element Vector{Float64}:
 0.0
 0.0
symplectic: 2×2 Matrix{Float64}:
 0.215425   0.0451226
 0.0451226  0.30567
```
"""
function squeeze(::Type{Td}, ::Type{Ts}, r::N, theta::N) where {Td,Ts,N<:Real}
    disp = Td(zeros(2))
    cr, sr = cosh(r), sinh(r)
    s = sinh(r) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    symplectic = Ts(s)
    return GaussianUnitary(disp, symplectic)
end
squeeze(::Type{T}, r::N, theta::N) where {T,N<:Real} = squeeze(T, T, r, theta)
function squeeze(r::N, theta::N) where {N<:Real}
    disp = zeros(2)
    cr, sr = cosh(r), sinh(r)
    symplectic = sinh(r) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    return GaussianUnitary(disp, symplectic)
end

"""
    twosqueeze([Tm=Vector{Float64}, Ts=Matrix{Float64}], r<:Real, theta<:Real)
    twosqueeze([Tm=Vector{Float64}, Ts=Matrix{Float64}], r<:Real, theta<:Real, noise::Ts)

Gaussian operator that squeezes a two-mode vacuum state into a two-mode squeezed state,
known as the two-mode squeezing operator. The amplitude and phase squeezing parameters 
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
julia> twosqueeze(0.25, pi/4)
GaussianUnitary
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
function twosqueeze(::Type{Td}, ::Type{Ts}, r::N, theta::N) where {Td,Ts,N<:Real}
    disp = Td(zeros(4))
    v1 = cosh(r) * Matrix{Float64}(I, 2, 2)
    v2 = sinh(r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    symplectic = Ts([v1 -v2; -v2 v1])
    return GaussianUnitary(disp, symplectic)
end
twosqueeze(::Type{T}, r::N, theta::N) where {T,N<:Real} = twosqueeze(T, T, r, theta)
function twosqueeze(r::N, theta::N) where {N<:Real}
    disp = zeros(4)
    v1 = cosh(r) * Matrix{Float64}(I, 2, 2)
    v2 = sinh(r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    symplectic = [v1 -v2; -v2 v1]
    return GaussianUnitary(disp, symplectic)
end

"""
    phaseshift([Tm=Vector{Float64}, Ts=Matrix{Float64}], theta<:Real)
    phaseshift([Tm=Vector{Float64}, Ts=Matrix{Float64}], theta<:Real, noise::Ts)

Gaussian operator that rotates the phase of a given Gaussian mode by `theta`,
as the phase shift operator. Noise can be added to the operation
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
julia> phaseshift(3pi/4)
GaussianUnitary
displacement: 2-element Vector{Float64}:
 0.0
 0.0
symplectic: 2×2 Matrix{Float64}:
 -0.707107   0.707107
 -0.707107  -0.707107
```
"""
function phaseshift(::Type{Td}, ::Type{Ts}, theta::N) where {Td,Ts,N<:Real}
    disp = Td(zeros(2))
    symplectic = Ts([cos(theta) sin(theta); -sin(theta) cos(theta)])
    return GaussianUnitary(disp, symplectic)
end
phaseshift(::Type{T}, theta::N) where {T,N<:Real} = phaseshift(T, T, theta)
function phaseshift(theta::N) where {N<:Real}
    disp = zeros(2)
    symplectic = [cos(theta) sin(theta); -sin(theta) cos(theta)]
    return GaussianUnitary(disp, symplectic)
end

"""
    beamsplitter([Tm=Vector{Float64}, Ts=Matrix{Float64}], transmit<:Real)
    beamsplitter([Tm=Vector{Float64}, Ts=Matrix{Float64}], transmit<:Real, noise::Ts)

Gaussian operator that serves as the beam splitter transformation of a
two-mode Gaussian state, known as the beam splitter operator. The transmittivity
of the operator is given by `transmit`. Noise can be added to the operation
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
julia> beamsplitter(0.75)
GaussianUnitary
displacement: 4-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
symplectic: 4×4 Matrix{Float64}:
  0.866025   0.0       0.5       0.0
  0.0        0.866025  0.0       0.5
 -0.5       -0.0       0.866025  0.0
 -0.0       -0.5       0.0       0.866025
```
"""
function beamsplitter(::Type{Td}, ::Type{Ts}, transmit::N) where {Td,Ts,N<:Real}
    disp = Td(zeros(4))
    I2 = Matrix{Float64}(I, 2, 2)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    symplectic = Ts([a1*I2 a2*I2; -a2*I2 a1*I2])
    return GaussianUnitary(disp, symplectic)
end
beamsplitter(::Type{T}, transmit::N) where {T,N<:Real} = beamsplitter(T, T, transmit)
function beamsplitter(transmit::N) where {N<:Real}
    disp = zeros(4)
    I2 = Matrix{Float64}(I, 2, 2)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    symplectic = [a1*I2 a2*I2; -a2*I2 a1*I2]
    return GaussianUnitary(disp, symplectic)
end

##
# Operations on Gaussian unitaries
##

function tensor(::Type{Td}, ::Type{Ts}, op1::GaussianUnitary, op2::GaussianUnitary) where {Td,Ts}
    disp′, symplectic′ = _tensor_fields(op1, op2)
    return GaussianUnitary(Td(disp′), Ts(symplectic′))
end
tensor(::Type{T}, op1::GaussianUnitary, op2::GaussianUnitary) where {T} = tensor(T, T, op1, op2)
function tensor(op1::GaussianUnitary, op2::GaussianUnitary)
    disp′, symplectic′ = _tensor_fields(op1, op2)
    return GaussianUnitary(disp′, symplectic′)
end
function _tensor_fields(op1::GaussianUnitary, op2::GaussianUnitary)
    disp1, disp2 = op1.disp, op2.disp
    length1, length2 = length(disp1), length(disp2)
    slengths = length1 + length2
    symp1, symp2 = op1.symplectic, op2.symplectic
    # initialize direct sum of displacement vectors
    disp′ = zeros(slengths)
    @inbounds for i in eachindex(disp1)
        disp′[i] = disp1[i]
    end
    @inbounds for i in eachindex(disp2)
        disp′[i+length1] = disp2[i]
    end
    # initialize direct sum of symplectic matrices
    symplectic′ = zeros(slengths, slengths)
    axes1 = axes(symp1)
    @inbounds for i in axes1[1], j in axes1[2]
        symplectic′[i,j] = symp1[i,j]
    end
    axes2 = axes(symp2)
    @inbounds for i in axes2[1], j in axes2[2]
        symplectic′[i+length1,j+length1] = symp2[i,j]
    end
    return disp′, symplectic′
end