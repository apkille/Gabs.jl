##
# Predefined Gaussian states
##

"""
    vacuumstate([Tm=Vector{Float64}, Tc=Matrix{Float64}], basis::SymplecticBasis)

Gaussian state with zero photons, known as the vacuum state. The symplectic representation is defined by `basis`.

## Mathematical description of a vacuum state

A vacuum state ``|0\\rangle`` is characterized by the mean vector ``\\mathbf{\\bar{x}}`` and covariance
matrix ``\\mathbf{V}``, expressed respectively as follows:

```math
\\mathbf{\\bar{x}} = \\mathbf{0}, \\quad \\mathbf{V} = \\mathbf{I}.
```

## Example

```jldoctest
julia> vacuumstate(QuadPairBasis(1))
GaussianState for 1 mode in QuadPairBasis representation.
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
"""
function vacuumstate(::Type{Tm}, ::Type{Tc}, basis::SymplecticBasis{N}) where {Tm,Tc,N<:Int}
    mean, covar = _vacuumstate(basis)
    return GaussianState(basis, Tm(mean), Tc(covar))
end
vacuumstate(::Type{T}, basis::SymplecticBasis{N}) where {T,N<:Int} = vacuumstate(T, T, basis)
function vacuumstate(basis::SymplecticBasis{N}) where {N<:Int}
    mean, covar = _vacuumstate(basis)
    return GaussianState(basis, mean, covar)
end
function _vacuumstate(basis::SymplecticBasis{N}) where {N<:Int}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    covar = Matrix{Float64}(I, 2*nmodes, 2*nmodes)
    return mean, covar
end

"""
    thermalstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] basis::SymplecticBasis, photons<:Int)

Gaussian state at thermal equilibrium, known as the thermal state. The symplectic representation
is defined by `basis`. The mean photon number of the state is given by `photons`.

## Mathematical description of a thermal state

A thermal state ``|\\bar{n}\\rangle``, where ``\\bar{n}`` is the mean number of photons,
is characterized by the mean vector ``\\mathbf{\\bar{x}}`` and covariance
matrix ``\\mathbf{V}``, expressed respectively as follows:

```math
\\mathbf{\\bar{x}} = \\mathbf{0}, \\quad \\mathbf{V} = \\left(\\bar{n} + \\frac{1}{2}\\right)\\mathbf{I}.
```

## Example

```jldoctest
julia> thermalstate(QuadPairBasis(1), 4)
GaussianState for 1 mode in QuadPairBasis representation.
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 4.5  0.0
 0.0  4.5
```
"""
function thermalstate(::Type{Tm}, ::Type{Tc}, basis::SymplecticBasis{N}, photons::P) where {Tm,Tc,N<:Int,P}
    mean, covar = _thermalstate(basis, photons)
    return GaussianState(basis, Tm(mean), Tc(covar))
end
thermalstate(::Type{T}, basis::SymplecticBasis{N}, photons::P) where {T,N<:Int,P} = thermalstate(T, T, basis, photons)
function thermalstate(basis::SymplecticBasis{N}, photons::P) where {N<:Int,P}
    mean, covar = _thermalstate(basis, photons)
    return GaussianState(basis, mean, covar)
end
function _thermalstate(basis::Union{QuadPairBasis{N},QuadBlockBasis{N}}, photons::P) where {N<:Int,P<:Int}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    covar = Matrix{Float64}((photons + 1/2) * I, 2*nmodes, 2*nmodes)
    return mean, covar
end
function _thermalstate(basis::Union{QuadPairBasis{N},QuadBlockBasis{N}}, photons::P) where {N<:Int,P<:Vector}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    covar = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        val = photons[i] + 1/2
        covar[2*i-1, 2*i-1] = val
        covar[2*i, 2*i] = val
    end
    return mean, covar
end

"""
    coherentstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] basis::SymplecticBasis, alpha<:Number)

Gaussian state that is the quantum analogue of a monochromatic electromagnetic field, known
as the coherent state. The symplectic representation is defined by `basis`.
The complex amplitude of the state is given by `alpha`.

## Mathematical description of a coherent state

A coherent state ``|\\alpha\\rangle``, where ``\\alpha`` is the complex amplitude,
is characterized by the mean vector ``\\mathbf{\\bar{x}}`` and covariance
matrix ``\\mathbf{V}``, expressed respectively as follows:

```math
\\mathbf{\\bar{x}} = \\sqrt{2}\\left(\\text{Re}(\\alpha), \\text{Im}(\\alpha)\\right)^{\\text{T}},
\\quad \\mathbf{V} = \\mathbf{I}.
```

## Example

```jldoctest
julia> coherentstate(QuadPairBasis(1), 1.0+im)
GaussianState for 1 mode in QuadPairBasis representation.
mean: 2-element Vector{Float64}:
 1.4142135623730951
 1.4142135623730951
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
"""
function coherentstate(::Type{Tm}, ::Type{Tc}, basis::SymplecticBasis{N}, alpha::A) where {Tm,Tc,N<:Int,A}
    mean, covar = _coherentstate(basis, alpha)
    return GaussianState(basis, Tm(mean), Tc(covar))
end
coherentstate(::Type{T}, basis::SymplecticBasis{N}, alpha::A) where {T,N<:Int,A} = coherentstate(T, T, basis, alpha)
function coherentstate(basis::SymplecticBasis{N}, alpha::A) where {N<:Int,A}
    mean, covar = _coherentstate(basis, alpha)
    return GaussianState(basis, mean, covar)
end
function _coherentstate(basis::QuadPairBasis{N}, alpha::A) where {N<:Int,A<:Number}
    nmodes = basis.nmodes
    mean = repeat([sqrt(2)*real(alpha), sqrt(2)*imag(alpha)], nmodes)
    covar = Matrix{Float64}(I, 2*nmodes, 2*nmodes)
    return mean, covar
end
function _coherentstate(basis::QuadPairBasis{N}, alpha::A) where {N<:Int,A<:Vector}
    nmodes = basis.nmodes
    mean = sqrt(2) * reinterpret(Float64, alpha)
    covar = Matrix{Float64}(I, 2*nmodes, 2*nmodes)
    return mean, covar
end
function _coherentstate(basis::QuadBlockBasis{N}, alpha::A) where {N<:Int,A<:Number}
    nmodes = basis.nmodes
    mean = repeat([sqrt(2)*real(alpha), sqrt(2)*imag(alpha)], inner = nmodes)
    covar = Matrix{Float64}(I, 2*nmodes, 2*nmodes)
    return mean, covar
end
function _coherentstate(basis::QuadBlockBasis{N}, alpha::A) where {N<:Int,A<:Vector}
    nmodes = basis.nmodes
    re = reinterpret(Float64, alpha)
    mean = vcat(@view(re[1:2:end]), @view(re[2:2:end]))
    mean .*= sqrt(2)
    covar = Matrix{Float64}(I, 2*nmodes, 2*nmodes)
    return mean, covar
end

"""
    squeezedstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] basis::SymplecticBasis, r<:Real, theta<:Real)

Gaussian state with quantum uncertainty in its phase and amplitude, known as
the squeezed state. The symplectic representation is defined by `basis`. The amplitude and phase squeezing parameters are given by `r`
and `theta`, respectively.

## Mathematical description of a squeezed state

A squeezed state ``|r, \\theta\\rangle``, where ``r`` is the amplitude squeezing
parameter and ``\\theta`` is the phase squeezing parameter,
is characterized by the mean vector ``\\mathbf{\\bar{x}}`` and covariance
matrix ``\\mathbf{V}``, expressed respectively as follows:

```math
\\mathbf{\\bar{x}} = \\mathbf{0},
\\quad \\mathbf{V} = \\frac{1}{2}\\left(\\cosh(2r)\\mathbf{I} - \\sinh(2r)\\mathbf{R}(\\theta)\\right),
```

where ``\\mathbf{R}(\\theta)`` is the rotation matrix.

## Example

```jldoctest
julia> squeezedstate(QuadPairBasis(1), 0.5, pi/4)
GaussianState for 1 mode in QuadPairBasis representation.
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 0.356044  0.415496
 0.415496  1.18704
```
"""
function squeezedstate(::Type{Tm}, ::Type{Tc}, basis::SymplecticBasis{N}, r::R, theta::R) where {Tm,Tc,N<:Int,R}
    mean, covar = _squeezedstate(basis, r, theta)
    return GaussianState(basis, Tm(mean), Tc(covar))
end
squeezedstate(::Type{T}, basis::SymplecticBasis{N}, r::R, theta::R) where {T,N<:Int,R} = squeezedstate(T, T, basis, r, theta)
function squeezedstate(basis::SymplecticBasis{N}, r::R, theta::R) where {N<:Int,R}
    mean, covar = _squeezedstate(basis, r, theta)
    return GaussianState(basis, mean, covar)
end
function _squeezedstate(basis::QuadPairBasis{N}, r::R, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    covar = zeros(2*nmodes, 2*nmodes)
    cr, sr = cosh(2*r), sinh(2*r)
    ct, st = cos(theta), sin(theta)
    @inbounds for i in Base.OneTo(nmodes)
        covar[2*i-1, 2*i-1] = (1/2) * (cr - sr*ct)
        covar[2*i-1, 2*i] = (1/2) * sr * st
        covar[2*i, 2*i-1] = (1/2) * sr * st
        covar[2*i, 2*i] = (1/2) * (cr + sr*ct)
    end
    return mean, covar
end
function _squeezedstate(basis::QuadPairBasis{N}, r::R, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    covar = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        cr, sr = cosh(2*r[i]), sinh(2*r[i])
        ct, st = cos(theta[i]), sin(theta[i])
        covar[2*i-1, 2*i-1] = (1/2) * (cr - sr*ct)
        covar[2*i-1, 2*i] = (1/2) * sr * st
        covar[2*i, 2*i-1] = (1/2) * sr * st
        covar[2*i, 2*i] = (1/2) * (cr + sr*ct)
    end
    return mean, covar
end
function _squeezedstate(basis::QuadBlockBasis{N}, r::R, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    covar = zeros(2*nmodes, 2*nmodes)
    cr, sr = cosh(2*r), sinh(2*r)
    ct, st = cos(theta), sin(theta)
    @inbounds for i in Base.OneTo(nmodes)
        covar[i, i] = (1/2) * (cr - sr*ct)
        covar[i, i+nmodes] = (1/2) * sr * st
        covar[i+nmodes, i] = (1/2) * sr * st
        covar[i+nmodes, i+nmodes] = (1/2) * (cr + sr*ct)
    end
    return mean, covar
end
function _squeezedstate(basis::QuadBlockBasis{N}, r::R, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    covar = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        cr, sr = cosh(2*r[i]), sinh(2*r[i])
        ct, st = cos(theta[i]), sin(theta[i])
        covar[i, i] = (1/2) * (cr - sr*ct)
        covar[i, i+nmodes] = (1/2) * sr * st
        covar[i+nmodes, i] = (1/2) * sr * st
        covar[i+nmodes, i+nmodes] = (1/2) * (cr + sr*ct)
    end
    return mean, covar
end

"""
    eprstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] basis::SymplecticBasis, r<:Real, theta<:Real)

Gaussian state that is a two-mode squeezed state, known as the Einstein-Podolski-Rosen (EPR) state. The symplectic
representation is defined by `basis`. The amplitude and phase squeezing parameters are given by `r` and `theta`, respectively.

## Mathematical description of an EPR state

An EPR state ``|r, \\theta\\rangle_{\\text{EPR}}``, where ``r`` is the amplitude squeezing
parameter and ``\\theta`` is the phase squeezing parameter,
is characterized by the mean vector ``\\mathbf{\\bar{x}}`` and covariance
matrix ``\\mathbf{V}``, expressed respectively as follows:

```math
\\mathbf{\\bar{x}} = \\mathbf{0},
\\quad \\mathbf{V} = \\frac{1}{2}\\begin{pmatrix}
                                    \\cosh(2r)\\mathbf{I} & -\\sinh(2r)\\mathbf{R}(\\theta) \\\\
                                    -\\sinh(2r)\\mathbf{R}(\\theta) & \\cosh(2r)\\mathbf{I} \\\\
                                    \\end{pmatrix},
```

where ``\\mathbf{R}(\\theta)`` is the rotation matrix.

## Example

```jldoctest
julia> eprstate(QuadPairBasis(2), 0.5, pi/4)
GaussianState for 2 modes in QuadPairBasis representation.
mean: 4-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
covariance: 4×4 Matrix{Float64}:
 0.77154    0.0       0.415496   0.415496
 0.0        0.77154   0.415496  -0.415496
 0.415496   0.415496  0.77154    0.0
 0.415496  -0.415496  0.0        0.77154
```
"""
function eprstate(::Type{Tm}, ::Type{Tc}, basis::SymplecticBasis{N}, r::R, theta::R) where {Tm,Tc,N<:Int,R}
    mean, covar = _eprstate(basis, r, theta)
    return GaussianState(basis, Tm(mean), Tc(covar))
end
eprstate(::Type{T}, basis::SymplecticBasis{N}, r::R, theta::R) where {T,N<:Int,R} = eprstate(T, T, basis, r, theta)
function eprstate(basis::SymplecticBasis{N}, r::R, theta::R) where {N<:Int,R}
    mean, covar = _eprstate(basis, r, theta)
    return GaussianState(basis, mean, covar)
end
function _eprstate(basis::QuadPairBasis{N}, r::R, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    cr, sr = (1/2)*cosh(2*r), (1/2)*sinh(2*r)
    ct, st = cos(theta), sin(theta)
    covar = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(Int(nmodes/2))
        covar[4*i-3, 4*i-3] = cr
        covar[4*i-3, 4*i-1] = sr * ct
        covar[4*i-3, 4*i] = sr * st

        covar[4*i-2, 4*i-2] = cr
        covar[4*i-2, 4*i-1] = sr * st
        covar[4*i-2, 4*i] = -sr * ct

        covar[4*i-1, 4*i-3] = sr * ct
        covar[4*i-1, 4*i-2] = sr * st
        covar[4*i-1, 4*i-1] = cr

        covar[4*i, 4*i-3] = sr * st
        covar[4*i, 4*i-2] = -sr * ct
        covar[4*i, 4*i] = cr
    end
    return mean, covar
end
function _eprstate(basis::QuadPairBasis{N}, r::R, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    covar = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(Int(nmodes/2))
        cr, sr = (1/2)*cosh(2*r[i]), (1/2)*sinh(2*r[i])
        ct, st = cos(theta[i]), sin(theta[i])

        covar[4*i-3, 4*i-3] = cr
        covar[4*i-3, 4*i-1] = sr * ct
        covar[4*i-3, 4*i] = sr * st

        covar[4*i-2, 4*i-2] = cr
        covar[4*i-2, 4*i-1] = sr * st
        covar[4*i-2, 4*i] = -sr * ct

        covar[4*i-1, 4*i-3] = sr * ct
        covar[4*i-1, 4*i-2] = sr * st
        covar[4*i-1, 4*i-1] = cr

        covar[4*i, 4*i-3] = sr * st
        covar[4*i, 4*i-2] = -sr * ct
        covar[4*i, 4*i] = cr
    end
    return mean, covar
end
function _eprstate(basis::QuadBlockBasis{N}, r::R, theta::R) where {N<:Int,R<:Real}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    cr, sr = (1/2)*cosh(2*r), (1/2)*sinh(2*r)
    ct, st = cos(theta), sin(theta)
    covar = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(Int(nmodes/2))
        covar[2*i-1, 2*i-1] = cr
        covar[2*i-1, 2*i] = -sr * ct
        covar[2*i, 2*i-1] = -sr * ct
        covar[2*i, 2*i] = cr

        covar[2*i-1, 2*i+nmodes] = -sr * st
        covar[2*i, 2*i+nmodes-1] = -sr * st

        covar[2*i+nmodes-1, 2*i+nmodes-1] = cr
        covar[2*i+nmodes-1, 2*i+nmodes] = sr * ct
        covar[2*i+nmodes, 2*i+nmodes-1] = sr * ct
        covar[2*i+nmodes, 2*i+nmodes] = cr

        covar[2*i+nmodes-1, 2*i] = -sr * st
        covar[2*i+nmodes, 2*i-1] = -sr * st
    end
    return mean, covar
end
function _eprstate(basis::QuadBlockBasis{N}, r::R, theta::R) where {N<:Int,R<:Vector}
    nmodes = basis.nmodes
    mean = zeros(2*nmodes)
    covar = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(Int(nmodes/2))
        cr, sr = (1/2)*cosh(2*r[i]), (1/2)*sinh(2*r[i])
        ct, st = cos(theta[i]), sin(theta[i])

        covar[2*i-1, 2*i-1] = cr
        covar[2*i-1, 2*i] = -sr * ct
        covar[2*i, 2*i-1] = -sr * ct
        covar[2*i, 2*i] = cr

        covar[2*i-1, 2*i+nmodes] = -sr * st
        covar[2*i, 2*i+nmodes-1] = -sr * st

        covar[2*i+nmodes-1, 2*i+nmodes-1] = cr
        covar[2*i+nmodes-1, 2*i+nmodes] = sr * ct
        covar[2*i+nmodes, 2*i+nmodes-1] = sr * ct
        covar[2*i+nmodes, 2*i+nmodes] = cr

        covar[2*i+nmodes-1, 2*i] = -sr * st
        covar[2*i+nmodes, 2*i-1] = -sr * st
    end
    return mean, covar
end

##
# Operations on Gaussian states
##

"""
    tensor(state1::GaussianState, state2::GaussianState)

tensor product of Gaussian states, which can also be called with `⊗`.

## Example
```jldoctest
julia> basis = QuadPairBasis(1);

julia> coherentstate(basis, 1.0+im) ⊗ thermalstate(basis, 2)
GaussianState for 2 modes in QuadPairBasis representation.
mean: 4-element Vector{Float64}:
 1.4142135623730951
 1.4142135623730951
 0.0
 0.0
covariance: 4×4 Matrix{Float64}:
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  2.5  0.0
 0.0  0.0  0.0  2.5
```
"""
function tensor(::Type{Tm}, ::Type{Tc}, state1::GaussianState, state2::GaussianState) where {Tm,Tc}
    typeof(state1.basis) == typeof(state2.basis) || throw(ArgumentError(SYMPLECTIC_ERROR))
    mean, covar = _tensor(state1, state2)
    return GaussianState(state1.basis + state2.basis, Tm(mean), Tc(covar))
end
tensor(::Type{T}, state1::GaussianState, state2::GaussianState) where {T} = tensor(T, T, state1, state2)
function tensor(state1::GaussianState, state2::GaussianState)
    typeof(state1.basis) == typeof(state2.basis) || throw(ArgumentError(SYMPLECTIC_ERROR))
    mean, covar = _tensor(state1, state2)
    return GaussianState(state1.basis + state2.basis, mean, covar)
end
function _tensor(state1::GaussianState{B1,M1,V1}, state2::GaussianState{B2,M2,V2}) where {B1<:QuadPairBasis,B2<:QuadPairBasis,M1,M2,V1,V2}
    mean1, mean2 = state1.mean, state2.mean
    basis1, basis2 = state1.basis, state2.basis
    nmodes1, nmodes2 = basis1.nmodes, basis2.nmodes
    nmodes = nmodes1 + nmodes2
    block1, block2 = Base.OneTo(2*nmodes1), Base.OneTo(2*nmodes2)
    # initialize direct sum of mean vectors
    mean′ = zeros(2*nmodes)
    @inbounds for i in block1
        mean′[i] = mean1[i]
    end
    @inbounds for i in block2
        mean′[i+2*nmodes1] = mean2[i]
    end
    # initialize direct sum of covariance matrices
    covar′ = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in block1, j in block1
        covar′[i,j] = covar1[i,j]
    end
    @inbounds for i in block2, j in block2
        covar′[i+2*nmodes,j+2*nmodes] = covar2[i,j]
    end
    # extract output array types
    mean′′ = _promote_output_vector(typeof(mean1), typeof(mean2), mean′)
    covar′′ = _promote_output_matrix(typeof(covar1), typeof(covar2), covar′)
    return mean′′, covar′′
end
function _tensor(state1::GaussianState{B1,M1,V1}, state2::GaussianState{B2,M2,V2}) where {B1<:QuadBlockBasis,B2<:QuadBlockBasis,M1,M2,V1,V2}
    mean1, mean2 = state1.mean, state2.mean
    basis1, basis2 = state1.basis, state2.basis
    nmodes1, nmodes2 = basis1.nmodes, basis2.nmodes
    nmodes = nmodes1 + nmodes2
    block1, block2 = Base.OneTo(nmodes1), Base.OneTo(nmodes2)
    # initialize direct sum of mean vectors
    mean1, mean2 = state1.mean, state2.mean
    mean′ = zeros(2*nmodes)
    @inbounds for i in block1
        mean′[i] = mean1[i]
        mean′[i+nmodes] = mean1[i+nmodes1]
    end
    @inbounds for i in block2
        mean′[i+nmodes1] = mean2[i]
        mean′[i+nmodes+nmodes1] = mean2[i+nmodes2]
    end
    # initialize direct sum of covariance matrices
    covar1, covar2 = state1.covar, state2.covar
    covar′ = zeros(2*nmodes, 2*nmodes)
    covar1, covar2 = state1.covar, state2.covar
    covar′ = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in block1, j in block1
        covar′[i,j] = covar1[i,j]
        covar′[i,j+nmodes] = covar1[i,j+nmodes1]
        covar′[i+nmodes,j] = covar1[i+nmodes1,j]
        covar′[i+nmodes,j+nmodes] = covar1[i+nmodes1,j+nmodes1]
    end
    @inbounds for i in block2, j in block2
        covar′[i+nmodes1,j+nmodes1] = covar2[i,j]
        covar′[i+nmodes1,j+nmodes+nmodes1] = covar2[i,j+nmodes2]
        covar′[i+nmodes+nmodes1,j+nmodes1] = covar2[i+nmodes2,j]
        covar′[i+nmodes+nmodes1,j+nmodes+nmodes1] = covar2[i+nmodes2,j+nmodes2]
    end
    # extract output array types
    mean′′ = _promote_output_vector(typeof(mean1), typeof(mean2), mean′)
    covar′′ = _promote_output_matrix(typeof(covar1), typeof(covar2), covar′)
    return mean′′, covar′′
end

"""
    ptrace([Tm=Vector{Float64}, Tc=Matrix{Float64},] state::GaussianState, idx<:Int)
    ptrace([Tm=Vector{Float64}, Tc=Matrix{Float64},] state::GaussianState, indices<:AbstractVector)

Partial trace of a Gaussian state over a subsystem indicated by `idx`, or multiple subsystems
indicated by `indices`.

## Example
```jldoctest
julia> basis = QuadPairBasis(1);

julia> state = coherentstate(basis, 1.0+im) ⊗ thermalstate(basis, 2) ⊗ squeezedstate(basis, 3.0, pi/4)
GaussianState for 3 modes in QuadPairBasis representation.
mean: 6-element Vector{Float64}:
 1.4142135623730951
 1.4142135623730951
 0.0
 0.0
 0.0
 0.0
covariance: 6×6 Matrix{Float64}:
 1.0  0.0  0.0  0.0   0.0       0.0
 0.0  1.0  0.0  0.0   0.0       0.0
 0.0  0.0  2.5  0.0   0.0       0.0
 0.0  0.0  0.0  2.5   0.0       0.0
 0.0  0.0  0.0  0.0  29.5414   71.3164
 0.0  0.0  0.0  0.0  71.3164  172.174

julia> ptrace(state, 2)
GaussianState for 1 mode in QuadPairBasis representation.
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 2.5  0.0
 0.0  2.5

julia> ptrace(state, [1, 3])
GaussianState for 2 modes in QuadPairBasis representation.
mean: 4-element Vector{Float64}:
 1.4142135623730951
 1.4142135623730951
 0.0
 0.0
covariance: 4×4 Matrix{Float64}:
 1.0  0.0   0.0       0.0
 0.0  1.0   0.0       0.0
 0.0  0.0  29.5414   71.3164
 0.0  0.0  71.3164  172.174
```
"""
function ptrace(::Type{Tm}, ::Type{Tc}, state::GaussianState, idx::N) where {Tm,Tc,N<:Int}
    mean′, covar′ = _ptrace(state, idx)
    return GaussianState(typeof(state.basis)(1), Tm(mean′), Tc(covar′))
end
ptrace(::Type{T}, state::GaussianState, idx::N) where {T,N<:Int} = ptrace(T, T, state, idx)
function ptrace(state::GaussianState, idx::N) where {N<:Int}
    mean′, covar′ = _ptrace(state, idx)
    return GaussianState(typeof(state.basis)(1), mean′, covar′)
end
function ptrace(::Type{Tm}, ::Type{Tc}, state::GaussianState, indices::N) where {Tm,Tc,N<:AbstractVector}
    mean, covar = _ptrace(state, indices)
    return GaussianState(typeof(state.basis)(length(indices)), Tm(mean), Tc(covar))
end
ptrace(::Type{T}, state::GaussianState, indices::N) where {T,N<:AbstractVector} = ptrace(T, T, state, indices)
function ptrace(state::GaussianState, indices::T) where {T<:AbstractVector}
    mean, covar = _ptrace(state, indices)
    return GaussianState(typeof(state.basis)(length(indices)), mean, covar)
end
function _ptrace(state::GaussianState{B,M,V}, idx::T) where {B<:QuadPairBasis,M,V,T<:Int}
    idxV = 2*idx-1:(2*idx)
    mean = state.mean
    covar = state.covar
    # initialize partial trace of mean vector
    mean′ = mean[idxV]
    # initialize partial trace of covariance matrix
    covar′ = covar[idxV, idxV]
    # extract output array types
    mean′′ = _promote_output_vector(typeof(mean), mean′, 2)
    covar′′ = _promote_output_matrix(typeof(covar), covar′, 2)
    return mean′′, covar′′
end
function _ptrace(state::GaussianState{B,M,V}, indices::T) where {B<:QuadPairBasis,M,V,T<:AbstractVector}
    idxlength = length(indices)
    mean = state.mean
    covar = state.covar
    # initialize partial trace of mean vector
    mean′ = zeros(2*idxlength)
    @inbounds for i in eachindex(indices)
        idx = indices[i]
        mean′[2*i-1] = mean[2*idx-1]
        mean′[2*i] = mean[2*idx]
    end
    # initialize partial trace of covariance matrix
    covar′ = zeros(2*idxlength, 2*idxlength)
    @inbounds for i in eachindex(indices)
        idx = indices[i]
        covar′[2*i-1, 2*i-1] = covar[2*idx-1, 2*idx-1]
        covar′[2*i-1, 2*i] = covar[2*idx-1, 2*idx]
        covar′[2*i, 2*i-1] = covar[2*idx, 2*idx-1]
        covar′[2*i, 2*i] = covar[2*idx, 2*idx]
    end 
    mean′′ = _promote_output_vector(typeof(mean), mean′, 2*idxlength)
    covar′′ = _promote_output_matrix(typeof(covar), covar′, 2*idxlength)
    return mean′′, covar′′
end
function _ptrace(state::GaussianState{B,M,V}, idx::T) where {B<:QuadBlockBasis,M,V,T<:Int}
    basis = state.basis
    nmodes = basis.nmodes
    mean = state.mean
    covar = state.covar
    # initialize partial trace of mean vector
    mean′ = [mean[idx], mean[idx+nmodes]]
    # initialize partial trace of covariance matrix
    covar′ = [covar[idx,idx] covar[idx,idx+nmodes]; covar[idx+nmodes,idx] covar[idx+nmodes,idx+nmodes]]
    # extract output array types
    mean′′ = _promote_output_vector(typeof(mean), mean′, 2)
    covar′′ = _promote_output_matrix(typeof(covar), covar′, 2)
    return mean′′, covar′′
end
function _ptrace(state::GaussianState{B,M,V}, indices::T) where {B<:QuadBlockBasis,M,V,T<:AbstractVector}
    basis = state.basis
    nmodes = basis.nmodes
    idxlength = length(indices)
    mean = state.mean
    covar = state.covar
    # initialize partial trace of mean vector
    mean′ = zeros(2*idxlength)
    @inbounds for i in eachindex(indices)
        idx = indices[i]
        mean′[i] = mean[idx]
        mean′[i+idxlength] = mean[idx+nmodes]
    end
    # initialize partial trace of covariance matrix
    covar′ = zeros(2*idxlength, 2*idxlength)
    @inbounds for i in Base.OneTo(idxlength)
        idx = indices[i]
        @inbounds for j in i:idxlength
            otheridx = indices[j]
            covar′[i,j] = covar[idx,otheridx]
            covar′[j,i] = covar[otheridx,idx]
            covar′[i+idxlength,j] = covar[idx+nmodes,otheridx]
            covar′[j,i+idxlength] = covar[otheridx,idx+nmodes]
            covar′[i+idxlength,j+idxlength] = covar[idx+nmodes, otheridx+nmodes]
            covar′[j+idxlength,i+idxlength] = covar[otheridx+nmodes, idx+nmodes]
        end
    end 
    mean′′ = _promote_output_vector(typeof(mean), mean′, 2*idxlength)
    covar′′ = _promote_output_matrix(typeof(covar), covar′, 2*idxlength)
    return mean′′, covar′′
end