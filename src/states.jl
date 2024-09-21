##
# Predefined Gaussian states
##

"""
    vacuumstate([Tm=Vector{Float64}, Tc=Matrix{Float64}])

Gaussian state with zero photons, known as the vacuum state.

## Mathematical description of a vacuum state

A vacuum state ``|0\\rangle`` is characterized by the mean vector ``\\mathbf{\\bar{x}}`` and covariance
matrix ``\\mathbf{V}``, expressed respectively as follows:

```math
\\mathbf{\\bar{x}} = \\mathbf{0}, \\quad \\mathbf{V} = \\mathbf{I}.
```

## Example

```jldoctest
julia> vacuumstate()
GaussianState for 1 mode.
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
"""
function vacuumstate(::Type{Tm}, ::Type{Tc}) where {Tm,Tc}
    mean = Tm(zeros(2))
    covar = Tc(Matrix{Float64}(I, 2, 2))
    return GaussianState(mean, covar, 1)
end
vacuumstate(::Type{T}) where {T} = vacuumstate(T, T)
function vacuumstate()
    mean = zeros(2)
    covar = Matrix{Float64}(I, 2, 2)
    return GaussianState(mean, covar, 1)
end

"""
    thermalstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] photons<:Int)

Gaussian state at thermal equilibrium, known as the thermal state. The mean photon number of the
state is given by `photons`.

## Mathematical description of a thermal state

A thermal state ``|\\bar{n}\\rangle``, where ``\\bar{n}`` is the mean number of photons,
is characterized by the mean vector ``\\mathbf{\\bar{x}}`` and covariance
matrix ``\\mathbf{V}``, expressed respectively as follows:

```math
\\mathbf{\\bar{x}} = \\mathbf{0}, \\quad \\mathbf{V} = \\left(\\bar{n} + \\frac{1}{2}\\right)\\mathbf{I}.
```

## Example

```jldoctest
julia> thermalstate(4)
GaussianState for 1 mode.
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 4.5  0.0
 0.0  4.5
```
"""
function thermalstate(::Type{Tm}, ::Type{Tc}, photons::N) where {Tm,Tc,N<:Int}
    mean = Tm(zeros(2))
    covar = (photons + 1/2) * Tc(Matrix{Float64}(I, 2, 2))
    return GaussianState(mean, covar, 1)
end
thermalstate(::Type{T}, photons::N) where {T, N<:Int} = thermalstate(T, T, photons)
function thermalstate(photons::N) where {N<:Int}
    mean = zeros(2)
    covar = (photons + 1/2) * Matrix{Float64}(I, 2, 2)
    return GaussianState(mean, covar, 1)
end

"""
    coherentstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] alpha<:Number)

Gaussian state that is the quantum analogue of a monochromatic electromagnetic field, known
as the coherent state. The complex amplitude of the state is given by `alpha`.

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
function coherentstate(::Type{Tm}, ::Type{Tc}, alpha::N) where {Tm,Tc,N<:Number}
    mean = sqrt(2) * Tm([real(alpha), imag(alpha)])
    covar = Tc(Matrix{Float64}(I, 2, 2))
    return GaussianState(mean, covar, 1)
end
coherentstate(::Type{T}, alpha::N) where {T, N<:Number} = coherentstate(T, T, alpha)
function coherentstate(alpha::N) where {N<:Number}
    mean = sqrt(2) * [real(alpha), imag(alpha)]
    covar = Matrix{Float64}(I, 2, 2)
    return GaussianState(mean, covar, 1)
end

"""
    squeezedstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] r<:Real, theta<:Real)

Gaussian state with quantum uncertainty in its phase and amplitude, known as
the squeezed state. The amplitude and phase squeezing parameters are given by `r`
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
julia> squeezedstate(0.5, pi/4)
GaussianState for 1 mode.
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 0.356044  0.415496
 0.415496  1.18704
```
"""
function squeezedstate(::Type{Tm}, ::Type{Tc}, r::N, theta::N) where {Tm,Tc,N<:Real}
    mean = Tm(zeros(2))
    cr, sr = cosh(2*r), sinh(2*r)
    v = (1/2) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    covar = Tc(v)
    return GaussianState(mean, covar, 1)
end
squeezedstate(::Type{T}, r::N, theta::N) where {T,N<:Real} = squeezedstate(T, T, r, theta)
function squeezedstate(r::N, theta::N) where {N<:Real}
    mean = zeros(2)
    cr, sr = cosh(2*r), sinh(2*r)
    covar = (1/2) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    return GaussianState(mean, covar, 1)
end

"""
    eprstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] r<:Real, theta<:Real)

Gaussian state that is a two-mode squeezed state, known as the Einstein-Podolski-Rosen (EPR) state.
The amplitude and phase squeezing parameters are given by `r` and `theta`, respectively.

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
julia> eprstate(0.5, pi/4)
GaussianState for 2 modes.
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
function eprstate(::Type{Tm}, ::Type{Tc}, r::N, theta::N) where {Tm,Tc,N<:Real}
    mean = Tm(zeros(4))
    v1 = (1/2) * cosh(2*r) * Matrix{Float64}(I, 2, 2)
    v2 = (1/2) * sinh(2*r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    covar = Tc([v1 v2; v2 v1])
    return GaussianState(mean, covar, 2)
end
eprstate(::Type{T}, r::N, theta::N) where {T,N<:Real} = eprstate(T, T, r, theta)
function eprstate(r::N, theta::N) where {N<:Real}
    mean = zeros(4)
    v1 = (1/2) * cosh(2*r) * Matrix{Float64}(I, 2, 2)
    v2 = (1/2) * sinh(2*r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    covar = [v1 v2; v2 v1]
    return GaussianState(mean, covar, 2)
end

##
# Operations on Gaussian states
##

"""
    tensor(state1::GaussianState, state2::GaussianState)

tensor product of Gaussian states, which can also be called with `⊗`.

## Example
```jldoctest
julia> coherentstate(1.0+im) ⊗ thermalstate(2)
GaussianState for 2 modes.
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
    mean′, covar′ = _tensor_fields(state1, state2)
    return GaussianState(Tm(mean′), Tc(covar′), state1.nmodes + state2.nmodes)
end
tensor(::Type{T}, state1::GaussianState, state2::GaussianState) where {T} = tensor(T, T, state1, state2)
function tensor(state1::GaussianState, state2::GaussianState)
    mean′, covar′ = _tensor_fields(state1, state2)
    return GaussianState(mean′, covar′, state1.nmodes + state2.nmodes)
end
function _tensor_fields(state1::GaussianState, state2::GaussianState)
    mean1, mean2 = state1.mean, state2.mean
    length1, length2 = 2*state1.nmodes, 2*state2.nmodes
    slengths = length1 + length2
    covar1, covar2 = state1.covar, state2.covar
    # initialize direct sum of mean vectors
    mean′ = zeros(length1+length2)
    @inbounds for i in eachindex(mean1)
        mean′[i] = mean1[i]
    end
    @inbounds for i in eachindex(mean2)
        mean′[i+length1] = mean2[i]
    end
    # initialize direct sum of covariance matrices
    covar′ = zeros(slengths, slengths)
    axes1 = axes(covar1)
    @inbounds for i in axes1[1], j in axes1[2]
        covar′[i,j] = covar1[i,j]
    end
    axes2 = axes(covar2)
    @inbounds for i in axes2[1], j in axes2[2]
        covar′[i+length1,j+length1] = covar2[i,j]
    end
    # extract output array types
    mean′′ = _promote_output_vector(typeof(mean1), typeof(mean2), mean′)
    covar′′ = _promote_output_matrix(typeof(covar1), typeof(covar2), covar′)
    return mean′′, covar′′
end

"""
    ptrace([Tm=Vector{Float64}, Tc=Matrix{Float64},] state::GaussianState, idx<:Int)
    ptrace([Tm=Vector{Float64}, Tc=Matrix{Float64},] state::GaussianState, indices<:AbstractVector)

Partial trace of a Gaussian state over a subsytem indicated by `idx`, or multiple subsystems
indicated by `indices`.

## Example
```jldoctest
julia> state = coherentstate(1.0+im) ⊗ thermalstate(2) ⊗ squeezedstate(3.0, pi/4)
GaussianState for 3 modes.
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
GaussianState for 1 mode.
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
 2.5  0.0
 0.0  2.5

julia> ptrace(state, [1, 3])
GaussianState for 2 modes.
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
    mean′, covar′ = _ptrace_fields(state, idx)
    return GaussianState(Tm(mean′), Tc(covar′), 1)
end
ptrace(::Type{T}, state::GaussianState, idx::N) where {T,N<:Int} = ptrace(T, T, state, idx)
function ptrace(state::GaussianState, idx::N) where {N<:Int}
    mean′, covar′ = _ptrace_fields(state, idx)
    return GaussianState(mean′, covar′, 1)
end
function ptrace(::Type{Tm}, ::Type{Tc}, state::GaussianState, indices::N) where {Tm,Tc,N<:AbstractVector}
    mean′, covar′ = _ptrace_fields(state, indices)
    return GaussianState(Tm(mean′), Tc(covar′), length(indices))
end
ptrace(::Type{T}, state::GaussianState, indices::N) where {T,N<:AbstractVector} = ptrace(T, T, state, indices)
function ptrace(state::GaussianState, indices::T) where {T<:AbstractVector}
    mean′, covar′ = _ptrace_fields(state, indices)
    return GaussianState(mean′, covar′, length(indices))
end
function _ptrace_fields(state::GaussianState, idx::T) where {T<:Int}
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
function _ptrace_fields(state::GaussianState, indices::T) where {T<:AbstractVector}
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