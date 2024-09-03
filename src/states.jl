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
GaussianState
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
    return GaussianState(mean, covar)
end
vacuumstate(::Type{T}) where {T} = vacuumstate(T, T)
function vacuumstate()
    mean = zeros(2)
    covar = Matrix{Float64}(I, 2, 2)
    return GaussianState(mean, covar)
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
GaussianState
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
    return GaussianState(mean, covar)
end
thermalstate(::Type{T}, photons::N) where {T, N<:Int} = thermalstate(T, T, photons)
function thermalstate(photons::N) where {N<:Int}
    mean = zeros(2)
    covar = (photons + 1/2) * Matrix{Float64}(I, 2, 2)
    return GaussianState(mean, covar)
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
\\quad \\mathbf{V} = \\left(\\bar{n} + \\frac{1}{2}\\right)\\mathbf{I}.
```

## Example

```jldoctest
julia> coherentstate(1.0+im)
GaussianState
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
    return GaussianState(mean, covar)
end
coherentstate(::Type{T}, alpha::N) where {T, N<:Number} = coherentstate(T, T, alpha)
function coherentstate(alpha::N) where {N<:Number}
    mean = sqrt(2) * [real(alpha), imag(alpha)]
    covar = Matrix{Float64}(I, 2, 2)
    return GaussianState(mean, covar)
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
GaussianState
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
    return GaussianState(mean, covar)
end
squeezedstate(::Type{T}, r::N, theta::N) where {T,N<:Real} = squeezedstate(T, T, r, theta)
function squeezedstate(r::N, theta::N) where {N<:Real}
    mean = zeros(2)
    cr, sr = cosh(2*r), sinh(2*r)
    covar = (1/2) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    return GaussianState(mean, covar)
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
GaussianState
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
    return GaussianState(mean, covar)
end
eprstate(::Type{T}, r::N, theta::N) where {T,N<:Real} = eprstate(T, T, r, theta)
function eprstate(r::N, theta::N) where {N<:Real}
    mean = zeros(4)
    v1 = (1/2) * cosh(2*r) * Matrix{Float64}(I, 2, 2)
    v2 = (1/2) * sinh(2*r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    covar = [v1 v2; v2 v1]
    return GaussianState(mean, covar)
end