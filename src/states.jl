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