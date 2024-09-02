function displace(::Type{Td}, ::Type{Tt}, alpha::N) where {Td,Tt,N<:Number}
    disp = sqrt(2) * Td([real(alpha), imag(alpha)])
    transform = Tt(Matrix{Float64}(I, 2, 2))
    noise = Tt(zeros(2, 2))
    return GaussianChannel(disp, transform, noise)
end
displace(::Type{T}, alpha::N) where {T,N<:Number} = displace(T, T, alpha)
function displace(alpha::N) where {N<:Number}
    disp = sqrt(2) * [real(alpha), imag(alpha)]
    transform = Matrix{Float64}(I, 2, 2)
    noise = zeros(2, 2)
    return GaussianChannel(disp, transform, noise)
end

function squeeze(::Type{Td}, ::Type{Tt}, r::N, theta::N) where {Td,Tt,N<:Real}
    disp = Td(zeros(2))
    cr, sr = cosh(r), sinh(r)
    t = sinh(r) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    transform = Tt(t)
    noise = Tt(zeros(2, 2))
    return GaussianChannel(disp, transform, noise)
end
squeeze(::Type{T}, r::N, theta::N) where {T,N<:Real} = squeeze(T, T, r, theta)
function squeeze(r::N, theta::N) where {N<:Real}
    disp = zeros(2)
    cr, sr = cosh(r), sinh(r)
    transform = sinh(r) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    noise = zeros(2, 2)
    return GaussianChannel(disp, transform, noise)
end

function twosqueeze(::Type{Td}, ::Type{Tt}, r::N, theta::N) where {Td,Tt,N<:Real}
    disp = Td(zeros(4))
    v1 = cosh(r) * Matrix{Float64}(I, 2, 2)
    v2 = sinh(r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    transform = Tt([v1 -v2; -v2 v1])
    noise = Tt(zeros(4, 4))
    return GaussianChannel(disp, transform, noise)
end
twosqueeze(::Type{T}, r::N, theta::N) where {T,N<:Real} = twosqueeze(T, T, r, theta)
function twosqueeze(r::N, theta::N) where {N<:Real}
    disp = zeros(4)
    v1 = cosh(r) * Matrix{Float64}(I, 2, 2)
    v2 = sinh(r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    transform = [v1 -v2; -v2 v1]
    noise = zeros(4, 4)
    return GaussianChannel(disp, transform, noise)
end

function phaseshift(::Type{Td}, ::Type{Tt}, theta::N) where {Td,Tt,N<:Real}
    disp = Td(zeros(2))
    transform = Tt([cos(theta) sin(theta); -sin(theta) cos(theta)])
    noise = Tt(zeros(2, 2))
    return GaussianChannel(disp, transform, noise)
end
phaseshift(::Type{T}, theta::N) where {T,N<:Real} = phaseshift(T, T, theta)
function phaseshift(theta::N) where {N<:Real}
    disp = zeros(2)
    transform = [cos(theta) sin(theta); -sin(theta) cos(theta)]
    noise = zeros(2, 2)
    return GaussianChannel(disp, transform, noise)
end

function beamsplitter(::Type{Td}, ::Type{Tt}, transmit::N) where {Td,Tt,N<:Real}
    disp = Td(zeros(4))
    I2 = Matrix{Float64}(I, 2, 2)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    transform = Tt([a1*I2 a2*I2; -a2*I2 a1*I2])
    noise = Tt(zeros(4, 4))
    return GaussianChannel(disp, transform, noise)
end
beamsplitter(::Type{T}, transmit::N) where {T,N<:Real} = beamsplitter(T, T, transmit)
function beamsplitter(transmit::N) where {N<:Real}
    disp = zeros(4)
    I2 = Matrix{Float64}(I, 2, 2)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    transform = [a1*I2 a2*I2; -a2*I2 a1*I2]
    noise = zeros(4, 4)
    return GaussianChannel(disp, transform, noise)
end