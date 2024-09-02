function displace(::Type{Td}, ::Type{Tt}, alpha::A, noise::N) where {Td,Tt,A<:Number,N}
    disp = sqrt(2) * Td([real(alpha), imag(alpha)])
    transform = Tt(Matrix{Float64}(I, 2, 2))
    return GaussianChannel(disp, transform, Tt(noise))
end
displace(::Type{T}, alpha::A, noise::N) where {T,A<:Number,N} = displace(T, T, alpha, noise)
function displace(alpha::A, noise::N) where {A<:Number,N}
    disp = sqrt(2) * [real(alpha), imag(alpha)]
    transform = Matrix{Float64}(I, 2, 2)
    return GaussianChannel(disp, transform, noise)
end

function squeeze(::Type{Td}, ::Type{Tt}, r::R, theta::R, noise::N) where {Td,Tt,R<:Real,N}
    disp = Td(zeros(2))
    cr, sr = cosh(r), sinh(r)
    t = sinh(r) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    transform = Tt(t)
    return GaussianChannel(disp, transform, Tt(noise))
end
squeeze(::Type{T}, r::R, theta::R, noise::N) where {T,R<:Real,N} = squeeze(T, T, r, theta, noise)
function squeeze(r::R, theta::R, noise::N) where {R<:Real,N}
    disp = zeros(2)
    cr, sr = cosh(r), sinh(r)
    transform = sinh(r) * [cr-sr*cos(theta) sr*sin(theta); sr*sin(theta) cr+sr*cos(theta)]
    return GaussianChannel(disp, transform, noise)
end

function twosqueeze(::Type{Td}, ::Type{Tt}, r::R, theta::R, noise::N) where {Td,Tt,R<:Real,N}
    disp = Td(zeros(4))
    v1 = cosh(r) * Matrix{Float64}(I, 2, 2)
    v2 = sinh(r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    transform = Tt([v1 -v2; -v2 v1])
    return GaussianChannel(disp, transform, Tt(noise))
end
twosqueeze(::Type{T}, r::R, theta::R, noise::N) where {T,R<:Real,N} = twosqueeze(T, T, r, theta, noise)
function twosqueeze(r::R, theta::R, noise::N) where {R<:Real,N}
    disp = zeros(4)
    v1 = cosh(r) * Matrix{Float64}(I, 2, 2)
    v2 = sinh(r) * [cos(theta) sin(theta); sin(theta) -cos(theta)]
    transform = [v1 -v2; -v2 v1]
    return GaussianChannel(disp, transform, noise)
end

function phaseshift(::Type{Td}, ::Type{Tt}, theta::R, noise::N) where {Td,Tt,R<:Real,N}
    disp = Td(zeros(2))
    transform = Tt([cos(theta) sin(theta); -sin(theta) cos(theta)])
    return GaussianChannel(disp, transform, Tt(noise))
end
phaseshift(::Type{T}, theta::R, noise::N) where {T,R<:Real,N} = phaseshift(T, T, theta, noise)
function phaseshift(theta::R, noise::N) where {R<:Real,N}
    disp = zeros(2)
    transform = [cos(theta) sin(theta); -sin(theta) cos(theta)]
    return GaussianChannel(disp, transform, noise)
end

function beamsplitter(::Type{Td}, ::Type{Tt}, transmit::R, noise::N) where {Td,Tt,R<:Real,N}
    disp = Td(zeros(4))
    I2 = Matrix{Float64}(I, 2, 2)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    transform = Tt([a1*I2 a2*I2; -a2*I2 a1*I2])
    return GaussianChannel(disp, transform, Tt(noise))
end
beamsplitter(::Type{T}, transmit::R, noise::N) where {T,R<:Real,N} = beamsplitter(T, T, transmit, noise)
function beamsplitter(transmit::R, noise::N) where {R<:Real,N}
    disp = zeros(4)
    I2 = Matrix{Float64}(I, 2, 2)
    a1, a2 = sqrt(transmit), sqrt(1 - transmit)
    transform = [a1*I2 a2*I2; -a2*I2 a1*I2]
    return GaussianChannel(disp, transform, noise)
end