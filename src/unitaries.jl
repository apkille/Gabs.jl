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