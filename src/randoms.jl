"""
    randstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] nmodes<:Int; pure=false)

Calculate a random Gaussian state.
"""
function randstate(::Type{Tm}, ::Type{Tc}, nmodes::N; pure = false) where {Tm,Tc,N<:Int}
    mean, covar = _randstate_fields(nmodes, pure)
    return GaussianState(Tm(mean), Tc(covar), nmodes)
end
randstate(::Type{T}, nmodes::N; pure = false) where {T,N<:Int} = randstate(T,T,nmodes,pure=pure)
function randstate(nmodes::N; pure = false) where {N<:Int}
    mean, covar = _randstate_fields(nmodes, pure)
    return GaussianState(mean, covar, nmodes)
end
function _randstate_fields(nmodes::N, pure) where {N<:Int}
    mean = randn(2*nmodes)
    covar = zeros(2*nmodes, 2*nmodes)
    symp = randsymplectic(nmodes)
    # generate pure Gaussian state
    if pure
        mul!(covar, symp, symp')
        return mean, covar
    end
    # create buffer for matrix multiplication
    buf = zeros(2*nmodes, 2*nmodes)
    # William decomposition for mixed Gaussian states
    sympeigs = abs.(rand(nmodes)) .+ 1.0
    will = diagm(repeat(sympeigs, inner = 2))
    mul!(covar, symp, mul!(buf, will, symp'))
    return mean, covar
end

"""
    randunitary([Td=Vector{Float64}, Ts=Matrix{Float64},] nmodes<:Int; passive=false)

Calculate a random Gaussian unitary operator.
"""
function randunitary(::Type{Td}, ::Type{Ts}, nmodes::N; passive = false) where {Td,Ts,N}
    disp, symp = _randunitary_fields(nmodes, passive)
    return GaussianUnitary(Td(disp), Ts(symp))
end
randunitary(::Type{T}, nmodes::N; passive = false) where {T,N<:Int} = randunitary(T,T,nmodes; passive = passive)
function randunitary(nmodes::N; passive = false) where {N<:Int}
    disp, symp = _randunitary_fields(nmodes, passive)
    return GaussianUnitary(disp, symp)
end
function _randunitary_fields(nmodes::N, passive) where {N<:Int}
    disp = rand(2*nmodes)
    symp = randsymplectic(nmodes, passive = passive)
    return disp, symp
end

"""
    randchannel([Td=Vector{Float64}, Tt=Matrix{Float64},] nmodes<:Int)

Calculate a random Gaussian channel.
"""
function randchannel(::Type{Td}, ::Type{Tt}, nmodes::N) where {Td,Tt,N<:Int}
    disp, transform, noise = _randchannel(nmodes)
    return GaussianChannel(Td(disp), Tt(transform), Tt(noise), nmodes)
end
randchannel(::Type{T}, nmodes::N) where {T,N<:Int} = randchannel(T,T,nmodes)
function randchannel(nmodes::N) where {N<:Int}
    disp, transform, noise = _randchannel(nmodes)
    return GaussianChannel(disp, transform, noise, nmodes)
end
function _randchannel(nmodes::N) where {N<:Int}
    disp = randn(2*nmodes)
    # generate symplectic matrix describing the evolution of the system with N modes
    # and environment with 2N modes
    symp = randsymplectic(3*nmodes)
    transform, B = symp[1:2*nmodes, 1:2*nmodes], @view(symp[1:2*nmodes, 2*nmodes+1:6*nmodes])
    # generate noise matrix from evolution of environment
    noise = zeros(2*nmodes, 2*nmodes)
    mul!(noise, B, B')
    return disp, transform, noise
end

"""
    randsymplectic([T=Matrix{Float64},] nmodes<:Int, passive=false)

Calculate a random symplectic matrix of size `2*nmodes x 2*nmodes`.
"""
function randsymplectic(nmodes::N; passive = false) where {N<:Int}
    # generate random orthogonal symplectic matrix
    O = _rand_orthogonal_symplectic(nmodes)
    if passive
        return O
    end
    # instead generate random active symplectic matrix
    O′ = _rand_orthogonal_symplectic(nmodes)
    # direct sum of symplectic matrices for single-mode squeeze transformations
    rs = rand(nmodes)
    squeezes = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes)
        val = rs[i]
        squeezes[2*i-1, 2*i-1] = val
        squeezes[2*i, 2*i] = 1/val
    end
    return O * squeezes * O′
end
function randsymplectic(::Type{T}, nmodes::N; passive = false) where {T, N<:Int} 
    symp = randsymplectic(nmodes, passive = passive)
    return T(symp)
end

# Generates random orthogonal symplectic matrix by blocking real
# and imaginary parts of a random unitary matrix
function _rand_orthogonal_symplectic(nmodes::N) where {N<:Int}
    U = _rand_unitary(nmodes)
    O = zeros(2*nmodes, 2*nmodes)
    @inbounds for i in Base.OneTo(nmodes), j in Base.OneTo(nmodes)
        val = U[i,j]
        O[2*i-1,2*j-1] = real(val)
        O[2*i, 2*j-1] = -imag(val)
        O[2*i-1, 2*j] = imag(val)
        O[2*i, 2*j] = real(val)
    end
    return O
end
# Generates unitary matrix randomly distributed over Haar measure;
# see https://arxiv.org/abs/math-ph/0609050 for algorithm.
# This approach is faster and creates less allocations than rand(Haar(2), size) from RandomMatrices.jl
function _rand_unitary(size::N) where {N<:Int}
    M = rand(ComplexF64, size, size) ./ sqrt(2.0)
    q, r = qr(M)
    d = diagm([r[i, i] / abs(r[i, i]) for i in Base.OneTo(size)])
    return q * d
end