"""
    randstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] nmodes<:Int; pure=false)

Calculate a random Gaussian state.
"""
function randstate(::Type{Tm}, ::Type{Tc}, repr::SymplecticRepr{N}; pure = false) where {Tm,Tc,N<:Int}
    mean, covar = _randstate(repr, pure)
    return GaussianState(repr, Tm(mean), Tc(covar))
end
randstate(::Type{T}, repr::SymplecticRepr{N}; pure = false) where {T,N<:Int} = randstate(T,T,repr,pure=pure)
function randstate(repr::SymplecticRepr{N}; pure = false) where {N<:Int}
    mean, covar = _randstate(repr, pure)
    return GaussianState(repr, mean, covar)
end
function _randstate(repr::CanonicalForm{N}, pure) where {N<:Int}
    nmodes = repr.nmodes
    mean = randn(2*nmodes)
    covar = zeros(2*nmodes, 2*nmodes)
    symp = randsymplectic(repr)
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
function randunitary(::Type{Td}, ::Type{Ts}, repr::SymplecticRepr{N}; passive = false) where {Td,Ts,N<:Int}
    disp, symp = _randunitary(repr, passive)
    return GaussianUnitary(repr, Td(disp), Ts(symp))
end
randunitary(::Type{T}, repr::SymplecticRepr{N}; passive = false) where {T,N<:Int} = randunitary(T,T,repr; passive = passive)
function randunitary(repr::SymplecticRepr{N}; passive = false) where {N<:Int}
    disp, symp = _randunitary(repr, passive)
    return GaussianUnitary(repr, disp, symp)
end
function _randunitary(repr::SymplecticRepr{N}, passive) where {N<:Int}
    nmodes = repr.nmodes
    disp = rand(2*nmodes)
    symp = randsymplectic(repr, passive = passive)
    return disp, symp
end

"""
    randchannel([Td=Vector{Float64}, Tt=Matrix{Float64},] nmodes<:Int)

Calculate a random Gaussian channel.
"""
function randchannel(::Type{Td}, ::Type{Tt}, repr::SymplecticRepr{N}) where {Td,Tt,N<:Int}
    disp, transform, noise = _randchannel(repr)
    return GaussianChannel(repr, Td(disp), Tt(transform), Tt(noise))
end
randchannel(::Type{T}, repr::SymplecticRepr{N}) where {T,N<:Int} = randchannel(T,T,repr)
function randchannel(repr::SymplecticRepr{N}) where {N<:Int}
    disp, transform, noise = _randchannel(repr)
    return GaussianChannel(repr, disp, transform, noise)
end
function _randchannel(repr::SymplecticRepr{N}) where {N<:Int}
    nmodes = repr.nmodes
    disp = randn(2*nmodes)
    # generate symplectic matrix describing the evolution of the system with N modes
    # and environment with 2N modes
    symp = randsymplectic(typeof(repr)(3*nmodes))
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
function randsymplectic(repr::CanonicalForm{N}; passive = false) where {N<:Int}
    nmodes = repr.nmodes
    # generate random orthogonal symplectic matrix
    O = _rand_orthogonal_symplectic(repr)
    if passive
        return O
    end
    # instead generate random active symplectic matrix
    O′ = _rand_orthogonal_symplectic(repr)
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
function randsymplectic(::Type{T}, repr::CanonicalForm{N}; passive = false) where {T, N<:Int} 
    symp = randsymplectic(repr, passive = passive)
    return T(symp)
end

# Generates random orthogonal symplectic matrix by blocking real
# and imaginary parts of a random unitary matrix
function _rand_orthogonal_symplectic(repr::CanonicalForm{N}) where {N<:Int}
    nmodes = repr.nmodes
    U = _rand_unitary(repr)
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
# This approach is faster and creates less allocations than rand(Haar(2), nmodes) from RandomMatrices.jl
function _rand_unitary(repr::CanonicalForm{N}) where {N<:Int}
    nmodes = repr.nmodes
    M = rand(ComplexF64, nmodes, nmodes) ./ sqrt(2.0)
    q, r = qr(M)
    d = diagm([r[i, i] / abs(r[i, i]) for i in Base.OneTo(nmodes)])
    return q * d
end