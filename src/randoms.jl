"""
    randstate([Tm=Vector{Float64}, Tc=Matrix{Float64},] nmodes<:Int)

Calculate a random Gaussian state.
"""
function randstate(::Type{Tm}, ::Type{Tc}, nmodes::N) where {Tm,Tc,N<:Int}
    mean = randn(2*nmodes)
    covar = randn(2*nmodes, 2*nmodes)
    return GaussianState(Tm(mean), Tc(covar), nmodes)
end
randstate(::Type{T}, nmodes::N) where {T,N<:Int} = randstate(T,T,nmodes)
function randstate(nmodes::N) where {N<:Int}
    mean = randn(2*nmodes)
    covar = randn(2*nmodes, 2*nmodes)
    return GaussianState(mean, covar, nmodes)
end

"""
    randchannel([Td=Vector{Float64}, Tt=Matrix{Float64},] nmodes<:Int)

Calculate a random Gaussian channel.
"""
function randchannel(::Type{Td}, ::Type{Tt}, nmodes::N) where {Td,Tt,N<:Int}
    disp = randn(2*nmodes)
    transform = randn(2*nmodes, 2*nmodes)
    noise = randn(2*nmodes, 2*nmodes)
    return GaussianChannel(Td(disp), Tt(transform), Tt(noise), nmodes)
end
randchannel(::Type{T}, nmodes::N) where {T,N<:Int} = randchannel(T,T,nmodes)
function randchannel(nmodes::N) where {N<:Int}
    disp = randn(2*nmodes)
    transform = randn(2*nmodes, 2*nmodes)
    noise = randn(2*nmodes, 2*nmodes)
    return GaussianChannel(disp, transform, noise, nmodes)
end