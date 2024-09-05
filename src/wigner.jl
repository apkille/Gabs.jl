"""
    symplecticform([T = Matrix{Float64},] modes<:Int)

Compute the symplectic form matrix of size 2N x 2N, where N is given by `modes`.
"""
function symplecticform(modes::N) where {N<:Int}
    Omega = zeros(2modes, 2modes)
    @inbounds for i in 1:modes, j in modes:2modes
        if isequal(i, j-modes)
            Omega[i,j] = 1.0
        end
    end
    @inbounds for i in modes:2modes, j in 1:modes
        if isequal(i-modes,j)
            Omega[i, j] = -1.0
        end
    end
    return Omega
end
symplecticform(::Type{T}, modes::N) where {T, N<:Int} = T(symplecticform(modes))

"""
    wigner(state::GaussianState, x)

Compute the Wigner function of an N-mode Gaussian state at `x`, a vector of size 2N.
"""
function wigner(state::GaussianState, x::T) where {T}
    mean = state.mean
    isequal(length(mean), length(x)) || throw(ArgumentError(WIGNER_ERROR))

    V = state.covar
    modes = length(mean)/2
    diff = x .- mean
    arg = -(1/2) * transpose(diff) * inv(V) * diff

    return exp(arg)/((2pi)^modes * sqrt(det(V)))
end

"""
    wignerchar(state::GaussianState, xi)

Compute the Wigner characteristic function of an N-mode Gaussian state at `xi`, a vector of size 2N.
"""
function wignerchar(state::GaussianState, xi::T) where {T}
    mean = state.mean
    isequal(length(mean), length(xi)) || throw(ArgumentError(WIGNER_ERROR))

    V = state.covar
    modes = Int(length(mean)/2)
    Omega = symplecticform(modes)

    arg1 = -(1/2) * transpose(xi) * (Omega*V*transpose(Omega))*xi
    arg2 = im * transpose(Omega*mean) * xi

    return exp(arg1 .- arg2)
end