"""
    wigner(state::GaussianState, x)

Compute the Wigner function of an N-mode Gaussian state at `x`, a vector of size 2N.
"""
function wigner(state::GaussianState, x::T) where {T}
    basis = state.basis
    nmodes = basis.nmodes
    mean = state.mean
    isequal(length(mean), length(x)) || throw(ArgumentError(WIGNER_ERROR))

    V = state.covar
    diff = x .- mean
    arg = -(1/2) * transpose(diff) * inv(V) * diff

    return exp(arg)/((2pi)^nmodes * sqrt(det(V)))
end

"""
    wignerchar(state::GaussianState, xi)

Compute the Wigner characteristic function of an N-mode Gaussian state at `xi`, a vector of size 2N.
"""
function wignerchar(state::GaussianState, xi::T) where {T}
    basis = state.basis
    nmodes = basis.nmodes
    mean = state.mean
    isequal(length(mean), length(xi)) || throw(ArgumentError(WIGNER_ERROR))

    V = state.covar
    Omega = symplecticform(basis)

    arg1 = -(1/2) * transpose(xi) * (Omega*V*transpose(Omega))*xi
    arg2 = im * transpose(Omega*mean) * xi

    return exp(arg1 .- arg2)
end