"""
    purity(state::GaussianState)

Calculate the purity of a Gaussian state, defined by `1/sqrt((2/ħ) det(V))`.
"""
purity(x::GaussianState) = (b = x.basis; (x.ħ/2)^(b.nmodes)/sqrt(det(x.covar)))

"""
    entropy_vn(state::GaussianState, tol=1e-15)

Calculate the Von Neumann entropy of a Gaussian state, defined as

```math
S(\\rho) = -Tr(\\rho \\log(\\rho)) = \\sum_i f(v_i)
```

such that ``\\log`` denotes the natural logarithm, ``v_i`` is the symplectic
spectrum of ``\\mathbf{V}/\\hbar``, and the ``f`` is taken to be

```math
f(x) = (x + 1/2) \\log(x + 1/2) - (x - 1/2) \\log(x - 1/2)
```
wherein it is understood that ``0 \\log(0) \\equiv 0``.

# Arguments
* `state`: Gaussian state whose Von Neumann entropy is to be calculated.
* `tol=1e-15`: Tolerance above the logarithmic singularity.
"""
function entropy_vn(state::GaussianState, tol=1e-15)
    spectrum = sympspectrum(state) ./ state.ħ
    return reduce(+, _entropy_vn.(filter(x -> x-(1/2) > tol, spectrum)))
end

_entropy_vn(x) = (x+(1/2)) * log(x+(1/2)) - (x-(1/2)) * log(x-(1/2))

"""
    fidelity(state1::GaussianState, state2::GaussianState, tol=1e-15)

Calculate the joint fidelity of two Gaussian states, defined as

```math
F(\\rho, \\sigma) = Tr(\\sqrt{\\sqrt{\\rho} \\sigma \\sqrt{\\rho}}).
```

See: Banchi, Braunstein, and Pirandola, Phys. Rev. Lett. 115, 260501 (2015)

# Arguments
* `state1`, `state2`: Gaussian states whose joint fidelity is to be calculated.
* `tol=1e-15`: Tolerance above the square root singularity.
"""
function fidelity(state1::GaussianState, state2::GaussianState, tol=1e-15)
    state1.basis == state2.basis || throw(ArgumentError(SYMPLECTIC_ERROR))
    state1.ħ == state2.ħ || throw(ArgumentError(HBAR_ERROR))
    delta = state2.mean - state1.mean
    V_sum = state1.covar + state2.covar
    Omega = symplecticform(state1.basis)
    # slightly different from Banachi, Braunstein, and Pirandola
    V_aux = (V_sum \ ((Omega ./ 4) + (state2.covar * Omega * state1.covar)))
    spectrum = filter(x -> x > 0, imag.(eigvals(V_aux))) .* (2/state.ħ)
    F_tot = sqrt(reduce(*, _fidelity.(filter(x -> x-1 > tol, spectrum))))
    prefactor = F_tot / (det(V_sum))^(1/4)
    return prefactor * exp(- (transpose(delta) * (V_sum \ delta)) / 4)
end

_fidelity(x) = x + sqrt(x^2 - 1)

"""
    logarithmic_negativity(state::GaussianState, indices=1, tol=1e-15)

Calculate the logarithmic negativity of a Gaussian state partition, defined as

```math
N(\\rho) = \\log\\|\\rho^{T_B}\\|_1 = - \\sum_i \\log(\\tilde{v}_<^i)
```

such that ``\\log`` denotes the natural logarithm, ``\\tilde{v}_<^i`` is the
symplectic spectrum of ``\\mathbf{\\tilde{V}}/\\hbar`` which is ``< 1/2``.

Therein, ``\\mathbf{\\tilde{V}} = \\mathbf{T} \\mathbf{V} \\mathbf{T}`` where
```math
\\forall k : \\mathbf{T} q_k = q_k
\\forall k \\in \\mathrm{B} : \\mathbf{T} p_k = -p_k
\\forall k \\notin \\mathrm{B} : \\mathbf{T} p_k = p_k
```

# Arguments
* `state`: Gaussian state whose logarithmic negativity is to be calculated.
* `indices`: integer or collection thereof, specifying the binary partition.
* `tol=1e15`: Tolerance above the logarithmic singularity.
"""
function logarithmic_negativity(state::GaussianState, indices=1, tol=1e-15)
    tilde = _tilde(state)
    M = symplecticform(state.basis) * tilde
    spectrum = filter(x -> x > 0, imag.(eigvals(M))) ./ state.ħ
    return -reduce(+, log.(filter(x -> x > tol && x < 1/2, spectrum)))
end

function _tilde(state::GaussianState{B,M,V}, indices=1) where {B<:QuadPairBasis,M,V}
    T = state.covar
    nmodes = state.basis.nmodes
    for i in indices
        # first loop is cache friendly, second one thrashes
        @inbounds for j in Base.OneTo(2*nmodes)
	    T[j, 2*i] *= -1
        end
        @inbounds for j in Base.OneTo(2*nmodes)
            T[2*i, j] *= -1
        end
    end
    return T
end
function _tilde(state::GaussianState{B,M,V}, indices=1) where {B<:QuadBlockBasis,M,V}
    T = state.covar
    nmodes = state.basis.nmodes
    for i in indices
        # first loop is cache friendly, second one thrashes
        @inbounds for j in Base.OneTo(2*nmodes)
	    T[j, nmodes + i] *= -1
        end
        @inbounds for j in Base.OneTo(2*nmodes)
            T[nmodes + i, j] *= -1
        end
    end
    return T
end
