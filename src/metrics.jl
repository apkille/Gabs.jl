"""
    purity(state::GaussianState)

Calculate the purity of a Gaussian state, defined by `1/sqrt((2/ħ) det(V))`.
"""
purity(x::GaussianState) = (b = x.basis; (x.ħ/2)^(b.nmodes)/sqrt(det(x.covar)))

"""
    entropy_vn(state::GaussianState; refine = x -> true)

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
* `refine: (Optional) Function(::Number)::Bool used for narrowing down the
	   range of input values passed to f(x).

## Example

```julia
entropy_vn_smoothed(state, tol) = entropy_vn(state; refine = x -> (x - (1/2)) > tol)
```
"""
function entropy_vn(state::GaussianState; refine = x -> true)
    T = symplecticform(state.basis) * state.covar
    T = filter(x -> x > 1/2 && refine(x), imag.(eigvals(T)) ./ state.ħ)
    return reduce(+, _entropy_vn.(T))
end

# this is the same as f(x)
_entropy_vn(x) = x < 19 ?
    (x + (1/2)) * log(x + (1/2)) - (x - (1/2)) * log(x - (1/2)) :
    log(x) + 1 - (1/(24 * x^2)) - (1/(320 * x^4)) - (1/(2688 * x^6))

"""
    fidelity(state1::GaussianState, state2::GaussianState; refine = x -> true)

Calculate the joint fidelity of two Gaussian states, defined as

```math
F(\\rho, \\sigma) = Tr(\\sqrt{\\sqrt{\\rho} \\sigma \\sqrt{\\rho}}).
```

See: Banchi, Braunstein, and Pirandola, Phys. Rev. Lett. 115, 260501 (2015)

# Arguments
* `state1`, `state2`: Gaussian states whose joint fidelity is to be calculated.
* `refine: (Optional) Function(::Number)::Bool used for narrowing down the
	   range of input values passed to `x -> x + sqrt(x^2 - 1)`.

## Example
```julia
fidelity_smoothed(state1, state2, tol) = fidelity(state1, state2; refine = x -> (x - 1) > tol)
```
"""
function fidelity(state1::GaussianState, state2::GaussianState; refine = x -> true)
    state1.basis == state2.basis || throw(ArgumentError(SYMPLECTIC_ERROR))
    state1.ħ == state2.ħ || throw(ArgumentError(HBAR_ERROR))
    A = state2.mean - state1.mean
    B = state1.covar + state2.covar
    # many nasty factors of ħ ahead, tread carefully
    output = state1.ħ^(state1.basis.nmodes/2) * exp(- (transpose(A) * (B \ A)) / 4) / (det(B))^(1/4)
    A = symplecticform(state1.basis)
    # slightly different from Banachi, Braunstein, and Pirandola
    B = (B \ ((A .* ((state1.ħ^2)/4)) + (state2.covar * A * state1.covar)))
    B = filter(x -> x > 1 && refine(x), imag.(eigvals(B)) .* (2/state1.ħ))
    return output * sqrt(reduce(*, _fidelity.(B)))
end

# this is the same as x + sqrt(x^2 - 1) when x > 0, but overflows gradually
_fidelity(x) = x^2 < floatmax(typeof(x)) ? x + sqrt(x^2 - 1) : 2 * x

"""
    logarithmic_negativity(state::GaussianState, indices; refine = x -> true)

Calculate the logarithmic negativity of a Gaussian state partition, defined as

```math
N(\\rho) = \\log\\|\\rho^{T_B}\\|_1 = - \\sum_i \\log(2 \\tilde{v}_i^<)
```

such that ``\\log`` denotes the natural logarithm, ``\\tilde{v}_i^<`` is the
symplectic spectrum of ``\\mathbf{\\tilde{V}}/\\hbar`` which is ``< 1/2``.

Therein, ``\\mathbf{\\tilde{V}} = \\mathbf{T} \\mathbf{V} \\mathbf{T}`` where
```math
\\forall k : \\mathbf{T} q_k = q_k
\\forall k \\in \\mathrm{B} : \\mathbf{T} p_k = -p_k
\\forall k \\notin \\mathrm{B} : \\mathbf{T} p_k = p_k
```

# Arguments
* `state`: Gaussian state whose logarithmic negativity is to be calculated.
* `indices`: Integer or collection thereof, specifying the binary partition.
* `refine: (Optional) Function(::Number)::Bool used for narrowing down the
	   range of input values passed to log(x).

## Example

```julia
function logarithmic_negativity_smoothed(state, indices, tol1, tol2)
    logarithmic_negativity(state, indices; refine = x -> (x > tol1 && (1 - x) > tol2))
end
```
"""
function logarithmic_negativity(state::GaussianState, indices; refine = x -> true)
    T = _tilde(state, indices)
    T = symplecticform(state.basis) * T
    T = filter(x -> x < 1 && refine(x), imag.(eigvals(T)) .* (2/state.ħ))
    T = reduce(+, log.(T))
    # in case the reduction happened over an empty set
    return T < 0 ? -T : T
end

function _tilde(state::GaussianState{B,M,V}, indices) where {B<:QuadPairBasis,M,V}
    nmodes = state.basis.nmodes
    indices = collect(indices)
    all(x -> x >= 1 && x <= nmodes, indices) || throw(ArgumentError(INDEX_ERROR))
    T = copy(state.covar)
    @inbounds for i in indices
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
function _tilde(state::GaussianState{B,M,V}, indices) where {B<:QuadBlockBasis,M,V}
    nmodes = state.basis.nmodes
    indices = collect(indices)
    all(x -> x >= 1 && x <= nmodes, indices) || throw(ArgumentError(INDEX_ERROR))
    T = copy(state.covar)
    @inbounds for i in indices
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
