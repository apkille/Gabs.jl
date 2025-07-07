# src/nongaussian_states.jl
"""
    catstate_even(basis::SymplecticBasis, α::Number; squeeze_params=nothing, ħ=2)

Create an even Schrödinger cat state `|cat+⟩ = (|α⟩ + |-α⟩)/√N` where `N` is the normalization constant.

The even cat state is a quantum superposition of two coherent states (or squeezed states) with 
opposite phases. The normalization accounts for the overlap between the two component states.

# Arguments
- `basis::SymplecticBasis`: The symplectic basis for the state
- `α::Number`: Complex amplitude of the coherent states
- `squeeze_params=nothing`: Optional tuple `(r, θ)` for squeezed cat states where `r` is the squeeze parameter and `θ` is the squeeze angle
- `ħ=2`: Reduced Planck constant

# Mathematical Description:
For coherent states: `|cat+⟩ = (|α⟩ + |-α⟩)/√(2(1 + exp(-2|α|²)))`
For squeezed states: The squeezed vacuum is first created, then displaced to `±α`.
"""
function catstate_even(basis::SymplecticBasis, α::Number; squeeze_params=nothing, ħ=2)
    if squeeze_params === nothing
        state_plus = coherentstate(basis, α, ħ=ħ)
        state_minus = coherentstate(basis, -α, ħ=ħ)
    else
        r, θ = squeeze_params
        vac = vacuumstate(basis, ħ=ħ)
        squeeze_op = squeeze(basis, r, θ, ħ=ħ)
        squeezed_vac = squeeze_op * vac
        displace_plus = displace(basis, α, ħ=ħ)
        displace_minus = displace(basis, -α, ħ=ħ)
        state_plus = displace_plus * squeezed_vac
        state_minus = displace_minus * squeezed_vac
    end
    overlap = exp(-2 * abs2(α))
    norm_factor = 1 / sqrt(2 * (1 + overlap))
    coeffs = [norm_factor, norm_factor]
    states = [state_plus, state_minus]
    return GaussianLinearCombination(basis, coeffs, states)
end

"""
    catstate_odd(basis::SymplecticBasis, α::Number; squeeze_params=nothing, ħ=2)

Create an odd Schrödinger cat state `|cat-⟩ = (|α⟩ - |-α⟩)/√N` where `N` is the normalization constant.

The odd cat state is a quantum superposition of two coherent states (or squeezed states) with 
opposite phases and a relative minus sign.

# Arguments
- `basis::SymplecticBasis`: The symplectic basis for the state
- `α::Number`: Complex amplitude of the coherent states
- `squeeze_params=nothing`: Optional tuple `(r, θ)` for squeezed cat states
- `ħ=2`: Reduced Planck constant

# Mathematical Description
For coherent states: `|cat-⟩ = (|α⟩ - |-α⟩)/√(2(1 - exp(-2|α|²)))`

"""
function catstate_odd(basis::SymplecticBasis, α::Number; squeeze_params=nothing, ħ=2)
    if squeeze_params === nothing
        state_plus = coherentstate(basis, α, ħ=ħ)
        state_minus = coherentstate(basis, -α, ħ=ħ)
    else
        r, θ = squeeze_params
        vac = vacuumstate(basis, ħ=ħ)
        squeeze_op = squeeze(basis, r, θ, ħ=ħ)
        squeezed_vac = squeeze_op * vac
        displace_plus = displace(basis, α, ħ=ħ)
        displace_minus = displace(basis, -α, ħ=ħ)
        state_plus = displace_plus * squeezed_vac
        state_minus = displace_minus * squeezed_vac
    end
    overlap = exp(-2 * abs2(α))
    norm_factor = 1 / sqrt(2 * (1 - overlap))
    coeffs = [norm_factor, -norm_factor]
    states = [state_plus, state_minus]
    return GaussianLinearCombination(basis, coeffs, states)
end

"""
    catstate(basis::SymplecticBasis, α::Number, phase::Real=0; squeeze_params=nothing, ħ=2)

Create a general Schrödinger cat state `|cat⟩ = (|α⟩ + e^(iφ)|-α⟩)/√N` where `φ` is the relative phase.

This is the most general form of cat state, which reduces to even `(φ=0)` or odd `(φ=π)` cat states
for specific phase values.

# Arguments
- `basis::SymplecticBasis`: The symplectic basis for the state
- `α::Number`: Complex amplitude of the coherent states
- `phase::Real=0`: Relative phase between the two coherent state components
- `squeeze_params=nothing`: Optional tuple `(r, θ)` for squeezed cat states
- `ħ=2`: Reduced Planck constant

# Mathematical Description
`|cat⟩ = (|α⟩ + e^(iφ)|-α⟩)/√(2(1 + Re(e^(iφ)exp(-2|α|²))))`
"""
function catstate(basis::SymplecticBasis, α::Number, phase::Real=0; squeeze_params=nothing, ħ=2)
    if squeeze_params === nothing
        state_plus = coherentstate(basis, α, ħ=ħ)
        state_minus = coherentstate(basis, -α, ħ=ħ)
    else
        r, θ = squeeze_params
        vac = vacuumstate(basis, ħ=ħ)
        squeeze_op = squeeze(basis, r, θ, ħ=ħ)
        squeezed_vac = squeeze_op * vac
        displace_plus = displace(basis, α, ħ=ħ)
        displace_minus = displace(basis, -α, ħ=ħ)
        state_plus = displace_plus * squeezed_vac
        state_minus = displace_minus * squeezed_vac
    end
    overlap = exp(-2 * abs2(α))
    phase_factor = exp(1im * phase)
    norm_factor = 1 / sqrt(2 * (1 + real(phase_factor * overlap)))
    coeffs = [norm_factor, norm_factor * phase_factor]
    if abs(phase - π) < 1e-14 || abs(phase) < 1e-14
        coeffs = real.(coeffs)
    end
    states = [state_plus, state_minus]
    return GaussianLinearCombination(basis, coeffs, states)
end

"""
    gkpstate(basis::SymplecticBasis; lattice=:square, delta=0.1, nmax=5, ħ=2)

Create a Gottesman-Kitaev-Preskill (GKP) state as a finite-energy approximation using squeezed states.

GKP states are quantum error-correcting codes that protect against small displacement errors.
They are constructed as superpositions of squeezed states positioned at lattice points in phase space.

# Arguments
- `basis::SymplecticBasis`: The symplectic basis for the state
- `lattice=:square`: Lattice type (:square or :hexagonal)
- `delta=0.1`: Squeezing parameter for finite energy approximation (smaller = more squeezed)
- `nmax=5`: Maximum lattice index for truncation (controls number of peaks)
- `ħ=2`: Reduced Planck constant

# Notes
- Very small `delta` values (< 1e-6) may cause numerical instability
- Large `nmax` values (> 50) create many states and may impact performance
- For square lattice: creates `2*nmax + 1` total states

# Mathematical Description
For square lattice: `|GKP⟩ = Σₖ |xₖ⟩` where `xₖ = √(2πħ) × k` for integer `k ∈ [-nmax, nmax]`
Each `|xₖ⟩` is approximated by a squeezed state with squeezing delta in the conjugate direction.
"""
function gkpstate(basis::SymplecticBasis; lattice="square", delta=0.1, nmax=5, ħ=2)
    if lattice == "square"
        return _square(basis, delta, nmax, ħ)
    elseif lattice == "hexagonal"
        return _hexagonal(basis, delta, nmax, ħ)
    else
        throw(ArgumentError("lattice must be \"square\" or \"hexagonal\", got \"$lattice\""))
    end
end

"""
    _square(basis::SymplecticBasis, delta, nmax, ħ)

Internal function to create square lattice GKP state.
"""
function _square(basis::SymplecticBasis, delta, nmax, ħ)
    if delta <= 0
        throw(ArgumentError("delta must be positive, got $delta"))
    end
    if nmax <= 0
        throw(ArgumentError("nmax must be positive, got $nmax"))
    end
    lattice_spacing = sqrt(2 * π * ħ)
    lattice_points = [k * lattice_spacing for k in -nmax:nmax]
    n_points = length(lattice_points)
    states = Vector{GaussianState}(undef, n_points)
    coeffs = ones(Float64, n_points)
    for (i, x_pos) in enumerate(lattice_points)
        squeeze_op = squeeze(basis, delta, π/2, ħ=ħ)
        displace_op = displace(basis, x_pos, ħ=ħ)
        vac = vacuumstate(basis, ħ=ħ)
        squeezed_state = squeeze_op * vac
        final_state = displace_op * squeezed_state
        states[i] = final_state
    end
    lc = GaussianLinearCombination(basis, coeffs, states)
    norm_val = norm_factor(lc.states, lc.coeffs)
    lc.coeffs .= lc.coeffs .* norm_val
    return lc
end

"""
    _hexagonal(basis::SymplecticBasis, delta, nmax, ħ)

Internal function to create hexagonal lattice GKP state.
"""
function _hexagonal(basis::SymplecticBasis, delta, nmax, ħ)
    lattice_spacing = sqrt(2 * π * ħ / sqrt(3))  
    states = GaussianState[]
    coeffs = Float64[]
    for m in -nmax:nmax
        for n in -nmax:nmax
            if abs(m) + abs(n) > nmax
                continue
            end
            x_pos = lattice_spacing * (m + 0.5 * n)
            p_pos = lattice_spacing * (sqrt(3) / 2 * n)
            squeeze_op = squeeze(basis, delta, 0.0, ħ=ħ)  
            displacement = x_pos + 1im * p_pos
            displace_op = displace(basis, displacement, ħ=ħ)
            vac = vacuumstate(basis, ħ=ħ)
            squeezed_state = squeeze_op * vac
            final_state = displace_op * squeezed_state
            push!(states, final_state)
            push!(coeffs, 1.0)
        end
    end
    lc = GaussianLinearCombination(basis, coeffs, states)
    norm_val = norm_factor(lc.states, lc.coeffs)
    lc.coeffs .= lc.coeffs .* norm_val
    return lc
end

"""
    norm_factor(states::Vector{GaussianState}, coeffs::Vector{<:Number})

Calculate the normalization factor for a linear combination of Gaussian states.

This function computes the normalization constant needed to ensure `⟨ψ|ψ⟩ = 1` for a state
`|ψ⟩ = Σᵢ cᵢ|ψᵢ⟩` by calculating all overlap integrals between component states.

# Arguments
- `states::Vector{GaussianState}`: Vector of Gaussian states
- `coeffs::Vector{<:Number}`: Vector of coefficients

# Returns
- `Float64`: Normalization factor `N` such that the state `Σᵢ (cᵢ/N)|ψᵢ⟩` is normalized

# Notes
- Returns 1.0 if near-zero or negative normalization is detected (may indicate coefficient cancellation)
- Numerical instabilities can occur with very small coefficients or nearly orthogonal states

# Mathematical Description
`N² = Σᵢⱼ cᵢ* cⱼ ⟨ψᵢ|ψⱼ⟩`
"""
function norm_factor(states::Vector{<:GaussianState}, coeffs::Vector{<:Number})
    n = length(states)
    @assert length(coeffs) == n "Number of coefficients must match number of states"
    norm_squared = 0.0
    for i in 1:n
        for j in 1:n
            overlap = _overlap(states[i], states[j])
            norm_squared += real(conj(coeffs[i]) * coeffs[j] * overlap)
        end
    end
    if norm_squared <= 1e-15
        return 1.0
    end
    return 1.0 / sqrt(norm_squared)
end

"""
    _overlap(state1::GaussianState, state2::GaussianState)

Calculate the overlap `⟨ψ₁|ψ₂⟩` between two Gaussian states.

# Notes
- Returns 0.0 if numerical instabilities are encountered during calculation
- Handles edge cases like singular covariance matrices gracefully
"""
function _overlap(state1::GaussianState, state2::GaussianState)
    @assert state1.basis == state2.basis "States must have the same basis"
    @assert state1.ħ == state2.ħ "States must have the same ħ"
    if state1 === state2
        return ComplexF64(1.0)
    end
    if isapprox(state1.mean, state2.mean, atol=1e-12) && 
       isapprox(state1.covar, state2.covar, atol=1e-12)
        return ComplexF64(1.0)
    end
    μ1, μ2 = state1.mean, state2.mean
    V1, V2 = state1.covar, state2.covar
    Δμ = μ1 - μ2
    V_sum = V1 + V2
    try
        exp_factor = exp(-0.25 * dot(Δμ, V_sum \ Δμ))
        det_factor = sqrt(det(V1) * det(V2) / det(V_sum))
        overlap = exp_factor * det_factor
        return ComplexF64(overlap)
    catch e
        return ComplexF64(0.0)
    end
end

"""
    fidelity_approximation(ideal_gkp::GaussianLinearCombination, finite_gkp::GaussianLinearCombination)

Estimate the fidelity between an ideal GKP state and its finite-energy approximation.

The fidelity gives a measure of how well the finite approximation represents the ideal
infinite-energy GKP state.

# Arguments
- `ideal_gkp::GaussianLinearCombination`: Ideal GKP state (high nmax, small delta)
- `finite_gkp::GaussianLinearCombination`: Finite approximation

# Returns
- `Float64`: Fidelity estimate `F = |⟨ideal|finite⟩|²`

"""
function fidelity_approximation(ideal_gkp::GaussianLinearCombination, finite_gkp::GaussianLinearCombination)
    @assert ideal_gkp.basis == finite_gkp.basis "States must have the same basis"
    @assert ideal_gkp.ħ == finite_gkp.ħ "States must have the same ħ"
    overlap = 0.0 + 0.0im
    for (c1, s1) in ideal_gkp
        for (c2, s2) in finite_gkp
            overlap += conj(c1) * c2 * _overlap(s1, s2)
        end
    end
    return abs2(overlap)
end

"""
    catstate_even(basis::SymplecticBasis, αs::AbstractVector; squeeze_params=nothing, ħ=2)

Create multi-mode even cat states as tensor products of single-mode cat states.

# Arguments
- `basis::SymplecticBasis`: The symplectic basis for the multi-mode system
- `αs::AbstractVector`: Vector of complex amplitudes for each mode
- `squeeze_params=nothing`: Optional vector of tuples `(r, θ)` for each mode
- `ħ=2`: Reduced Planck constant
"""
function catstate_even(basis::SymplecticBasis, αs::AbstractVector; squeeze_params=nothing, ħ=2)
    nmodes = basis.nmodes
    @assert length(αs) == nmodes "Number of amplitudes must match number of modes"
    if squeeze_params !== nothing
        @assert length(squeeze_params) == nmodes "Number of squeeze parameters must match number of modes"
    end
    single_mode_basis = typeof(basis)(1)
    cat_states = []
    for i in 1:nmodes
        α = αs[i]
        squeeze_param = squeeze_params === nothing ? nothing : squeeze_params[i]
        cat = catstate_even(single_mode_basis, α, squeeze_params=squeeze_param, ħ=ħ)
        push!(cat_states, cat)
    end
    result = cat_states[1]
    for i in 2:nmodes
        result = _tensor(result, cat_states[i])
    end
    return result
end

"""
    catstate_odd(basis::SymplecticBasis, αs::AbstractVector; squeeze_params=nothing, ħ=2)

Create multi-mode odd cat states as tensor products of single-mode cat states.
"""
function catstate_odd(basis::SymplecticBasis, αs::AbstractVector; squeeze_params=nothing, ħ=2)
    nmodes = basis.nmodes
    @assert length(αs) == nmodes "Number of amplitudes must match number of modes"
    if squeeze_params !== nothing
        @assert length(squeeze_params) == nmodes "Number of squeeze parameters must match number of modes"
    end
    single_mode_basis = typeof(basis)(1)
    cat_states = []
    for i in 1:nmodes
        α = αs[i]
        squeeze_param = squeeze_params === nothing ? nothing : squeeze_params[i]
        cat = catstate_odd(single_mode_basis, α, squeeze_params=squeeze_param, ħ=ħ)
        push!(cat_states, cat)
    end
    result = cat_states[1]
    for i in 2:nmodes
        result = _tensor(result, cat_states[i])
    end
    return result
end

"""
    catstate(basis::SymplecticBasis, αs::AbstractVector, phases::AbstractVector=zeros(length(αs)); squeeze_params=nothing, ħ=2)

Create multi-mode general cat states as tensor products of single-mode cat states.
"""
function catstate(basis::SymplecticBasis, αs::AbstractVector, phases::AbstractVector=zeros(length(αs)); squeeze_params=nothing, ħ=2)
    nmodes = basis.nmodes
    @assert length(αs) == nmodes "Number of amplitudes must match number of modes"
    @assert length(phases) == nmodes "Number of phases must match number of modes"
    if squeeze_params !== nothing
        @assert length(squeeze_params) == nmodes "Number of squeeze parameters must match number of modes"
    end
    single_mode_basis = typeof(basis)(1)
    cat_states = []
    for i in 1:nmodes
        α = αs[i]
        phase = phases[i]
        squeeze_param = squeeze_params === nothing ? nothing : squeeze_params[i]
        cat = catstate(single_mode_basis, α, phase, squeeze_params=squeeze_param, ħ=ħ)
        push!(cat_states, cat)
    end
    result = cat_states[1]
    for i in 2:nmodes
        result = _tensor(result, cat_states[i])
    end
    return result
end

"""
    _tensor(lc1::GaussianLinearCombination, lc2::GaussianLinearCombination)

Internal function to compute tensor product of two GaussianLinearCombination objects.
"""
function _tensor(lc1::GaussianLinearCombination, lc2::GaussianLinearCombination)
    @assert typeof(lc1.basis) == typeof(lc2.basis) "Linear combinations must have compatible bases"
    @assert lc1.ħ == lc2.ħ "Linear combinations must have the same ħ"
    new_basis = lc1.basis ⊕ lc2.basis
    result_size = length(lc1) * length(lc2)
    CoeffType = promote_type(eltype(lc1.coeffs), eltype(lc2.coeffs))
    new_coeffs = Vector{CoeffType}(undef, result_size)
    new_states = Vector{GaussianState}(undef, result_size)
    idx = 1
    for (c1, s1) in lc1
        for (c2, s2) in lc2
            new_coeffs[idx] = c1 * c2
            new_states[idx] = s1 ⊗ s2
            idx += 1
        end
    end
    return GaussianLinearCombination(new_basis, new_coeffs, new_states)
end