# Linear combinations of Gaussian states

"""
    GaussianLinearCombination{B<:SymplecticBasis,C,S}

Represents a linear combination of Gaussian states of the form Σᵢ cᵢ|ψᵢ⟩ where cᵢ are coefficients 
and |ψᵢ⟩ are Gaussian states, all sharing the same symplectic basis and ħ value.

## Fields
- `basis::B`: Symplectic basis shared by all states
- `coeffs::Vector{C}`: Complex coefficients for the linear combination  
- `states::Vector{S}`: Vector of Gaussian states
- `ħ::Number`: Reduced Planck's constant (must be same for all states)

## Examples

julia> basis = QuadPairBasis(1)
QuadPairBasis(1)

julia> state1 = coherentstate(basis, 1.0)
GaussianState for 1 mode.
  symplectic basis: QuadPairBasis
mean: 2-element Vector{Float64}:
 2.0
 0.0
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0

julia> lcgs1 = GaussianLinearCombination(state1)
GaussianLinearCombination with 1 terms for 1 mode.
  symplectic basis: QuadPairBasis
  ħ = 2
  [1] 1.0 * GaussianState

julia> state2 = coherentstate(basis, -1.0)
GaussianState for 1 mode.
  symplectic basis: QuadPairBasis
mean: 2-element Vector{Float64}:
 -2.0
  0.0
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0

julia> cat_state = 0.5 * GaussianLinearCombination(state1) + 0.5 * GaussianLinearCombination(state2)
GaussianLinearCombination with 2 terms for 1 mode.
  symplectic basis: QuadPairBasis
  ħ = 2
  [1] 0.5 * GaussianState
  [2] 0.5 * GaussianState

julia> lcgs2 = GaussianLinearCombination([0.6, 0.8], [state1, state2])
GaussianLinearCombination with 2 terms for 1 mode.
  symplectic basis: QuadPairBasis
  ħ = 2
  [1] 0.6 * GaussianState
  [2] 0.8 * GaussianState
"""
mutable struct GaussianLinearCombination{B<:SymplecticBasis,C,S}
    basis::B
    coeffs::Vector{C}
    states::Vector{S}
    ħ::Number
    function GaussianLinearCombination(basis::B, coeffs::Vector{C}, states::Vector{S}) where {B<:SymplecticBasis,C,S}
        length(coeffs) == length(states) || throw(DimensionMismatch("Number of coefficients ($(length(coeffs))) must match number of states ($(length(states)))"))
        isempty(states) && throw(ArgumentError("Cannot create an empty linear combination"))
        ħ = first(states).ħ
        @inbounds for (i, state) in enumerate(states)
            state isa GaussianState || throw(ArgumentError("Element $i is not a GaussianState: got $(typeof(state))"))
            if state.basis != basis
                throw(ArgumentError("State $i has incompatible basis: expected $(typeof(basis))($(basis.nmodes)), got $(typeof(state.basis))($(state.basis.nmodes))"))
            end
            if state.ħ != ħ
                throw(ArgumentError("State $i has different ħ: expected $ħ, got $(state.ħ)"))
            end
        end
        return new{B,C,S}(basis, coeffs, states, ħ)
    end
end
function GaussianLinearCombination(state::GaussianState{B,M,V}) where {B,M,V}
    coeff_type = float(real(promote_type(eltype(M), eltype(V))))
    return GaussianLinearCombination(state.basis, [one(coeff_type)], [state])
end
function GaussianLinearCombination(pairs::Vector{<:Tuple})
    isempty(pairs) && throw(ArgumentError("Cannot create an empty linear combination"))
    coeffs = [convert(Number, p[1]) for p in pairs]
    states = [p[2] for p in pairs]
    @inbounds for (i, state) in enumerate(states)
        state isa GaussianState || throw(ArgumentError("Element $i: second element must be a GaussianState"))
    end
    basis = first(states).basis
    return GaussianLinearCombination(basis, coeffs, states)
end
function GaussianLinearCombination(coeffs::Vector{<:Number}, states::Vector{<:GaussianState})
    isempty(states) && throw(ArgumentError("Cannot create an empty linear combination"))
    basis = first(states).basis
    return GaussianLinearCombination(basis, coeffs, states)
end
function GaussianLinearCombination(pairs::Pair{<:Number,<:GaussianState}...)
    isempty(pairs) && throw(ArgumentError("Cannot create an empty linear combination"))
    coeffs = [p.first for p in pairs]
    states = [p.second for p in pairs]
    basis = first(states).basis
    return GaussianLinearCombination(basis, coeffs, states)
end

"""
    +(lc1::GaussianLinearCombination, lc2::GaussianLinearCombination)

Add two linear combinations of Gaussian states. Both must have the same symplectic basis and ħ value.
"""
function Base.:+(lc1::GaussianLinearCombination{B}, lc2::GaussianLinearCombination{B}) where {B<:SymplecticBasis}
    lc1.basis == lc2.basis || throw(ArgumentError(SYMPLECTIC_ERROR))
    lc1.ħ == lc2.ħ || throw(ArgumentError(HBAR_ERROR))
    coeffs = vcat(lc1.coeffs, lc2.coeffs)
    states = vcat(lc1.states, lc2.states)
    return GaussianLinearCombination(lc1.basis, coeffs, states)
end
function Base.:+(lc1::GaussianLinearCombination{B1}, lc2::GaussianLinearCombination{B2}) where {B1<:SymplecticBasis,B2<:SymplecticBasis}
    throw(ArgumentError(SYMPLECTIC_ERROR))
end

"""
    -(lc1::GaussianLinearCombination, lc2::GaussianLinearCombination)

Subtract two linear combinations of Gaussian states.
"""
function Base.:-(lc1::GaussianLinearCombination{B}, lc2::GaussianLinearCombination{B}) where {B<:SymplecticBasis}
    return lc1 + (-1) * lc2
end
function Base.:-(lc1::GaussianLinearCombination{B1}, lc2::GaussianLinearCombination{B2}) where {B1<:SymplecticBasis,B2<:SymplecticBasis}
    throw(ArgumentError(SYMPLECTIC_ERROR))
end

"""
    -(lc::GaussianLinearCombination)

Negate a linear combination of Gaussian states.
"""
Base.:-(lc::GaussianLinearCombination) = (-1) * lc

"""
    *(α::Number, lc::GaussianLinearCombination)

Multiply a linear combination by a scalar from the left.
"""
function Base.:*(α::Number, lc::GaussianLinearCombination)
    new_coeffs = α .* lc.coeffs
    return GaussianLinearCombination(lc.basis, new_coeffs, copy(lc.states))
end

"""
    *(lc::GaussianLinearCombination, α::Number)

Multiply a linear combination by a scalar from the right.
"""
Base.:*(lc::GaussianLinearCombination, α::Number) = α * lc

"""
    Gabs.normalize!(lc::GaussianLinearCombination)

Normalize the coefficients of a linear combination in-place using L2 norm.
Note: Use Gabs.normalize! to avoid conflicts with LinearAlgebra.normalize!
"""
function normalize!(lc::GaussianLinearCombination)
    norm_val = sqrt(sum(abs2, lc.coeffs))
    if norm_val > 0
        lc.coeffs ./= norm_val
    end
    return lc
end

"""
    length(lc::GaussianLinearCombination)

Return the number of terms in the linear combination.
"""
Base.length(lc::GaussianLinearCombination) = length(lc.coeffs)

"""
    getindex(lc::GaussianLinearCombination, i::Integer)

Access the i-th term as a (coefficient, state) tuple.
"""
Base.getindex(lc::GaussianLinearCombination, i::Integer) = (lc.coeffs[i], lc.states[i])

"""
    iterate(lc::GaussianLinearCombination, state=1)

Iterate over (coefficient, state) pairs in the linear combination.
"""
function Base.iterate(lc::GaussianLinearCombination, state::Int=1)
    if state > length(lc)
        return nothing
    else
        return ((lc.coeffs[state], lc.states[state]), state + 1)
    end
end

"""
    Gabs.simplify!(lc::GaussianLinearCombination; atol::Real=1e-14)
Simplify a linear combination by:
1. Removing coefficients smaller than `atol`
2. Combining terms with identical states
3. Ensuring at least one term remains (with minimal coefficient vacuum state if all coefficients become negligible)
Returns the modified linear combination.
"""
function simplify!(lc::GaussianLinearCombination; atol::Real=1e-14)
    if isempty(lc.coeffs)
        return lc
    end
    keep_mask = abs.(lc.coeffs) .> atol
    if !any(keep_mask)
        vac = vacuumstate(lc.basis, ħ = lc.ħ)
        coeff_type = eltype(lc.coeffs)
        lc.coeffs = [coeff_type(atol)] 
        lc.states = [vac]
        return lc
    end
    coeffs = lc.coeffs[keep_mask]
    states = lc.states[keep_mask]
    unique_states = typeof(states[1])[]
    combined_coeffs = eltype(coeffs)[]
    sizehint!(unique_states, length(states))
    sizehint!(combined_coeffs, length(states))
    for (coeff, state) in zip(coeffs, states)
        existing_idx = findfirst(s -> isapprox(s, state, atol=1e-12), unique_states)
        if existing_idx === nothing
            push!(unique_states, state)
            push!(combined_coeffs, coeff)
        else
            combined_coeffs[existing_idx] += coeff
        end
    end
    final_mask = abs.(combined_coeffs) .> atol
    if !any(final_mask)
        vac = vacuumstate(lc.basis, ħ = lc.ħ)
        coeff_type = eltype(combined_coeffs)
        lc.coeffs = [coeff_type(atol)]
        lc.states = [vac]
    else
        lc.coeffs = combined_coeffs[final_mask]
        lc.states = unique_states[final_mask]
    end
    return lc
end

function Base.show(io::IO, mime::MIME"text/plain", lc::GaussianLinearCombination)
    basis_name = nameof(typeof(lc.basis))
    nmodes = lc.basis.nmodes
    print(io, "GaussianLinearCombination with $(length(lc)) terms")
    if nmodes == 1
        println(io, " for 1 mode.")
    else
        println(io, " for $(nmodes) modes.")
    end
    println(io, "  symplectic basis: ", basis_name)
    println(io, "  ħ = $(lc.ħ)")
    max_display = min(length(lc), 5)
    for i in 1:max_display
        coeff, state = lc[i]
        println(io, "  [$i] $(coeff) * GaussianState")
    end
    if length(lc) > max_display
        println(io, "  ⋮ ($(length(lc) - max_display) more terms)")
    end
end

Base.show(io::IO, lc::GaussianLinearCombination) = print(io, "GaussianLinearCombination($(length(lc)) terms)")

function Base.:(==)(lc1::GaussianLinearCombination, lc2::GaussianLinearCombination)
    lc1.basis == lc2.basis || return false
    lc1.ħ == lc2.ħ || return false
    length(lc1) == length(lc2) || return false
    return lc1.coeffs == lc2.coeffs && all(isequal(s1, s2) for (s1, s2) in zip(lc1.states, lc2.states))
end

function Base.isapprox(lc1::GaussianLinearCombination, lc2::GaussianLinearCombination; kwargs...)
    lc1.basis == lc2.basis || return false
    lc1.ħ == lc2.ħ || return false
    length(lc1) == length(lc2) || return false
    return isapprox(lc1.coeffs, lc2.coeffs; kwargs...) && 
           all(isapprox(s1, s2; kwargs...) for (s1, s2) in zip(lc1.states, lc2.states))
end

"""
    *(op::GaussianUnitary, lc::GaussianLinearCombination)

Apply a Gaussian unitary operation to a linear combination of Gaussian states.
The unitary is applied to each component state while preserving coefficients.
"""
function Base.:(*)(op::GaussianUnitary, lc::GaussianLinearCombination)
    op.basis == lc.basis || throw(ArgumentError(ACTION_ERROR))
    op.ħ == lc.ħ || throw(ArgumentError(HBAR_ERROR))
    
    new_states = [op * state for state in lc.states]
    return GaussianLinearCombination(lc.basis, copy(lc.coeffs), new_states)
end

"""
    *(op::GaussianChannel, lc::GaussianLinearCombination)

Apply a Gaussian channel to a linear combination of Gaussian states.
The channel is applied to each component state while preserving coefficients.
"""
function Base.:(*)(op::GaussianChannel, lc::GaussianLinearCombination)
    op.basis == lc.basis || throw(ArgumentError(ACTION_ERROR))
    op.ħ == lc.ħ || throw(ArgumentError(HBAR_ERROR))
    
    new_states = [op * state for state in lc.states]
    return GaussianLinearCombination(lc.basis, copy(lc.coeffs), new_states)
end

"""
    tensor(::Type{Tc}, ::Type{Ts}, lc1::GaussianLinearCombination, lc2::GaussianLinearCombination)

Compute tensor product of two linear combinations with specified output types.
Creates all pairwise tensor products: Σᵢⱼ cᵢcⱼ |ψᵢ⟩ ⊗ |ϕⱼ⟩.
"""
function tensor(::Type{Tm}, ::Type{Tc}, lc1::GaussianLinearCombination, lc2::GaussianLinearCombination) where {Tm,Tc}
    typeof(lc1.basis) == typeof(lc2.basis) || throw(ArgumentError(SYMPLECTIC_ERROR))
    lc1.ħ == lc2.ħ || throw(ArgumentError(HBAR_ERROR))
    new_basis = lc1.basis ⊕ lc2.basis
    n1, n2 = length(lc1), length(lc2)
    CoeffType = promote_type(eltype(lc1.coeffs), eltype(lc2.coeffs))
    new_coeffs = Vector{CoeffType}(undef, n1 * n2)
    new_states = Vector{GaussianState}(undef, n1 * n2)
    @inbounds for i in 1:n1
        for j in 1:n2
            idx = (i-1) * n2 + j
            new_coeffs[idx] = lc1.coeffs[i] * lc2.coeffs[j]
            new_states[idx] = tensor(Tm, Tc, lc1.states[i], lc2.states[j])
        end
    end
    return GaussianLinearCombination(new_basis, new_coeffs, new_states)
end

"""
    tensor(lc1::GaussianLinearCombination, lc2::GaussianLinearCombination)

Compute tensor product of two linear combinations of Gaussian states.
"""
function tensor(lc1::GaussianLinearCombination, lc2::GaussianLinearCombination)
    typeof(lc1.basis) == typeof(lc2.basis) || throw(ArgumentError(SYMPLECTIC_ERROR))
    lc1.ħ == lc2.ħ || throw(ArgumentError(HBAR_ERROR))
    new_basis = lc1.basis ⊕ lc2.basis
    n1, n2 = length(lc1), length(lc2)
    CoeffType = promote_type(eltype(lc1.coeffs), eltype(lc2.coeffs))
    new_coeffs = Vector{CoeffType}(undef, n1 * n2)
    new_states = Vector{GaussianState}(undef, n1 * n2)
    @inbounds for i in 1:n1
        for j in 1:n2
            idx = (i-1) * n2 + j
            new_coeffs[idx] = lc1.coeffs[i] * lc2.coeffs[j]
            new_states[idx] = lc1.states[i] ⊗ lc2.states[j]
        end
    end
    return GaussianLinearCombination(new_basis, new_coeffs, new_states)
end

function ptrace(::Type{Tm}, ::Type{Tc}, lc::GaussianLinearCombination, index::Int) where {Tm,Tc}
    length([index]) < lc.basis.nmodes || throw(ArgumentError(INDEX_ERROR))
    traced_states = [ptrace(Tm, Tc, state, [index]) for state in lc.states]
    result = GaussianLinearCombination(traced_states[1].basis, copy(lc.coeffs), traced_states)
    simplify!(result)
    return result
end

function ptrace(::Type{Tm}, ::Type{Tc}, lc::GaussianLinearCombination, indices::AbstractVector{<:Int}) where {Tm,Tc}
    length(indices) < lc.basis.nmodes || throw(ArgumentError(INDEX_ERROR))
    traced_states = [ptrace(Tm, Tc, state, indices) for state in lc.states]
    result = GaussianLinearCombination(traced_states[1].basis, copy(lc.coeffs), traced_states)
    simplify!(result)
    return result
end

"""
    ptrace(lc::GaussianLinearCombination, index::Int)

Compute partial trace of a linear combination over specified index.
"""
function ptrace(lc::GaussianLinearCombination, index::Int)
    length([index]) < lc.basis.nmodes || throw(ArgumentError(INDEX_ERROR))
    traced_states = [ptrace(state, [index]) for state in lc.states]
    result = GaussianLinearCombination(traced_states[1].basis, copy(lc.coeffs), traced_states)
    simplify!(result)
    return result
end

"""
    ptrace(lc::GaussianLinearCombination, indices::AbstractVector{<:Int})

Compute partial trace of a linear combination over specified indices.
Combines identical traced states automatically.
"""
function ptrace(lc::GaussianLinearCombination, indices::AbstractVector{<:Int})
    length(indices) < lc.basis.nmodes || throw(ArgumentError(INDEX_ERROR))
    traced_states = [ptrace(state, indices) for state in lc.states]
    result = GaussianLinearCombination(traced_states[1].basis, copy(lc.coeffs), traced_states)
    simplify!(result)
    return result
end

"""
    cross_wigner(state1::GaussianState, state2::GaussianState, x::AbstractVector)

Compute the off-diagonal Wigner kernel (cross-Wigner function) between two Gaussian states.

The cross-Wigner function is given by:

    W₁₂(x) = (1/(2π)ⁿ√det((V₁+V₂)/2)) × 
             exp[-½(x-μ̄)ᵀ((V₁+V₂)/2)⁻¹(x-μ̄)] × 
             exp[i(μ₁-μ₂)ᵀΩ(x-μ̄)/ħ]

where:
- μ̄ = (μ₁ + μ₂)/2 is the average of the two mean vectors
- μ₁, μ₂ are the mean vectors of the two states  
- V₁, V₂ are the covariance matrices of the two states
- Ω is the symplectic form matrix
- n is the number of modes
- ħ is the reduced Planck constant

This function captures the quantum interference between two Gaussian states and is 
essential for computing Wigner functions of superposition states. The cross-Wigner 
function appears in the interference terms when computing the Wigner function of 
a linear combination |ψ⟩ = c₁|ψ₁⟩ + c₂|ψ₂⟩:

    W(x) = |c₁|²W₁(x) + |c₂|²W₂(x) + 2Re(c₁*c₂*W₁₂(x))

The normalization is chosen to ensure the identity property W₁₁(x) = W₁(x).

# Arguments
- `state1::GaussianState`: First Gaussian state
- `state2::GaussianState`: Second Gaussian state  
- `x::AbstractVector`: Phase space point where to evaluate the function

# Returns
- `ComplexF64`: Complex value of the cross-Wigner function at point x

# Notes
- The function is Hermitian: W₁₂(x) = W₂₁*(x)
- For identical states: W₁₁(x) = W₁(x) (reduces to regular Wigner function)
- The implementation uses log-space arithmetic for numerical stability
"""
function cross_wigner(state1::GaussianState,
    state2::GaussianState,
    x::AbstractVector)
    μ1, μ2 = state1.mean, state2.mean
    V1, V2 = state1.covar, state2.covar
    n = length(μ1) ÷ 2
    ħ = state1.ħ
    Ω = symplecticform(state1.basis)
    μ̄ = 0.5*(μ1 .+ μ2)
    Δμ = μ1 .- μ2
    Vavg = 0.5 * (V1 .+ V2)
    lognorm = -n*log(2π) - 0.5*logdet(Vavg)
    norm = exp(lognorm)
    dx = x .- μ̄
    gauss = exp(-0.5 * dot(dx, Vavg \ dx))
    phase = cis(dot(Δμ, Ω*dx) / ħ)
    return norm * gauss * phase
end

"""
    wigner(lc::GaussianLinearCombination, x::AbstractVector)

Compute Wigner function of a linear combination including quantum interference.
W(x) = Σᵢ |cᵢ|² Wᵢ(x) + 2 Σᵢ<ⱼ Re(cᵢ*cⱼ W_cross(ψᵢ,ψⱼ,x))
"""
function wigner(lc::GaussianLinearCombination, x::AbstractVector)
    length(x) == length(lc.states[1].mean) || throw(ArgumentError(WIGNER_ERROR))
    result = 0.0
    for (c, state) in lc
        result += abs2(c) * wigner(state, x)
    end
    for i in 1:length(lc)
        ci, si = lc[i]
        for j in (i+1):length(lc)
            cj, sj = lc[j]
            cross_term = 2 * real(conj(ci) * cj * cross_wigner(si, sj, x))
            result += cross_term
        end
    end
    return result
end

"""
    cross_wignerchar(state1::GaussianState, state2::GaussianState, xi::AbstractVector)

Compute cross-Wigner characteristic function between two Gaussian states.
"""
function cross_wignerchar(state1::GaussianState, state2::GaussianState, xi::AbstractVector)
    state1.basis == state2.basis || throw(ArgumentError(SYMPLECTIC_ERROR))
    state1.ħ == state2.ħ || throw(ArgumentError(HBAR_ERROR))
    length(xi) == length(state1.mean) || throw(ArgumentError(WIGNER_ERROR))
    μ1, μ2 = state1.mean, state2.mean
    V1, V2 = state1.covar, state2.covar
    μ12 = (μ1 + μ2) / 2
    V12 = (V1 + V2) / 2
    Omega = symplecticform(state1.basis)
    arg1 = -0.5 * dot(xi, (Omega * V12 * transpose(Omega)) * xi)
    arg2 = 1im * dot(Omega * μ12, xi)
    return exp(arg1 - arg2)
end

"""
    wignerchar(lc::GaussianLinearCombination, xi::AbstractVector)

Compute Wigner characteristic function of a linear combination including interference.
"""
function wignerchar(lc::GaussianLinearCombination, xi::AbstractVector)
    length(xi) == length(lc.states[1].mean) || throw(ArgumentError(WIGNER_ERROR))
    result = 0.0 + 0.0im
    for (c, state) in lc
        result += abs2(c) * wignerchar(state, xi)
    end
    for i in 1:length(lc)
        ci, si = lc[i]
        for j in (i+1):length(lc)
            cj, sj = lc[j]
            cross_term = 2 * real(conj(ci) * cj * cross_wignerchar(si, sj, xi))
            result += cross_term
        end
    end
    return result
end

function purity(lc::GaussianLinearCombination)
    # A GaussianLinearCombination represents a pure superposition state
    # |ψ⟩ = Σᵢ cᵢ|ψᵢ⟩, so purity = Tr(ρ²) = 1
    return 1.0
end

function entropy_vn(lc::GaussianLinearCombination)
    # A GaussianLinearCombination represents a pure superposition state
    # |ψ⟩ = Σᵢ cᵢ|ψᵢ⟩, so S(ρ) = -Tr(ρ log ρ) = 0
    return 0.0
end

"""
    measurement_probability(lc::GaussianLinearCombination, measurement::GaussianState, indices)

Calculate measurement probability using Born rule: P = |⟨measurement|Tr_complement(lc)⟩|².
"""
function measurement_probability(lc::GaussianLinearCombination, measurement::GaussianState, indices::Union{Int, AbstractVector{<:Int}})
    indices_vec = indices isa Int ? [indices] : collect(indices)
    expected_modes = length(indices_vec)
    measurement.basis.nmodes == expected_modes || throw(ArgumentError(GENERALDYNE_ERROR))
    lc.ħ == measurement.ħ || throw(ArgumentError(HBAR_ERROR))
    norm_squared = 0.0
    for i in 1:length(lc)
        ci = lc.coeffs[i]
        si = lc.states[i]
        for j in 1:length(lc)
            cj = lc.coeffs[j]
            sj = lc.states[j]
            overlap = _gaussian_overlap(si, sj)
            norm_squared += real(conj(ci) * cj * overlap)
        end
    end
    if norm_squared < 1e-15
        return 0.0
    end
    normalized_coeffs = lc.coeffs ./ sqrt(norm_squared)
    complement_indices = setdiff(1:lc.basis.nmodes, indices_vec)
    if isempty(complement_indices)
        lc_measured_states = lc.states
        lc_measured_coeffs = normalized_coeffs
    else
        lc_measured_states = [ptrace(state, complement_indices) for state in lc.states]
        lc_measured_coeffs = normalized_coeffs
    end
    overlap = 0.0 + 0.0im
    for i in 1:length(lc_measured_states)
        c = lc_measured_coeffs[i]
        state = lc_measured_states[i]
        state_overlap = _gaussian_overlap(measurement, state)
        overlap += c * state_overlap
    end
    prob = abs2(overlap)
    return clamp(prob, 0.0, 1.0)
end

function tensor(::Type{T}, lc1::GaussianLinearCombination, lc2::GaussianLinearCombination) where {T}
    if T <: AbstractMatrix
        return tensor(Vector{eltype(T)}, T, lc1, lc2)
    else
        return tensor(T, T, lc1, lc2)
    end
end

function ptrace(::Type{T}, lc::GaussianLinearCombination, indices::Union{Int, AbstractVector{<:Int}}) where {T}
    if T <: AbstractMatrix
        return ptrace(Vector{eltype(T)}, T, lc, indices)
    else
        return ptrace(T, T, lc, indices)
    end
end

"""
    coherence_measure(lc::GaussianLinearCombination)

Calculate how much the overlaps between component states affect the 
coherence of the quantum superposition. Returns value between 0 and 1,
where 1 means perfectly orthogonal states (maximum coherence) and
values < 1 indicate overlapping states that reduce effective coherence.
"""
function coherence_measure(lc::GaussianLinearCombination)
    if length(lc) == 1
        return 1.0
    end
    n = length(lc)
    overlap_matrix = Matrix{ComplexF64}(undef, n, n)
    for i in 1:n
        for j in 1:n
            overlap_matrix[i, j] = _gaussian_overlap(lc.states[i], lc.states[j])
        end
    end
    eigenvals = real(eigvals(Hermitian(overlap_matrix)))
    eigenvals = eigenvals[eigenvals .> 1e-15]
    if isempty(eigenvals)
        return 0.0
    end
    participation_ratio = (sum(eigenvals))^2 / sum(eigenvals.^2)
    effective_coherence = participation_ratio / n
    return min(effective_coherence, 1.0)
end

#additional functions - not a part of project i coded them because of a misunderstadning, but could be useful ?:
"""
    representation_purity(lc::GaussianLinearCombination)

Calculate the effective purity of the linear combination representation, accounting 
for overlaps between component states.

This function computes how "pure" the representation appears when considering the 
overlaps between component Gaussian states. Unlike quantum purity (which is always 1 
for pure superposition states), this measures the representational efficiency.

The formula used is:
    P_rep = |Σᵢⱼ cᵢ* cⱼ ⟨ψᵢ|ψⱼ⟩|² / |Σᵢⱼ cᵢ* cⱼ ⟨ψᵢ|ψⱼ⟩|²

Returns values between 0 and 1, where:
- 1.0: All component states are orthogonal (most efficient representation)
- < 1.0: States have overlaps (redundant representation)

# Note
This is NOT the quantum mechanical purity of the state (which is always 1 for 
pure superpositions), but rather a measure of the representation's efficiency.

# Returns
- `Float64`: Representation purity between 0 and 1
"""
function representation_purity(lc::GaussianLinearCombination)
    n = length(lc)
    if n == 1
        return 1.0
    end
    norm_squared = complex(0.0)
    for i in 1:n
        ci = lc.coeffs[i]
        state_i = lc.states[i]
        for j in 1:n
            cj = lc.coeffs[j] 
            state_j = lc.states[j]
            overlap_ij = _gaussian_overlap(state_i, state_j)
            norm_squared += conj(ci) * cj * overlap_ij
        end
    end
    norm_real = real(norm_squared)
    if norm_real < 1e-15
        return 0.0
    end
    tr_rho_squared = complex(0.0)
    for i in 1:n
        ci = lc.coeffs[i]
        state_i = lc.states[i]
        for j in 1:n
            cj = lc.coeffs[j]
            state_j = lc.states[j]
            for k in 1:n
                ck = lc.coeffs[k]
                state_k = lc.states[k]
                for l in 1:n
                    cl = lc.coeffs[l]
                    state_l = lc.states[l]
                    overlap_ik = _gaussian_overlap(state_i, state_k)
                    overlap_lj = _gaussian_overlap(state_l, state_j)
                    tr_rho_squared += conj(ci) * cj * conj(ck) * cl * overlap_ik * overlap_lj
                end
            end
        end
    end
    purity_value = real(tr_rho_squared) / abs2(norm_squared)
    return clamp(purity_value, 0.0, 1.0)
end

"""
    representation_entropy(lc::GaussianLinearCombination)

Calculate the informational entropy of the linear combination representation,
considering overlaps between component states.

This function measures the "informational complexity" of representing the quantum 
state as a linear combination of Gaussian states. It accounts for redundancy 
introduced by overlapping component states.

The calculation constructs an effective density matrix from state overlaps:
    ρᵢⱼ = cᵢ* cⱼ ⟨ψᵢ|ψⱼ⟩
and computes S = -Tr(ρ log ρ) from its eigenvalues.

Returns values ≥ 0, where:
- 0.0: All component states are orthogonal (minimal redundancy)
- > 0.0: States overlap (representational redundancy)

# Note  
This is NOT the von Neumann entropy of the quantum state (which is always 0 for
pure superpositions), but rather a measure of the representation's complexity.

# Applications
- Analyzing numerical stability of linear combinations
- Studying representational efficiency  
- Understanding classical simulation complexity

# Returns
- `Float64`: Representation entropy ≥ 0
"""
function representation_entropy(lc::GaussianLinearCombination)
    n = length(lc)
    if n == 1
        return 0.0  
    end
    total_norm = sqrt(sum(abs2, lc.coeffs))
    if total_norm < 1e-15
        return 0.0
    end
    normalized_coeffs = lc.coeffs ./ total_norm
    ρ = Matrix{ComplexF64}(undef, n, n)
    for i in 1:n
        ci = normalized_coeffs[i]
        si = lc.states[i]
        for j in 1:n
            cj = normalized_coeffs[j]
            sj = lc.states[j]
            overlap = _gaussian_overlap(si, sj)
            ρ[i, j] = conj(ci) * cj * overlap
        end
    end
    if n > 100
        @warn "Computing Von Neumann entropy for large system (n=$n). " *
            "This may be computationally expensive."
    end
    eigenvals = real(eigvals(Hermitian(ρ)))
    eigenvals = eigenvals[eigenvals .> 1e-15] 
    total = sum(eigenvals)
    if total > 1e-15
        eigenvals ./= total
    end
    entropy = 0.0
    for λ in eigenvals
        if λ > 1e-15
            entropy -= λ * log(λ)
        end
    end
    return max(entropy, 0.0)
end


