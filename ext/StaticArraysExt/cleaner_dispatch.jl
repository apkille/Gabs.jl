# Gaussian states with cleaner dispatch

# Vacuum state
function vacuumstate(::Type{SArray}, basis::SymplecticBasis{N}; ħ=2) where {N<:Int}
    n = basis.nmodes
    T = typeof(ħ/2)
    vacuumstate(SVector{2n,T}, SMatrix{2n,2n,T}, basis; ħ=ħ)
end
function vacuumstate(::Type{SVector}, ::Type{SMatrix}, basis::SymplecticBasis{N}; ħ=2) where {N<:Int}
    n = basis.nmodes
    T = typeof(ħ/2)
    vacuumstate(SVector{2n,T}, SMatrix{2n,2n,T}, basis; ħ=ħ)
end
function vacuumstate(::Type{SVector{M,T1}}, ::Type{SMatrix{M,M,T2}}, basis::SymplecticBasis{N}; ħ=2) where {M,N<:Int,T1,T2}
    n = basis.nmodes
    M == 2n || error("Size mismatch: SVector{$M}/SMatrix{$M,$M} != 2n (n=$n)")
    T = promote_type(T1, T2, typeof(ħ/2))
    GaussianState(basis, SVector{2n,T}(zeros(T, 2n)), SMatrix{2n,2n,T}((ħ/2) * I); ħ=ħ)
end

# Thermal state
function thermalstate(::Type{SArray}, basis::SymplecticBasis{N}, photons; ħ=2) where {N<:Int}
    n = basis.nmodes
    T = promote_type(typeof(ħ/2), eltype(photons))
    thermalstate(SVector{2n,T}, SMatrix{2n,2n,T}, basis, photons; ħ=ħ)
end
function thermalstate(::Type{SVector}, ::Type{SMatrix}, basis::SymplecticBasis{N}, photons; ħ=2) where {N<:Int}
    n = basis.nmodes
    T = promote_type(typeof(ħ/2), eltype(photons))
    thermalstate(SVector{2n,T}, SMatrix{2n,2n,T}, basis, photons; ħ=ħ)
end
function thermalstate(::Type{SVector{M,T1}}, ::Type{SMatrix{M,M,T2}}, basis::SymplecticBasis{N}, photons; ħ=2) where {M,N<:Int,T1,T2}
    n = basis.nmodes
    M == 2n || error("Size mismatch: SVector{$M}/SMatrix{$M,$M} != 2n (n=$n)")
    T = promote_type(T1, T2, typeof(ħ/2), eltype(photons))
    mean = zeros(SVector{2n,T})
    covar = (2 * photons + 1) * (ħ/2) * one(SMatrix{2n,2n,T})
    GaussianState(basis, mean, covar; ħ=ħ)
end

# Coherent state

function coherentstate(::Type{SArray}, basis::SymplecticBasis{N}, alpha; ħ=2) where {N<:Int}
    n = basis.nmodes
    T = promote_type(typeof(ħ/2), real(eltype(alpha)))
    coherentstate(SVector{2n,T}, SMatrix{2n,2n,T}, basis, alpha; ħ=ħ)
end
function coherentstate(::Type{SVector}, ::Type{SMatrix}, basis::SymplecticBasis{N}, alpha; ħ=2) where {N<:Int}
    n = basis.nmodes
    T = promote_type(typeof(ħ/2), real(eltype(alpha)))
    coherentstate(SVector{2n,T}, SMatrix{2n,2n,T}, basis, alpha; ħ=ħ)
end
function _complex_to_real_vec(alpha::Number, ħ, T, n)
    real_part = sqrt(2ħ) * real(alpha)
    imag_part = sqrt(2ħ) * imag(alpha)
    SVector{2n,T}([repeat([real_part, imag_part], n)...])
end
function _complex_to_real_vec(alpha::AbstractVector, ħ, T, n)
    length(alpha) == n || error("Number of complex amplitudes ($(length(alpha))) must match number of modes ($n)")
    SVector{2n,T}(sqrt(2ħ) * vcat(real.(alpha), imag.(alpha)))
end
function coherentstate(::Type{SVector{M,T1}}, ::Type{SMatrix{M,M,T2}}, basis::SymplecticBasis{N}, alpha; ħ=2) where {M,N<:Int,T1,T2}
    n = basis.nmodes
    M == 2n || error("Size mismatch: SVector{$M}/SMatrix{$M,$M} != 2n (n=$n)")
    T = promote_type(T1, T2, typeof(ħ/2), real(eltype(alpha)))
    mean = _complex_to_real_vec(alpha, ħ, T, n)
    covar = (ħ/2) * one(SMatrix{2n,2n,T})
    GaussianState(basis, mean, covar; ħ=ħ)
end

# Squeezed state
function squeezedstate(::Type{SArray}, basis::SymplecticBasis{N}, r, theta; ħ=2) where {N<:Int}
    n = basis.nmodes
    T = promote_type(typeof(ħ/2), eltype(r), eltype(theta))
    squeezedstate(SVector{2n,T}, SMatrix{2n,2n,T}, basis, r, theta; ħ=ħ)
end
function squeezedstate(::Type{SVector}, ::Type{SMatrix}, basis::SymplecticBasis{N}, r, theta; ħ=2) where {N<:Int}
    n = basis.nmodes
    T = promote_type(typeof(ħ/2), eltype(r), eltype(theta))
    squeezedstate(SVector{2n,T}, SMatrix{2n,2n,T}, basis, r, theta; ħ=ħ)
end
function squeezedstate(::Type{SVector{M,T1}}, ::Type{SMatrix{M,M,T2}}, basis::SymplecticBasis{N}, r, theta; ħ=2) where {M,N<:Int,T1,T2}
    n = basis.nmodes
    M == 2n || error("Size mismatch: SVector{$M}/SMatrix{$M,$M} != 2n (n=$n)")
    T = promote_type(T1, T2, typeof(ħ/2), eltype(r), eltype(theta))
    mean = zeros(SVector{2n,T})
    covar = _squeezed_covar(SMatrix{2n,2n,T}, r, theta, ħ)
    GaussianState(basis, mean, covar; ħ=ħ)
end

function tensor(::Type{SVector}, ::Type{SMatrix}, state1::GaussianState, state2::GaussianState)
    M1, V1 = typeof(state1.mean), typeof(state1.covar)
    M2, V2 = typeof(state2.mean), typeof(state2.covar)
    out_size = 2*(state1.basis.nmodes + state2.basis.nmodes)
    tensor(
        SVector{out_size, promote_type(eltype(M1), eltype(M2))},
        SMatrix{out_size, out_size, promote_type(eltype(V1), eltype(V2))},
        state1,
        state2
    )
end
function tensor(::Type{SArray}, state1::GaussianState, state2::GaussianState)
    tensor(SVector, SMatrix, state1, state2)
end
function _output_size(basis::SymplecticBasis{N}, indices) where {N}
    nmodes = basis.nmodes
    notindices = setdiff(1:nmodes, indices)
    notidxlength = length(notindices)
    2 * notidxlength
end

function _output_types(state::GaussianState{B,M,V}) where {B,M,V}
    (M, V)
end

# ptrace
function ptrace(::Type{SVector}, ::Type{SMatrix}, state::GaussianState, indices)
    M, V = _output_types(state)
    out_size = _output_size(state.basis, indices)
    ptrace(SVector{out_size, eltype(M)}, SMatrix{out_size, out_size, eltype(V)}, state, indices)
end
function ptrace(::Type{SArray}, state::GaussianState, indices)
    M, V = _output_types(state)
    out_size = _output_size(state.basis, indices)
    ptrace(SVector{out_size, eltype(M)}, SMatrix{out_size, out_size, eltype(V)}, state, indices)
end
