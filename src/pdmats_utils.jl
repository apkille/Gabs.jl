"""
PDMats integration utilities for GABS.jl
"""

using PDMats
using LinearAlgebra
using SparseArrays: blockdiag

function chol_lower(a::Cholesky)
    return a.uplo === 'L' ? LowerTriangular(a.factors) : LowerTriangular(a.factors')
end

function vacuumstate_pd(basis::SymplecticBasis{N}; ħ = 2) where {N<:Int}
    mean = zeros(2*basis.nmodes)
    scal_mat = ScalMat(2*basis.nmodes, ħ/2)
    return PDGaussianState(basis, mean, scal_mat, ħ=ħ)
end

function thermalstate_pd(basis::SymplecticBasis{N}, photons::P; ħ = 2) where {N<:Int,P<:Number}
    mean = zeros(2*basis.nmodes)
    scal_mat = ScalMat(2*basis.nmodes, (2 * photons + 1) * (ħ/2))
    return PDGaussianState(basis, mean, scal_mat, ħ=ħ)
end

function thermalstate_pd(basis::QuadPairBasis{N}, photons::Vector{P}; ħ = 2) where {N<:Int,P<:Real}
    mean = zeros(2*basis.nmodes)
    diag_vals = Vector{Float64}(undef, 2*basis.nmodes)
    for i in 1:basis.nmodes
        val = (2 * photons[i] + 1) * (ħ/2)
        diag_vals[2*i-1] = val
        diag_vals[2*i] = val
    end
    diag_mat = PDiagMat(diag_vals)
    return PDGaussianState(basis, mean, diag_mat, ħ=ħ)
end

function coherentstate_pd(basis::SymplecticBasis{N}, alpha::A; ħ = 2) where {N<:Int,A}
    mean, _ = _coherentstate(basis, alpha; ħ = ħ)
    scal_mat = ScalMat(2*basis.nmodes, ħ/2)
    return PDGaussianState(basis, mean, scal_mat, ħ=ħ)
end

function whiten(state::PDGaussianState, x::AbstractVecOrMat)
    size(x, 1) == 2*state.basis.nmodes || throw(DimensionMismatch("Dimension mismatch in whitening operation"))
    return PDMats.whiten(state.covar, x)
end

function unwhiten(state::PDGaussianState, x::AbstractVecOrMat)
    size(x, 1) == 2*state.basis.nmodes || throw(DimensionMismatch("Dimension mismatch in unwhitening operation"))
    return PDMats.unwhiten(state.covar, x)
end

function whiten(state::GaussianState, x::AbstractVecOrMat)
    size(x, 1) == 2*state.basis.nmodes || throw(DimensionMismatch("Dimension mismatch in whitening operation"))
    
    chol_factor = cholesky(Symmetric(state.covar))
    return chol_lower(chol_factor) \ x
end

function unwhiten(state::GaussianState, x::AbstractVecOrMat)
    size(x, 1) == 2*state.basis.nmodes || throw(DimensionMismatch("Dimension mismatch in unwhitening operation"))
    
    chol_factor = cholesky(Symmetric(state.covar))
    return chol_lower(chol_factor) * x
end

function quad(state::PDGaussianState, x::AbstractVecOrMat)
    size(x, 1) == 2*state.basis.nmodes || throw(DimensionMismatch("Dimension mismatch in quadratic form"))
    return PDMats.quad(state.covar, x)
end

function invquad(state::PDGaussianState, x::AbstractVecOrMat)
    size(x, 1) == 2*state.basis.nmodes || throw(DimensionMismatch("Dimension mismatch in inverse quadratic form"))
    return PDMats.invquad(state.covar, x)
end

function quad(state::GaussianState, x::AbstractVecOrMat)
    size(x, 1) == 2*state.basis.nmodes || throw(DimensionMismatch("Dimension mismatch in quadratic form"))
    return quad(state.covar, x)
end

function invquad(state::GaussianState, x::AbstractVecOrMat)
    size(x, 1) == 2*state.basis.nmodes || throw(DimensionMismatch("Dimension mismatch in inverse quadratic form"))
    return invquad(state.covar, x)
end



using SparseArrays: blockdiag

function pd_blockdiag_cov(cov1::ScalMat, cov2::ScalMat)
    if cov1.value == cov2.value
        return ScalMat(cov1.dim + cov2.dim, cov1.value)
    else
        n1, n2 = cov1.dim, cov2.dim
        diagvec = vcat(fill(cov1.value, n1), fill(cov2.value, n2))
        return PDiagMat(diagvec)
    end
end

function pd_blockdiag_cov(cov1::PDiagMat, cov2::PDiagMat)
    new_diag = vcat(cov1.diag, cov2.diag)
    return PDiagMat(new_diag)
end

function pd_blockdiag_cov(cov1::AbstractPDMat, cov2::AbstractPDMat)
    return PDMat(blockdiag(Matrix(cov1), Matrix(cov2)))
end

function tensor(state1::PDGaussianState, state2::PDGaussianState)
    typeof(state1.basis) == typeof(state2.basis) || throw(ArgumentError(SYMPLECTIC_ERROR))
    state1.ħ == state2.ħ || throw(ArgumentError(HBAR_ERROR))
    
    mean1, mean2 = state1.mean, state2.mean
    basis1, basis2 = state1.basis, state2.basis
    nmodes1, nmodes2 = basis1.nmodes, basis2.nmodes
    nmodes = nmodes1 + nmodes2
    Mt = promote_type(eltype(mean1), eltype(mean2))
    mean_combined = zeros(Mt, 2 * nmodes)
    if basis1 isa QuadPairBasis
        for i in 1:(2*nmodes1)
            mean_combined[i] = mean1[i]
        end
        for i in 1:(2*nmodes2)
            mean_combined[i + 2*nmodes1] = mean2[i]
        end
    else 
        for i in 1:nmodes1
            mean_combined[i] = mean1[i]
            mean_combined[i+nmodes] = mean1[i+nmodes1]
        end
        for i in 1:nmodes2
            mean_combined[i+nmodes1] = mean2[i]
            mean_combined[i+nmodes+nmodes1] = mean2[i+nmodes2]
        end
    end

    cov_combined = pd_blockdiag_cov(state1.covar, state2.covar)
    new_basis = state1.basis ⊕ state2.basis
    return PDGaussianState(new_basis, mean_combined, cov_combined, ħ = state1.ħ)
end




function _ptrace(state::PDGaussianState{B, M, V}, idx::T) where {B<:QuadPairBasis, M, V<:AbstractPDMat, T<:Int}
    idxV = 2*idx - 1 : 2*idx
    mean_block = state.mean[idxV]
    new_cov = if state.covar isa ScalMat
        ScalMat(2, state.covar.value)
    elseif state.covar isa PDiagMat
        PDiagMat(state.covar.diag[idxV])
    else
        cov_dense = Matrix(state.covar)
        cov_block = @view cov_dense[idxV, idxV]
        PDMat(Symmetric(cov_block))
    end
    return mean_block, new_cov
end

function _ptrace(state::PDGaussianState{B, M, V}, idx::T) where {B<:QuadBlockBasis, M, V<:AbstractPDMat, T<:Int}
    nmodes = state.basis.nmodes
    mean_block = [ state.mean[idx], state.mean[idx + nmodes] ]
    new_cov = if state.covar isa ScalMat
        ScalMat(2, state.covar.value)
    elseif state.covar isa PDiagMat
        PDiagMat(state.covar.diag[[idx, idx + nmodes]])
    else
        cov_dense = Matrix(state.covar)
        A = cov_dense[idx, idx]
        B_ = cov_dense[idx, idx + nmodes]
        C = cov_dense[idx + nmodes, idx]
        D = cov_dense[idx + nmodes, idx + nmodes]
        PDMat(Symmetric([A B_; C D]))
    end
    return mean_block, new_cov
end


function ptrace(state::PDGaussianState, idx::Int)
    mean_block, cov_block = _ptrace(state, idx)
    new_basis = typeof(state.basis)(1)
    return PDGaussianState(new_basis, mean_block, cov_block, ħ=state.ħ)
end

function ptrace(state::PDGaussianState, indices::AbstractVector{Int})
    new_basis = typeof(state.basis)(length(indices))
    new_mean = zeros(eltype(state.mean), 2 * length(indices))
    cov_blocks = Vector{typeof(state.covar)}(undef, length(indices))
    for (i, idx) in enumerate(indices)
        m, c = _ptrace(state, idx)
        new_mean[2*i-1:2*i] .= m
        cov_blocks[i] = c
    end
    new_cov = pd_blockdiag_cov(cov_blocks...)
    return PDGaussianState(new_basis, new_mean, new_cov, ħ=state.ħ)
end

function apply!(state::PDGaussianState, op::GaussianUnitary)
    op.basis == state.basis || throw(DimensionMismatch(ACTION_ERROR))
    op.ħ == state.ħ || throw(ArgumentError(HBAR_ERROR))
    
    d, S = op.disp, op.symplectic
    state.mean .= S * state.mean .+ d
    
    temp_covar = S * Matrix(state.covar) * transpose(S)
    state.covar = PDMat(temp_covar)
    
    return state
end

function apply!(state::PDGaussianState, op::GaussianChannel)
    op.basis == state.basis || throw(DimensionMismatch(ACTION_ERROR))
    op.ħ == state.ħ || throw(ArgumentError(HBAR_ERROR))
    
    d, T, N = op.disp, op.transform, op.noise
    state.mean .= T * state.mean .+ d
    
    temp_covar = T * Matrix(state.covar) * transpose(T) .+ N
    state.covar = PDMat(temp_covar)
    
    return state
end