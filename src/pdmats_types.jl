"""
PDMats integration type definitions for GABS.jl
"""

using PDMats
using LinearAlgebra: Symmetric, tr, I, PosDefException

"""
Modified GaussianState structure that supports PDMats for more efficient operations
with covariance matrices.
"""
@kwdef struct PDGaussianState{B<:SymplecticBasis,M,V<:AbstractPDMat} <: StateVector{M,V}
    basis::B
    mean::M
    covar::V
    ħ::Number = 2
end

function PDGaussianState{B,M,V}(basis::B, mean::M, covar::V; ħ::Number = 2) where {B<:SymplecticBasis,M,V<:AbstractPDMat}
    size(covar, 1) == length(mean) == 2*(basis.nmodes) || throw(DimensionMismatch(STATE_ERROR))
    return PDGaussianState{B,M,V}(basis=basis, mean=mean, covar=covar, ħ=ħ)
end

function PDGaussianState(basis::B, mean::M, covar::V; ħ::Number = 2) where {B<:SymplecticBasis,M,V<:AbstractPDMat}
    size(covar, 1) == length(mean) == 2*(basis.nmodes) || throw(DimensionMismatch(STATE_ERROR))
    return PDGaussianState{B,M,V}(basis=basis, mean=mean, covar=covar, ħ=ħ)
end

function PDGaussianState(state::GaussianState)
    covar_sym = Symmetric(state.covar)
    pdmat = try
        PDMat(covar_sym)
    catch e
        if isa(e, PosDefException)
            eps_val = 1e-10 * tr(state.covar) / size(state.covar, 1)
            PDMat(covar_sym + eps_val * I)
        else
            rethrow(e)
        end
    end
    return PDGaussianState(state.basis, state.mean, pdmat, ħ=state.ħ)
end

function PDGaussianState(basis::B, mean::M, covar::AbstractMatrix; ħ::Number = 2) where {B<:SymplecticBasis,M}
    covar_sym = Symmetric(covar)
    pdmat = try
        PDMat(covar_sym)
    catch e
        if isa(e, PosDefException)
            eps_val = 1e-10 * tr(covar) / size(covar, 1)
            PDMat(covar_sym + eps_val * I)
        else
            rethrow(e)
        end
    end
    return PDGaussianState(basis, mean, pdmat, ħ=ħ)
end