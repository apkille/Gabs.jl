abstract type AbstractGeneralDyne end

struct Generaldyne{I} <: AbstractGeneralDyne
    system::GaussianState
    conditional::GaussianState
    indices::I
    function Generaldyne(sys::GaussianState, cond::GaussianState, ind::I) where {I}
        length(ind) == 2*(cond.nmodes) || throw(DimensionMismatch(GENERALDYNE_ERROR))
        return new{I}(sys, ind, cond)
    end
end
function probability(meas::Generaldyne)
    pure_mean, pure_covar = conditional.mean, conditional.covar
end

function outcome(meas::Generaldyne)
    sys_mean, sys_covar, sys_nmodes = system.mean, system.covar, system.nmodes
    cond_mean, cond_covar, cond_nmodes = conditional.mean, conditional.covar, conditional.nmodes
    nmodes′ = sys_nmodes - cond_nmodes
    ind = meas.indices

end

function _part_mean(sys::M, cond::C, ind::I) where {M,I}
    # initialize mean vectors of subsystems A and B
    sys_nmodes, nmodesB = sys.nmodes, cond.nmodes
    sys_covar = sys.covar
    nmodesA = sys_nmodes - nmodesB
    meanA, meanB = zeros(2*nmodesA), zeros(2*nmodesB)
    # loop through entire Gaussian system, writing quadratures to B if
    # index is specified in `ind` argument
    idxA, idxB = 1, 1
    @inbounds for i in Base.OneTo(2*sys_nmodes)
        if i in ind
            covarB[2*idxB-1:2*idxB, 2*idxB-1:2*idxB] = @view(sys_covar[2*i-1:2*i, 2*i-1:2*i])
            idxB += 1
        else
            covarA[2*idxA-1:2*idxA, 2*idxA-1:2*idxA] = @view(sys_covar[2*i-1:2*i, 2*i-1:2*i])
            idxA += 1
        end
    end
    return meanA, meanB
end
function _part_covar(sys::M, cond::C, ind::I) where {M,I}
    # initialize mean vectors of subsystems A and B
    sys_nmodes, nmodesB = sys.nmodes, cond.nmodes
    sys_covar = sys.covar
    nmodesA = sys_nmodes - nmodesB
    covarA, covarB = zeros(2*nmodesA, 2*nmodesA), zeros(2*nmodesB, 2*nmodesB)
    covarAB = zeros(2*nmodesA, 2*nmodesB)
    # loop through entire Gaussian system, writing quadratures to B if
    # index is specified in `ind` argument; if not, write to A
    idxA, idxB = 1, 1
    @inbounds for i in Base.OneTo(sys_nmodes)
        if i in ind
            meanB[2*idxB-1:2*idxB] = @view(sys_mean[2*i-1:2*i])
            idxB += 1
        else
            meanA[2*idxA-1:2*idxA] = @view(sys_mean[2*i-1:2*i])
            idxA += 1
        end
    end
    return meanA, meanB
end

struct Homodyne{I,A<:Real} <: AbstractGeneralDyne
    system::GaussianState
    indices::I
    angle::A
end
Homodyne(sys::GaussianState, ind::I) where {I} = Homodyne(sys, ind, 0.0)

struct Heterodyne{I,A<:Number} <: AbstractGeneralDyne
    system::GaussianState
    indices::I
    alpha::A
end