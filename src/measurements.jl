abstract type AbstractGeneralDyne end

struct Generaldyne{I} <: AbstractGeneralDyne
    system::GaussianState
    conditional::GaussianState
    indices::I
    function Generaldyne(sys::GaussianState, cond::GaussianState, ind::I) where {I}
        length(ind) == cond.nmodes || throw(DimensionMismatch(GENERALDYNE_ERROR))
        return new{I}(sys, cond, ind)
    end
end
function probability(meas::Generaldyne)

end
function outcome(meas::Generaldyne)
    sys, cond, ind = meas.system, meas.conditional, meas.indices
    meanA, meanB = _part_mean(sys, ind)
    meanbuf1, meanbuf2 = zeros(2*cond.nmodes), zeros(2*(sys.nmodes - cond.nmodes))
    covarA, covarB, covarAB = _part_covar(sys, cond, ind)
    meanA .= meanA .- mul!(meanbuf2, covarAB, (mul!(meanbuf1, inv(covarB .+ cond.covar), cond.mean .- meanB)))
    covarA .= covarA .- covarAB * ((covarB .+ cond.covar) \ transpose(covarAB))
    meanA′ = _promote_output_vector(typeof(cond.mean), meanA, 2*(sys.nmodes - cond.nmodes))
    covarA′ = _promote_output_matrix(typeof(cond.covar), covarA, 2*(sys.nmodes - cond.nmodes))
    return GaussianState(meanA′, covarA′, sys.nmodes - cond.nmodes)
end

function _part_mean(sys::M, ind::I) where {M,I}
    # block mean into its modes
    mean′ = BlockedArray(sys.mean, 2*ones(Int,sys.nmodes))
    meanA = mortar([view(mean′, Block(i)) for i in Base.OneTo(sys.nmodes) if !(i in ind)])
    meanB = mortar([view(mean′, Block(i)) for i in ind])
    return meanA, meanB
end
function _part_covar(sys::M, cond::C, ind::I) where {M,C,I}
    sizeA, sizeB = (sys.nmodes-cond.nmodes, sys.nmodes-cond.nmodes), (cond.nmodes, cond.nmodes)
    sizeAB = (sys.nmodes-cond.nmodes, cond.nmodes)
    # loop through entire Gaussian system, writing quadratures to B if
    # index is specified in `ind` argument
    # write diagonal elements to subsystems
    covar′ = BlockedArray(sys.covar, 2*ones(Int,sys.nmodes), 2*ones(Int,sys.nmodes))
    covarA = mortar(reshape([view(covar′, Block(i,j))
            for i in Base.OneTo(sys.nmodes), j in Base.OneTo(sys.nmodes)
                if !(i in ind) && !(j in ind)], sizeA))
    covarB = mortar(reshape([view(covar′, Block(i,j))
            for i in Base.OneTo(sys.nmodes), j in Base.OneTo(sys.nmodes)
                if i in ind && j in ind], sizeB))
    covarAB = mortar(reshape([view(covar′, Block(j, i))
            for i in 1:sys.nmodes for j in 1:i if (i in ind && !(j in ind)) || (!(i in ind) && j in ind)], sizeAB))
    return covarA, covarB, covarAB
end