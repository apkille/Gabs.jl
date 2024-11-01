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
function output(meas::Generaldyne)
    sys, cond, ind = meas.system, meas.conditional, meas.indices
    # partition mean and covar into subsystems A and B
    meanA, meanB = _part_mean(sys, ind)
    covarA, covarB, covarAB = _part_covar(sys, cond, ind)
    # Block array matmul and broadcasting is incredibly
    # slow, so convert types back to promoted `sys`` and `cond`` types
    meanA′, meanB′ = sys.mean isa Vector{Float64} ? Vector{Float64}.((meanA, meanB)) :
        (_promote_output_vector(typeof(cond.mean), meanA, 2*(sys.nmodes - cond.nmodes)), _promote_output_vector(typeof(cond.mean), meanB, 2*cond.nmodes))
    covarA′, covarB′, covarAB′ = sys.covar isa Matrix{Float64} ? Matrix{Float64}.((covarA, covarB, covarAB)) :
        (_promote_output_matrix(typeof(cond.covar), covarA, (2*(sys.nmodes - cond.nmodes), 2*(sys.nmodes - cond.nmodes))),
        _promote_output_matrix(typeof(cond.covar), covarB, (2*cond.nmodes, 2*cond.nmodes)),
        _promote_output_matrix(typeof(cond.covar), covarAB, (2*(sys.nmodes - cond.nmodes), 2*cond.nmodes)))
    # map subsystem A
    meanA′, covarA′ = _generaldyne_map(meanA′, meanB′, covarA′, covarB′, covarAB′, sys, cond)
    return GaussianState(meanA′, covarA′, sys.nmodes - cond.nmodes)
end

function _generaldyne_map(meanA, meanB, covarA, covarB, covarAB, sys, cond)
    # create alloc buffers for matrix multiplication
    meanbuf1, meanbuf2 = zeros(2*cond.nmodes), zeros(2*(sys.nmodes - cond.nmodes))
    covarbuf = zeros(2*(sys.nmodes - cond.nmodes), 2*(sys.nmodes - cond.nmodes))
    # maps subsystem A, which is not measured
    meanA .= meanA .+ mul!(meanbuf2, covarAB, (mul!(meanbuf1, inv(covarB .+ cond.covar), cond.mean .- meanB)))
    covarA .= covarA .- mul!(covarbuf, covarAB, (covarB .+ cond.covar) \ transpose(covarAB))
    return meanA, covarA
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