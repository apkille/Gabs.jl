function _generaldyne_map(meanA::SVector, meanB::SVector, covarA::SMatrix, covarB::SMatrix, covarAB::SMatrix, sys, cond)
    # maps subsystem A, which is not measured
    meanA′ = meanA + covarAB * (inv(covarB + cond.covar) * (cond.mean - meanB))
    covarA = covarA - covarAB * ((covarB + cond.covar) \ transpose(covarAB))
    return meanA, covarA
end