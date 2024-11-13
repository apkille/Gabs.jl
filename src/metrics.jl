"""
    purity(state::GaussianState)

Calculate the purity of a Gaussian state, defined by `1/sqrt(det(V))`.
"""
purity(state::GaussianState) = 1/sqrt(det(state.covar))