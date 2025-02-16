"""
    purity(state::GaussianState)

Calculate the purity of a Gaussian state, defined by `1/sqrt((2/ħ) det(V))`.
"""
purity(x::GaussianState) = (b = x.basis; (x.ħ/2)^(b.nmodes)/sqrt(det(x.covar)))