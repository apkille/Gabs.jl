"""
    purity(state::GaussianState)

Calculate the purity of a Gaussian state, defined by `1/sqrt((2/ħ) det(V))`.
"""
purity(x::GaussianState) = (b = x.basis; (x.ħ/2)^(b.nmodes)/sqrt(det(x.covar)))

"""
    purity(state::PDGaussianState)

Calculate the purity of a Gaussian state using optimized PDMats operations.
"""
function purity(x::PDGaussianState)
    b = x.basis
    (x.ħ/2)^(b.nmodes)/sqrt(det(x.covar))
end