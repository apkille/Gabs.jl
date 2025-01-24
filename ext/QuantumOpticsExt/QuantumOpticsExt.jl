module QuantumOpticsExt

using QuantumOpticsBase
using Gabs: GaussianState, wigner, wignerchar

import QuantumInterface: express, QuantumOpticsRepr

function express(state::GaussianState, repr::QuantumOpticsRepr)
    basis = state.basis
    N, C = basis.nmodes, repr.cutoff
    fb = FockBasis(C^N)
    op = Operator(fb, zeros(C^N, C^N))
    spec = sympspectrum(state)
    @inbounds for k in Base.OneTo(N)
        
    end
end

end