module Gabs

import LinearAlgebra
using LinearAlgebra: I

import QuantumInterface: StateVector, AbstractOperator, apply!, directsum, ⊕

export 
    # types
    GaussianState, GaussianChannel,
    # operations
    directsum, ⊕, apply, apply!,
    # predefined Gaussian states
    vacuumstate, thermalstate, coherentstate, squeezedstate, eprstate,
    # predefined Gaussian channels
    displace, squeeze, twosqueeze, phaseshift, beamsplitter

include("errors.jl")

include("types.jl")

include("states.jl")

include("channels.jl")

end
