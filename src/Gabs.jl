module Gabs

import LinearAlgebra
using LinearAlgebra: I, det

import QuantumInterface: StateVector, AbstractOperator, apply!, directsum, ⊕

export 
    # types
    GaussianState, GaussianUnitary, GaussianChannel,
    # operations
    directsum, ⊕, apply, apply!, ptrace,
    # predefined Gaussian states
    vacuumstate, thermalstate, coherentstate, squeezedstate, eprstate,
    # predefined Gaussian channels
    displace, squeeze, twosqueeze, phaseshift, beamsplitter,
    attenuator, amplifier,
    # wigner functions
    wigner, wignerchar,
    # symplectic form
    symplecticform

include("errors.jl")

include("types.jl")

include("states.jl")

include("unitaries.jl")

include("channels.jl")

include("wigner.jl")

end
