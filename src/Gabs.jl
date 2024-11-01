module Gabs

using BlockArrays: BlockedArray, BlockArray, Block, mortar

import LinearAlgebra
using LinearAlgebra: I, det, mul!

import QuantumInterface: StateVector, AbstractOperator, apply!, tensor, ⊗

export 
    # types
    GaussianState, GaussianUnitary, GaussianChannel, Generaldyne,
    # operations
    tensor, ⊗, apply!, ptrace, output,
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

include("utils.jl")

include("types.jl")

include("states.jl")

include("unitaries.jl")

include("channels.jl")

include("measurements.jl")

include("wigner.jl")

end
