module Gabs

using BlockArrays: BlockedArray, BlockArray, Block, mortar

import LinearAlgebra
using LinearAlgebra: I, det, mul!, diagm, diag, qr, eigvals

import QuantumInterface: StateVector, AbstractOperator, apply!, tensor, ⊗, directsum, ⊕

import SymplecticFactorizations: williamson, Williamson, polar, Polar, randsymplectic, symplecticform, issymplectic
using SymplecticFactorizations: williamson, Williamson, polar, Polar, BlockForm, PairForm

export 
    # types
    GaussianState, GaussianUnitary, GaussianChannel, Generaldyne,
    # symplectic representations
    QuadPairBasis, QuadBlockBasis, changebasis,
    # operations
    tensor, ⊗, directsum, ⊕, apply!, ptrace, output, prob,
    # predefined Gaussian states
    vacuumstate, thermalstate, coherentstate, squeezedstate, eprstate,
    # predefined Gaussian channels
    displace, squeeze, twosqueeze, phaseshift, beamsplitter,
    attenuator, amplifier,
    # random objects
    randstate, randunitary, randchannel, randsymplectic,
    # wigner functions
    wigner, wignerchar,
    # symplectic form and checks
    symplecticform, issymplectic, isgaussian, sympspectrum,
    # factorizations
    williamson, Williamson, polar, Polar,
    # metrics
    purity

include("errors.jl")

include("utils.jl")

include("symplectic.jl")

include("types.jl")

include("states.jl")

include("unitaries.jl")

include("channels.jl")

include("randoms.jl")

include("factorizations.jl")

include("measurements.jl")

include("wigner.jl")

include("metrics.jl")

end
