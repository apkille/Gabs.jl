module Gabs

import LinearAlgebra
using LinearAlgebra: I, det, mul!, diag, qr, eigvals, Diagonal, cholesky, Symmetric

import QuantumInterface: StateVector, AbstractOperator, apply!, tensor, ⊗, directsum, ⊕, entropy_vn, fidelity, logarithmic_negativity

import Random
using Random: randn!

import SymplecticFactorizations: williamson, Williamson, polar, Polar, blochmessiah, BlochMessiah, randsymplectic, symplecticform, issymplectic
using SymplecticFactorizations: williamson, Williamson, polar, Polar, blochmessiah, BlochMessiah, BlockForm, PairForm

export
    # types
    GaussianState, GaussianUnitary, GaussianChannel,
    # Gaussian measurements
    generaldyne, Generaldyne,
    # symplectic representations
    QuadPairBasis, QuadBlockBasis, changebasis,
    # operations
    tensor, ⊗, directsum, ⊕, apply!, ptrace,
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
    williamson, Williamson, polar, Polar, blochmessiah, BlochMessiah,
    # metrics
    purity, entropy_vn, fidelity, logarithmic_negativity

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

include("generaldyne.jl")

include("wigner.jl")

include("metrics.jl")

end
