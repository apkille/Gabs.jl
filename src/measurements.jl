abstract type AbstractGeneraldyne end

"""Defines a partial general-dyne measurement of a Gaussian system.

## Fields

- `system`: Initial Gaussian state.
- `outcome`: Measurement outcome Gaussian state.
- `indices`: Specific indices of `system` that are projected onto `outcome`.

## Mathematical description of a general-dyne measurement

A partial general-dyne measurement on a Gaussian state ``\\hat{\\rho}(\\mathbf{\\bar{x}}, \\mathbf{V})`` of a system partitioned
into subsystems ``A`` and ``B`` is a projection of subsystem ``B`` onto
a measurement outcome Gaussian state ``\\hat{\\rho}(\\mathbf{\\bar{x}}_{m}, \\mathbf{V}_{m})``.

## Example

```jldoctest
julia> basis = QuadPairBasis(1);

julia> vac = vacuumstate(basis); coh = coherentstate(basis, 1.0-im);

julia> state = vac ⊗ coh ⊗ vac ⊗ coh;

julia> Generaldyne(state, coh ⊗ vac ⊗ coh, [1, 3, 4])
Generaldyne on indices [1, 3, 4]
system: GaussianState for 4 modes.
 (xType: Vector{Float64} | vType: Matrix{Float64})
outcome: GaussianState for 3 modes.
 (xType: Vector{Float64} | vType: Matrix{Float64})
```
"""
struct Generaldyne{I} <: AbstractGeneraldyne
    system::GaussianState
    outcome::GaussianState
    indices::I
    function Generaldyne(sys::GaussianState, outcome::GaussianState, ind::I) where {I}
        sysbasis, outbasis = sys.basis, outcome.basis
        typeof(sysbasis) == typeof(outbasis) || throw(ArgumentError(SYMPLECTIC_ERROR))
        length(ind) == outbasis.nmodes || throw(DimensionMismatch(GENERALDYNE_ERROR))
        return new{I}(sys, outcome, ind)
    end
end

"""
    output(meas::Generaldyne)

Conditional mapping of the subsystem for the initial Gaussian system that is not measured 
in a `Generaldyne` object.

## Mathematical description of conditional evolution of a Gaussian system

Let ``\\hat{\\rho}(\\mathbf{\\bar{x}}, \\mathbf{V})`` be a Gaussian state of a system partitioned
into subsystems ``A`` and ``B``, such that
```math
\\mathbf{\\bar{x}} = \\begin{pmatrix} 
                        \\mathbf{\\bar{x}}_{\\text{A}} & \\mathbf{\\bar{x}}_{\\text{A}} 
                    \\end{pmatrix}^{\\text{T}},
\\quad \\mathbf{V} = \\begin{pmatrix} 
                        \\mathbf{V}_{\\text{A}} & \\mathbf{V}_{\\text{AB}} \\\\
                        \\mathbf{V}_{\\text{AB}}^{\\text{T}} & \\mathbf{V}_{\\text{B}}
                    \\end{pmatrix},
```
where ``(\\mathbf{\\bar{x}}_{\\text{A}}, \\mathbf{V}_{\\text{A}})`` and 
``(\\mathbf{\\bar{x}}_{\\text{B}}, \\mathbf{V}_{\\text{B}})`` are the mean and covariance matrices
of subsystems ``A`` and ``B``, respectively. Here ``\\mathbf{V}_{\\text{AB}}`` is a covariance matrix
representing the correlations between ``A`` and ``B``. If we project subsystem ``B`` onto
a measurement outcome ``\\hat{\\rho}(\\mathbf{\\bar{x}}_{m}, \\mathbf{V}_{m})``, then
subsystem ``A`` conditionally maps as follows:
```math
\\mathbf{\\bar{x}}_{\\text{A}} \\to \\mathbf{\\bar{x}}_{\\text{A}} + 
    \\frac{\\mathbf{V}_{\\text{AB}}}{\\mathbf{V}_{\\text{B}} + \\mathbf{V}_{m}}(\\mathbf{\\bar{x}}_{m} - \\mathbf{\\bar{x}}_{\\text{B}})^{\\text{T}},
\\quad\\mathbf{V}_{\\text{A}} \\to \\mathbf{V}_{\\text{A}} - 
    \\frac{\\mathbf{V}_{\\text{AB}}}{\\mathbf{V}_{\\text{B}} + \\mathbf{V}_{m}} \\mathbf{V}_{\\text{AB}}^{\\text{T}}.
```

## Example

```jldoctest
julia> basis = QuadPairBasis(1);

julia> vac = vacuumstate(basis); coh = coherentstate(basis, 1.0-im);

julia> state = vac ⊗ coh ⊗ vac ⊗ coh;

julia> gd = Generaldyne(state, coh ⊗ vac ⊗ coh, [1, 3, 4])
Generaldyne on indices [1, 3, 4]
system: GaussianState for 4 modes.
 (xType: Vector{Float64} | vType: Matrix{Float64})
outcome: GaussianState for 3 modes.
 (xType: Vector{Float64} | vType: Matrix{Float64})

julia> output(gd)
GaussianState for 1 mode.
  symplectic basis: QuadPairBasis
mean: 2-element Vector{Float64}:
  2.0
 -2.0
covariance: 2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
"""
function output(meas::Generaldyne)
    sys, outcome, ind = meas.system, meas.outcome, meas.indices
    sysbasis, outbasis = sys.basis, outcome.basis
    # partition mean and covar into subsystems A and B
    meanA, meanB = _part_mean(sys, ind)
    covarA, covarB, covarAB = _part_covar(sys, outcome, ind)
    # Block array matmul and broadcasting is incredibly
    # slow, so convert types back to promoted `sys`` and `cond`` types
    meanA′, meanB′ = sys.mean isa Vector{Float64} ? Vector{Float64}.((meanA, meanB)) :
        (_promote_output_vector(typeof(outcome.mean), meanA, 2*(sysbasis.nmodes - outbasis.nmodes)), _promote_output_vector(typeof(outcome.mean), meanB, 2*outbasis.nmodes))
    covarA′, covarB′, covarAB′ = sys.covar isa Matrix{Float64} ? Matrix{Float64}.((covarA, covarB, covarAB)) :
        (_promote_output_matrix(typeof(outcome.covar), covarA, (2*(sysbasis.nmodes - outbasis.nmodes), 2*(sysbasis.nmodes - outbasis.nmodes))),
        _promote_output_matrix(typeof(outcome.covar), covarB, (2*outbasis.nmodes, 2*outbasis.nmodes)),
        _promote_output_matrix(typeof(outcome.covar), covarAB, (2*(sysbasis.nmodes - outbasis.nmodes), 2*outbasis.nmodes)))
    # map subsystem A
    meanA′, covarA′ = _generaldyne_map(meanA′, meanB′, covarA′, covarB′, covarAB′, sys, outcome)
    return GaussianState(typeof(outbasis)(sysbasis.nmodes-outbasis.nmodes) , meanA′, covarA′)
end

"""
    prob(meas::Generaldyne)

Calculate the probability density of a general-dyne measurement.

## Mathematical description of probability density of a general-dyne measurement

Let ``\\hat{\\rho}(\\mathbf{\\bar{x}}, \\mathbf{V})`` be a Gaussian state of a system partitioned
into subsystems ``A`` and ``B``, such that
```math
\\mathbf{\\bar{x}} = \\begin{pmatrix} 
                        \\mathbf{\\bar{x}}_{\\text{A}} & \\mathbf{\\bar{x}}_{\\text{A}} 
                    \\end{pmatrix}^{\\text{T}},
\\quad \\mathbf{V} = \\begin{pmatrix} 
                        \\mathbf{V}_{\\text{A}} & \\mathbf{V}_{\\text{AB}} \\\\
                        \\mathbf{V}_{\\text{AB}}^{\\text{T}} & \\mathbf{V}_{\\text{B}}
                    \\end{pmatrix},
```
where ``(\\mathbf{\\bar{x}}_{\\text{A}}, \\mathbf{V}_{\\text{A}})`` and 
``(\\mathbf{\\bar{x}}_{\\text{B}}, \\mathbf{V}_{\\text{B}})`` are the mean and covariance matrices
of subsystems ``A`` and ``B``, respectively. Here ``\\mathbf{V}_{\\text{AB}}`` is a covariance matrix
representing the correlations between ``A`` and ``B``. If we project subsystem ``B`` onto
a measurement outcome ``\\hat{\\rho}(\\mathbf{\\bar{x}}_{m}, \\mathbf{V}_{m})``, then
the probability density in ``d\\mathbf{\\hat{x}}_{m}`` is

```math
p(\\mathbf{\\bar{x}}_{m}) = \\frac{e^{(\\mathbf{\\bar{x}}_{m} - \\mathbf{\\bar{x}}_{\\text{B}})^{\\text{T}} 
    \\frac{1}{\\mathbf{V}_{\\text{B}} + \\mathbf{V}_{m}} (\\mathbf{\\bar{x}}_{m} - \\mathbf{\\bar{x}}_{\\text{B}})}}{\\pi^{m} 
    \\sqrt{\\det(\\mathbf{V}_{\\text{B}} + \\mathbf{V}_{m})}},
```
where ``m`` is the number of modes in subsystem ``B``.

## Example

```jldoctest
julia> basis = QuadPairBasis(1);

julia> vac = vacuumstate(basis); coh = coherentstate(basis, 1.0-im);

julia> state = vac ⊗ coh ⊗ vac ⊗ coh;

julia> gd = Generaldyne(state, coh ⊗ vac ⊗ coh, [1, 3, 4])
Generaldyne on indices [1, 3, 4]
system: GaussianState for 4 modes.
 (xType: Vector{Float64} | vType: Matrix{Float64})
outcome: GaussianState for 3 modes.
 (xType: Vector{Float64} | vType: Matrix{Float64})

julia> prob(gd)
0.2201092644728679
```
"""
function prob(meas::Generaldyne)
    sys, outcome, ind = meas.system, meas.outcome, meas.indices
    sysbasis, outbasis = sys.basis, outcome.basis
    # partition mean and covar into subsystems
    mean′ = BlockedArray(sys.mean, 2*ones(Int,sysbasis.nmodes))
    meanB = _part_meanB(mean′, ind)
    sizeB = (outbasis.nmodes, outbasis.nmodes)
    covar′ = BlockedArray(sys.covar, 2*ones(Int,sysbasis.nmodes), 2*ones(Int,sysbasis.nmodes))
    covarB = _part_covarB(covar′, sysbasis.nmodes, ind, sizeB)
    # Block array matmul and broadcasting is incredibly
    # slow, so convert types back to promoted `sys`` and `cond`` types
    meanB′ = sys.mean isa Vector{Float64} ? Vector{Float64}(meanB) :
        _promote_output_vector(typeof(outcome.mean), meanB, 2*outbasis.nmodes)
    covarB′ = sys.covar isa Matrix{Float64} ? Matrix{Float64}(covarB) :
        _promote_output_matrix(typeof(outcome.covar), covarB, sizeB)
    return _prob_formula(meanB′, covarB′, outcome)
end

function Base.show(io::IO, mime::MIME"text/plain", x::Generaldyne)
    Base.summary(io, x)
    sys, outcome = x.system, x.outcome
    print(io, "\nsystem: ")
    Base.summary(io, sys)
    print(io, "\n (xType: $(typeof(sys.mean)) | vType: $(typeof(sys.covar)))") 
    print(io, "\noutcome: ")
    Base.summary(io, x.outcome)
    print(io, "\n (xType: $(typeof(outcome.mean)) | vType: $(typeof(outcome.covar)))")
end
function Base.summary(io::IO, x::Generaldyne)
    printstyled(io, :Generaldyne; color=:light_cyan)
    print(io, " on indices $(x.indices)")
end

function _prob_formula(mean, covar, outcome)
    outbasis = outcome.basis
    # create alloc buffers for matrix multiplication
    buf = zeros(2*outbasis.nmodes)
    meandiff = outcome.mean .- mean
    norm = pi^(outbasis.nmodes)*sqrt(det(covar .+ outcome.covar))
    return exp(transpose(meandiff) * mul!(buf, inv(covar .+ outcome.covar), meandiff))/norm
end
function _generaldyne_map(meanA, meanB, covarA, covarB, covarAB, sys, outcome)
    sysbasis, outbasis = sys.basis, outcome.basis
    # create alloc buffers for matrix multiplication
    meanbuf1, meanbuf2 = zeros(2*outbasis.nmodes), zeros(2*(sysbasis.nmodes - outbasis.nmodes))
    covarbuf = zeros(2*(sysbasis.nmodes - outbasis.nmodes), 2*(sysbasis.nmodes - outbasis.nmodes))
    # maps subsystem A, which is not measured
    meanA .= meanA .+ mul!(meanbuf2, covarAB, (mul!(meanbuf1, inv(covarB .+ outcome.covar), outcome.mean .- meanB)))
    covarA .= covarA .- mul!(covarbuf, covarAB, (covarB .+ outcome.covar) \ transpose(covarAB))
    return meanA, covarA
end
function _part_mean(sys::M, ind::I) where {M,I}
    sysbasis = sys.basis
    # block mean into its modes
    mean′ = BlockedArray(sys.mean, 2*ones(Int,sysbasis.nmodes))
    meanA = _part_meanA(mean′, sysbasis.nmodes, ind)
    meanB = _part_meanB(mean′, ind)
    return meanA, meanB
end
_part_meanA(mean::M, nmodes::N, ind::I) where {M,N,I} = mortar([view(mean, Block(i)) for i in Base.OneTo(nmodes) if !(i in ind)])
_part_meanB(mean::M, ind::I) where {M,I} = mortar([view(mean, Block(i)) for i in ind])
function _part_covar(sys::M, outcome::C, ind::I) where {M,C,I}
    sysbasis, outbasis = sys.basis, outcome.basis
    sizeA, sizeB = (sysbasis.nmodes-outbasis.nmodes, sysbasis.nmodes-outbasis.nmodes), (outbasis.nmodes, outbasis.nmodes)
    sizeAB = (sysbasis.nmodes-outbasis.nmodes, outbasis.nmodes)
    # loop through entire Gaussian system, writing quadratures to B if
    # index is specified in `ind` argument
    covar′ = BlockedArray(sys.covar, 2*ones(Int,sysbasis.nmodes), 2*ones(Int,sysbasis.nmodes))
    covarA = _part_covarA(covar′, sysbasis.nmodes, ind, sizeA)
    covarB = _part_covarB(covar′, sysbasis.nmodes, ind, sizeB)
    covarAB = _part_covarAB(covar′, sysbasis.nmodes, ind, sizeAB)
    return covarA, covarB, covarAB
end
function _part_covarA(covar::C, nmodes::N, ind::I, size::S) where {C,N,I,S}
    return mortar(reshape([view(covar, Block(i,j)) for i in Base.OneTo(nmodes), j in Base.OneTo(nmodes)
                if !(i in ind) && !(j in ind)], size))
end
function _part_covarB(covar::C, nmodes::N, ind::I, size::S) where {C,N,I,S}
    return mortar(reshape([view(covar, Block(i,j)) for i in Base.OneTo(nmodes), j in Base.OneTo(nmodes)
        if i in ind && j in ind], size))
end
function _part_covarAB(covar::C, nmodes::N, ind::I, size::S) where {C,N,I,S}
    return mortar(reshape([view(covar, Block(j, i)) for i in Base.OneTo(nmodes) for j in Base.OneTo(i) 
        if (i in ind && !(j in ind)) || (!(i in ind) && j in ind)], size))
end