abstract type AbstractGaussianMeasurement end

struct Generaldyne{R,S<:GaussianState} <: AbstractGaussianMeasurement
	result::R
	state::S
	function Generaldyne(r::R, s::S) where {R,S<:GaussianState}
		return new{R,S}(r, s)
	end
end

# iteration for destructuring into components
Base.iterate(F::Generaldyne) = (F.result, Val(:state))
Base.iterate(F::Generaldyne, ::Val{:state}) = (F.state, Val(:done))
Base.iterate(F::Generaldyne, ::Val{:done}) = nothing

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, H::Generaldyne{<:Any,<:GaussianState})
    Base.summary(io, H); println(io)
    println(io, "result:")
    Base.show(io, mime, H.result)
    println(io, "\noutput state:")
    Base.show(io, mime, H.state)
end

"""
    generaldyne(state::GaussianState, indices::Vector; select) -> Generaldyne
    generaldyne(state::GaussianState, index::Int; select) -> Generaldyne

Compute the projection of the subsystem of a Gaussian state `state` indicated by `indices`
on `select` and return a `Generaldyne` object. The keyword argument `select` can take the following forms:

- If `select` is a matrix, then the subsystem is projected onto a Gaussian state with a randomly sampled mean and covariance matrix `result`.
- If `select` is a Gaussian state, then the subsystem is projected onto `select`.

The `result` and mapped state `output` can be obtained from the Generaldyne object `M` via `M.result` and `M.output`.
Iterating the decomposition produces the components `result` and `output`.

# Examples
```jldoctest
julia> st = squeezedstate(QuadPairBasis(3), 1.0, pi/4);

julia> M = generaldyne(st, [1, 3], Matrix{Float64}(I, 4, 4))
Generaldyne{Vector{Float64}, GaussianState{QuadBlockBasis{Int64}, Vector{Float64}, Matrix{Float64}}}
result:
4-element Vector{Float64}:
  0.7763070778411039
 -4.333442961163451
 -0.6052604156263252
  1.041771385283006
output state:
GaussianState for 1 mode.
  symplectic basis: QuadBlockBasis
mean: 2-element Vector{Float64}:
 0.0
 0.0
covariance: 2×2 Matrix{Float64}:
  1.19762  -2.56458
 -2.56458   6.32677

julia> result, output = M; # destructuring via iteration

julia> result == M.result && output == M.output
true
```
"""
function generaldyne(state::GaussianState, indices::I; select::S) where {I,S<:Union{Matrix,GaussianState}}
	basis = state.basis
	indlength = length(indices)
	nmodes′ = basis.nmodes - indlength
	a, b, A, B, C = _part_state(state, indices)
	if select isa Matrix
		B .+= select
		symB = Symmetric(B)
		L = cholesky(symB).L
		resultmean = L * randn(2*indlength) + b
		meandiff = resultmean - b
		buf = C * inv(symB)
		a .+= buf * meandiff
		A .-= buf * C'
		result′ = GaussianState(typeof(basis)(indlength), resultmean, select, ħ = state.ħ)
	elseif select isa GaussianState
		B .+= select.covar
		symB = Symmetric(B)
		meandiff = select.mean - b
		buf = C * inv(symB)
		a .+= buf * meandiff
		A .-= buf * C'
		result′ = select
	else
		throw(ArgumentError("invalid `select`"))
	end
	state′ = GaussianState(typeof(basis)(nmodes′), a, A, ħ = state.ħ)
	return Generaldyne(result′, state′)
end

function Base.rand(::Type{Generaldyne}, state::GaussianState{<:QuadPairBasis,Tm,Tc}, indices::I; shots::Int = 1, select::Matrix) where {Tm,Tc,I}
	indlength = length(indices)
	basis = state.basis
	nmodes′ = basis.nmodes - indlength
	mean, covar = state.mean, state.covar
	idxlength = length(indices)
	b, B = zeros(2*indlength), zeros(2*indlength, 2*indlength)
	@inbounds for i in eachindex(indices)
		idx = indices[i]
		b[2i-1:2i] .= @view(mean[2idx-1:2idx])
		@inbounds for j in eachindex(indices)
			otheridx = indices[j]
			if idx == otheridx
				B[2i-1:2i, 2i-1:2i] .= @view(covar[2idx-1:2idx, 2idx-1:2idx])
			else
				B[2i-1:2i, 2j-1:2j] .= @view(covar[2idx-1:2idx, 2otheridx-1:2otheridx])
				B[2j-1:2j, 2i-1:2i] .= @view(covar[2otheridx-1:2otheridx, 2idx-1:2idx])
			end
		end
	end
	symB = Symmetric(B)
	L = cholesky(symB).L
	buf = zeros(2*indlength)
	results = zeros(shots, 2*indlength)
	@inbounds for i in Base.OneTo(shots)
		mul!(@view(results[i,:]), L, randn!(buf))
		@view(results[i, :]) .+= b
	end
	return results
end
function Base.rand(::Type{Generaldyne}, state::GaussianState{<:QuadBlockBasis,Tm,Tc}, indices::I; shots::Int = 1, select::Matrix) where {Tm,Tc,I}
	indlength = length(indices)
	basis = state.basis
	nmodes = basis.nmodes
	nmodes′ = nmodes - indlength
	mean, covar = state.mean, state.covar
	idxlength = length(indices)
	b, B = zeros(2*indlength), zeros(2*indlength, 2*indlength)
	@inbounds for i in eachindex(indices)
		idx = indices[i]
		b[i] = mean[idx]
		b[i+nmodes′] = mean[idx+nmodes]
		@inbounds for j in eachindex(indices)
			otheridx = indices[j]
			if idx == otheridx
				B[i, i] = covar[idx, idx]
				B[i+indlength, i] = covar[idx+nmodes, idx]
				B[i, i+indlength] = covar[idx, idx+nmodes]
				B[i+indlength, i+indlength] = covar[idx+nmodes, idx+nmodes]
			else
				B[i, j] = covar[idx, otheridx]
				B[i+indlength, j] = covar[idx+nmodes, otheridx]
				B[i, j+indlength] = covar[idx, otheridx+nmodes]
				B[i+indlength, j+indlength] = covar[idx+nmodes, otheridx+nmodes]

				B[j, i] = covar[otheridx, idx]
				B[j+indlength, i] = covar[otheridx+nmodes, idx]
				B[j, i+indlength] = covar[otheridx, idx+nmodes]
				B[j+indlength, i+indlength] = covar[otheridx+nmodes, idx+nmodes]
			end
		end
	end
	B .+= select
	symB = Symmetric(B)
	L = cholesky(symB).L
	buf = zeros(2*indlength)
	results = zeros(shots, 2*indlength)
	@inbounds for i in Base.OneTo(shots)
		mul!(@view(results[i,:]), L, randn!(buf))
		@view(results[i, :]) .+= b
	end
	return results
end

function _part_state(state::GaussianState{<:QuadPairBasis,M,V}, indices::I) where {M,V,I}
	indlength = length(indices)
	basis = state.basis
	notindices = setdiff(1:basis.nmodes, indices)
	nmodes = basis.nmodes
	nmodes′ = nmodes - indlength
	mean, covar = state.mean, state.covar
	idxlength = length(indices)
	# partition state into subsystems A and B
	A, B, C = zeros(2*nmodes′, 2*nmodes′), zeros(2*indlength, 2*indlength), zeros(2*nmodes′, 2*indlength)
	a, b = zeros(2*nmodes′), zeros(2*indlength)
	@inbounds for i in eachindex(notindices)
		idx = notindices[i]
		a[2i-1:2i] .= @view(mean[2idx-1:2idx])
		@inbounds for j in eachindex(notindices)
			otheridx = notindices[j]
			if idx == otheridx
				A[2i-1:2i, 2i-1:2i] .= @view(covar[2idx-1:2idx, 2idx-1:2idx])
			else
				A[2i-1:2i, 2j-1:2j] .= @view(covar[2idx-1:2idx, 2otheridx-1:2otheridx])
				A[2j-1:2j, 2i-1:2i] .= @view(covar[2otheridx-1:2otheridx, 2idx-1:2idx])
			end
		end
	end
	@inbounds for i in eachindex(indices)
		idx = indices[i]
		b[2i-1:2i] .= @view(mean[2idx-1:2idx])
		@inbounds for j in eachindex(indices)
			otheridx = indices[j]
			if idx == otheridx
				B[2i-1:2i, 2i-1:2i] .= @view(covar[2idx-1:2idx, 2idx-1:2idx])
			else
				B[2i-1:2i, 2j-1:2j] .= @view(covar[2idx-1:2idx, 2otheridx-1:2otheridx])
				B[2j-1:2j, 2i-1:2i] .= @view(covar[2otheridx-1:2otheridx, 2idx-1:2idx])
			end
		end
		@inbounds for j in eachindex(notindices)
			otheridx = notindices[j]
			C[2j-1:2j, 2i-1:2i] .= @view(covar[2otheridx-1:2otheridx, 2idx-1:2idx])
		end
	end
	return a, b, A, B, C
end
function _part_state(state::GaussianState{<:QuadBlockBasis,M,V}, indices::I) where {M,V,I}
	indlength = length(indices)
	basis = state.basis
	notindices = setdiff(1:basis.nmodes, indices)
	nmodes = basis.nmodes
	nmodes′ = nmodes - indlength
	mean, covar = state.mean, state.covar
	idxlength = length(indices)
	# partition state into subsystems A and B
	A, B, C = zeros(2*nmodes′, 2*nmodes′), zeros(2*indlength, 2*indlength), zeros(2*nmodes′, 2*indlength)
	a, b = zeros(2*nmodes′), zeros(2*indlength)
	@inbounds for i in eachindex(notindices)
		idx = notindices[i]
		a[i] = mean[idx]
		a[i+nmodes′] = mean[idx+nmodes]
		@inbounds for j in eachindex(notindices)
			otheridx = notindices[j]
			if idx == otheridx
				A[i, i] = covar[idx, idx]
				A[i+nmodes′, i] = covar[idx+nmodes, idx]
				A[i, i+nmodes′] = covar[idx, idx+nmodes]
				A[i+nmodes′, i+nmodes′] = covar[idx+nmodes, idx+nmodes]
			else
				A[i, j] = covar[idx, otheridx]
				A[i+nmodes′, j] = covar[idx+nmodes, otheridx]
				A[i, j+nmodes′] = covar[idx, otheridx+nmodes]
				A[i+nmodes′, j+nmodes′] = covar[idx+nmodes, otheridx+nmodes]

				A[j, i] = covar[otheridx, idx]
				A[j+nmodes′, i] = covar[otheridx+nmodes, idx]
				A[j, i+nmodes′] = covar[otheridx, idx+nmodes]
				A[j+nmodes′, i+nmodes′] = covar[otheridx+nmodes, idx+nmodes]
			end
		end
	end
	@inbounds for i in eachindex(indices)
		idx = indices[i]
		b[i] = mean[idx]
		b[i+nmodes′] = mean[idx+nmodes]
		@inbounds for j in eachindex(indices)
			otheridx = indices[j]
			if idx == otheridx
				B[i, i] = covar[idx, idx]
				B[i+indlength, i] = covar[idx+nmodes, idx]
				B[i, i+indlength] = covar[idx, idx+nmodes]
				B[i+indlength, i+indlength] = covar[idx+nmodes, idx+nmodes]
			else
				B[i, j] = covar[idx, otheridx]
				B[i+indlength, j] = covar[idx+nmodes, otheridx]
				B[i, j+indlength] = covar[idx, otheridx+nmodes]
				B[i+indlength, j+indlength] = covar[idx+nmodes, otheridx+nmodes]

				B[j, i] = covar[otheridx, idx]
				B[j+indlength, i] = covar[otheridx+nmodes, idx]
				B[j, i+indlength] = covar[otheridx, idx+nmodes]
				B[j+indlength, i+indlength] = covar[otheridx+nmodes, idx+nmodes]
			end
		end
		@inbounds for j in eachindex(notindices)
			otheridx = notindices[j]
			C[j, i] = covar[otheridx, idx]
			C[j+nmodes′, i] = covar[otheridx+nmodes, idx]
			C[j, i+indlength] = covar[otheridx, idx+nmodes]
			C[j+nmodes′, i+indlength] = covar[otheridx+nmodes, idx+nmodes]
		end
	end
	return a, b, A, B, C
end