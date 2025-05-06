abstract type AbstractGaussianMeasurement end

"""
    _part_state(state::GaussianState{<:QuadPairBasis,M,V}, indices::Vector) -> a, b, A, B, C
    _part_state(state::GaussianState{<:QuadPairBasis,M,V}, indices::Int) -> a, b, A, B, C

Low-level function that partitions `state` into subsystems A and B, 
the latter system's modes specified by `indices`. The vectors `a` and `b`
are mean vectors of systems A and B, respectively. The matrices `A` and `B`
are covariance matrices of systems A and B, respectively. Matrix `C` is the 
correlation matrix between A and B.
"""
function _part_state(state::GaussianState{<:QuadPairBasis,M,V}, indices::I) where {M,V,I}
	indlength = length(indices)
	basis = state.basis
	notindices = setdiff(1:basis.nmodes, indices)
	nmodes = basis.nmodes
	nmodes′ = nmodes - indlength
	mean, covar = state.mean, state.covar
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
		b[i+indlength] = mean[idx+nmodes]
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