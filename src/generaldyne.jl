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

function generaldyne(state::GaussianState{<:QuadBlockBasis,Tm,Tc}, indices::I, angle::N; result::S = nothing) where {Tm,Tc,I,N,S}
	indlength = length(indices)
	basis = state.basis
	notindices = setdiff(1:basis.nmodes, indices)
	nmodes′ = basis.nmodes - indlength
	mean, covar = state.mean, state.covar
	idxlength = length(indices)
	R = zeros(2*indlength, 2*indlength)
	D = Diagonal(repeat([1, 0], outer = indlength))
	A, B, C = zeros(2*nmodes′, 2*nmodes′), zeros(2*indlength, 2*indlength), zeros(2*nmodes′, 2*indlength)
	a, b = zeros(2*nmodes′), zeros(2*indlength)
	for i in eachindex(notindices)
		idx = notindices[i]
		a[2i-1:2i] .= @view(mean[2idx-1:2idx])
		for j in eachindex(notindices)
			otheridx = notindices[j]
			if idx == otheridx
				A[2i-1:2i, 2i-1:2i] .= @view(covar[2idx-1:2idx, 2idx-1:2idx])
			else
				A[2i-1:2i, 2j-1:2j] .= @view(covar[2idx-1:2idx, 2otheridx-1:2otheridx])
				A[2j-1:2j, 2i-1:2i] .= @view(covar[2otheridx-1:2otheridx, 2idx-1:2idx])
			end
		end
	end
	for i in eachindex(indices)
		idx = indices[i]
		R[2i-1, 2i-1] = cos(angle)
		R[2i-1, 2i] = -sin(angle)
		R[2i, 2i-1] = sin(angle)
		R[2i, 2i] = cos(angle)
		b[2i-1:2i] .= @view(mean[2idx-1:2idx])
		for j in eachindex(indices)
			otheridx = indices[j]
			if idx == otheridx
				B[2i-1:2i, 2i-1:2i] .= @view(covar[2idx-1:2idx, 2idx-1:2idx])
			else
				B[2i-1:2i, 2j-1:2j] .= @view(covar[2idx-1:2idx, 2otheridx-1:2otheridx])
				B[2j-1:2j, 2i-1:2i] .= @view(covar[2otheridx-1:2otheridx, 2idx-1:2idx])
			end
		end
		for j in eachindex(notindices)
			otheridx = notindices[j]
			C[2j-1:2j, 2i-1:2i] .= @view(covar[2otheridx-1:2otheridx, 2idx-1:2idx])
		end
	end
	copyto!(B, R * B * R')
	copyto!(C, C * R')
	copyto!(b, R * b)
	prod = Symmetric(D * B * D)
	mpprod = pinv(prod)
	buf = C * mpprod
	if isnothing(result)
		result = randn(MyNormal(b, prod))
		meandiff = bufvec
	else 
		meandiff = b - result
	end
	copyto!(a, a - buf * meandiff)
	copyto!(A, A - (buf * C'))
	state′ = GaussianState(QuadBlockBasis(nmodes′), a, A)
	return Homodyne(result, state′)
end