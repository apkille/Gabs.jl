module MakieExt

using Makie
using Gabs: GaussianState, wigner, wignerchar

Makie.used_attributes(::Type{<:Makie.Heatmap}, ::Any, ::Any, x::GaussianState) = (:dist,)

# Directly use heatmap to create wigner distribution
function Makie.convert_arguments(P::Type{<:Makie.Heatmap}, q, p, state::GaussianState; dist = :wigner)
    isequal(length(state.mean),2) || throw(ArgumentError(Gabs.HEATMAP_ERROR))
    if isequal(dist, :wigner)
        data = [wigner(state, [i,j]) for i in q, j in p]
    elseif isequal(dist, :wignerchar)
        data = [real(wignerchar(state, [i,j])) for i in q, j in p]
    else
        throw(ErrorException("`distribution` should be `:wigner` or `:wignerchar`"))
    end
    return convert_arguments(P, q, p, data)
end

end