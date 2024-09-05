module MakieExt

using Makie
using Gabs

function Makie.convert_arguments(P::Type{<:Makie.Heatmap}, q, p, s::GaussianState)
    wigner_data = [wigner(s, [i,j]) for i in q, j in p]
    Makie.convert_arguments(P, q, p, wigner_data)
end

end