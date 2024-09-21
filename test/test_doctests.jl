@testitem "Doctests" tags=[:doctests] begin
    using Documenter
    using Gabs

    DocMeta.setdocmeta!(Gabs, :DocTestSetup, :(using Gabs); recursive=true)
    doctest(Gabs)
end