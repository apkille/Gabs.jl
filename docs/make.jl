using Revise # for interactive work on docs
push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterCitations
using Gabs
using QuantumInterface

DocMeta.setdocmeta!(Gabs, :DocTestSetup, :(using Gabs, QuantumInterface); recursive=true)

function main()
    bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"), style=:authoryear)

    makedocs(
    plugins=[bib],
    doctest = false,
    clean = true,
    sitename = "Gabs.jl",
    format = Documenter.HTML(
        assets=["assets/init.js"]
    ),
    modules = [Gabs, QuantumInterface],
    checkdocs = :exports,
    warnonly = false,
    authors = "Andrew Kille",
    pages = [
        "Gabs.jl" => "index.md",
        "API" => "API.md",
    ]
    )

    deploydocs(
        repo = "github.com/apkille/Gabs.jl.git"
    )
end

main()
