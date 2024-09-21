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
        "Getting Started with Gabs" => "intro.md",
        "Manual" => "manual.md",
        "Tutorials" => "tutorials.md",
        "Gaussian Zoos" => "zoos.md",
        "API" => "API.md",
        "References" => "bibliography.md"
    ]
    )

    deploydocs(; repo = "github.com/apkille/Gabs.jl", devbranch = "master")
end

main()
