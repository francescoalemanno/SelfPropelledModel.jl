using Documenter, SelfPropelledModel

makedocs(
    modules = [SelfPropelledModel],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Francesco Alemanno",
    sitename = "SelfPropelledModel.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/francescoalemanno/SelfPropelledModel.jl.git",
    push_preview = true
)
