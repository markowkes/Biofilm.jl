using Documenter, Biofilm

makedocs(sitename="Biofilm Documentation")

deploydocs(
    repo = "github.com/markowkes/Biofilm.jl.git",
    versions = nothing,
)