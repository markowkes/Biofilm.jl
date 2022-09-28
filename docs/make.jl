using Documenter, Biofilm

makedocs(
    sitename="Biofilm Documentation",
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Case Parameters" => "parameters.md",
        "Examples" => "examples.md",
        "Postprocessing" => "postprocessing.md",
        "Solver Theory" => "theory.md",
        ]
    )

deploydocs(
    repo = "github.com/markowkes/Biofilm.jl.git",
    versions = nothing,
)