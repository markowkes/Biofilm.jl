using Documenter, Biofilm

makedocs(
    sitename="Biofilm Documentation",
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Postprocessing" => "postprocessing.md"
        #"Reference" => "reference.md",
        ]
    )

deploydocs(
    repo = "github.com/markowkes/Biofilm.jl.git",
    versions = nothing,
)