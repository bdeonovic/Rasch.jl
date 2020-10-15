using Rasch
using Documenter

makedocs(;
    modules=[Rasch],
    authors="Benjamin Deonovic",
    repo="https://github.com/bdeonovic/Rasch.jl/blob/{commit}{path}#L{line}",
    sitename="Rasch.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bdeonovic.github.io/Rasch.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bdeonovic/Rasch.jl",
)
