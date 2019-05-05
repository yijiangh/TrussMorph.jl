using Documenter, TrussMorph

makedocs(;
    modules=[TrussMorph],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/yijiangh/TrussMorph.jl/blob/{commit}{path}#L{line}",
    sitename="TrussMorph.jl",
    authors="Yijiang Huang",
    assets=[],
)

deploydocs(;
    repo="github.com/yijiangh/TrussMorph.jl",
)
