using Documenter
using Dynema

DocMeta.setdocmeta!(Dynema, :DocTestSetup, :(using Dynema); recursive = true)

makedocs(
    sitename = "Dynema.jl",
    authors = "Jose Alquicira-Hernandez",
    modules = [Dynema],
    workdir = @__DIR__,  # so `include("src/assets/....jl")` in @example/@setup blocks
                         # resolves the same way regardless of which page runs it
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://joseah.github.io/Dynema.jl",
        edit_link = "main",
        assets = ["assets/favicon.ico"],
        sidebar_sitename = false,        # logo already includes the "Dynema" wordmark
        ansicolor = false,               # strip ANSI color instead of inlining <span> styles
        size_threshold = 600 * 2^10,      # 600 KiB hard limit (default 200 KiB)
        size_threshold_warn = 300 * 2^10, # 300 KiB warn limit (default 100 KiB)
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Quick start (end to end)" => "tutorials/quick_start.md",
            "Command-line overview" => "tutorials/command_line.md",
            "Extra: Compressing and indexing a VCF: `dynema-prepare-vcf`" => "tutorials/prepare_vcf.md",
            "Extra: Building gene bed files from a GTF: `dynema-prepare-bed`" => "tutorials/prepare_bed.md",
        ],
        "API Reference" => "functions.md",
        "Internals" => "internals.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/joseah/Dynema.jl.git",
    devbranch = "main",
    push_preview = true,
)


