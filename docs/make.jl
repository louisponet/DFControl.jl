using Documenter, DFControl

makedocs(modules = [DFControl],
         pages = Any["Home" => "index.md",
                     "Library" => Any["job.md",
                     "input.md",
                     "structure.md",
                     "atom.md",
                     "fileio.md",
                     "defaults.md",
                     "utils.md"]],
         sitename = "DFControl.jl",
         format = :html
        )


deploydocs(
    repo = "github.com/louisponet/DFControl.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps = nothing,
    make = nothing
    )
