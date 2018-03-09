using Documenter, DFControl

makedocs(modules = [DFControl],
         pages = Any["Home" => "index.md",
                     "job.md",
                     "input.md",
                     "structure.md",
                     "atom.md",
                     "fileio.md",
                     "defaults.md",
                     "utils.md"],
         sitename = "DFControl.jl",
         format = :html
        )
