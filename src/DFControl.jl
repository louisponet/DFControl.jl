module DFControl
    using LinearAlgebra
    using Statistics
    using Media
    import Base.Iterators.flatten
    using RecipesBase
    include("FixedSizeArrays.jl")
    using .FixedSizeArrays
    using Parameters

    abstract type Package end
    struct Wannier90 <: Package end
    struct QE <: Package end
    struct Abinit <: Package end

    const depsdir = joinpath(dirname(@__DIR__), "deps")

    include("atom.jl")
    include("structure.jl")
    include("types.jl")

    include("input.jl")
    include("utils.jl")
    include("job.jl")

    export DFJob, Exec, DFInput
    include("API.jl")

    include("constants.jl")

    include("fileio.jl")
    include("plotting.jl")

    include("server.jl")
    export qstat

    include("defaults.jl")
    export setdefault_pseudodir
    export setdefault_server
    export configuredefault_pseudos
    export removedefault_pseudodir
    export removedefault_pseudos
    export setdefault_jobheader
    export @add_default
    export setdefault_input
    export removedefault_input
    export getdefault_pseudo
    export getdefault_server
    export getdefault_jobheader
    export getdefault_pseudodir
    export getdefault_pseudodirs
    export removedefault

    using Crayons
    const dfprintln = println
    const dfprint = print
    include("display/overloads.jl")
    using Requires
    function __init__()
        @require Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d" include("display/printing_juno.jl")
        init_defaults(default_file)
    end

    const pythonpath = Sys.iswindows() ? joinpath(depsdir, "python2", "python") : joinpath(dirname(@__DIR__), "deps", "python2", "bin", "python")
    const cif2cellpath = Sys.iswindows() ? joinpath(depsdir, "python2", "Scripts", "cif2cell") : joinpath(dirname(@__DIR__), "deps", "python2", "bin", "cif2cell")

end
