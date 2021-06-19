module DFControl
    const CONFIG_DIR = occursin("cache", first(Base.DEPOT_PATH)) ?
                       abspath(Base.DEPOT_PATH[2], "config","DFControl") :
                       abspath(Base.DEPOT_PATH[1], "config","DFControl")
                       
    config_path(path...) = joinpath(CONFIG_DIR, path...)
    
    using LinearAlgebra
    using Statistics
    using Media
    using Dates
    using DelimitedFiles
    using Reexport
    import Base.Iterators.flatten

    using RecipesBase
    @reexport using StaticArrays
    using Parameters
    using JLD2
    
    using spglib_jll
    const SPGLIB = spglib_jll.libsymspg

	@reexport using Unitful
	import Unitful: Length, @unit, FreeUnits, unit, 𝐋, FreeUnits
	include("units.jl")

	using NearestNeighbors
    using Crayons
    import REPL
    using REPL.TerminalMenus
    
    abstract type Package end
    struct Wannier90 <: Package end
    struct QE        <: Package end
    struct Abinit    <: Package end
    struct Elk       <: Package end
    export Wannier90, QE, Abinit, Elk

    const depsdir = joinpath(dirname(@__DIR__), "deps")

    include("typedefs.jl")
    export Vec3, Point3
    include("types.jl")
    export DFJob, Exec, DFCalculation, InputData
    export Structure, Atom, Pseudo, DFTU, Orbital, Projection
    
    include("atom.jl")
    include("structure.jl")
    
    include("execs.jl")
    include("calculation.jl")
    include("utils.jl")
    export yesterday, lastweek, lastmonth

    include("job.jl")
    include("versioning.jl")
    include("registry.jl")

    include("server.jl")
    include("API.jl")

    include("constants.jl")

    include("fileio.jl")
    include("plotting.jl")

    include("defaults.jl")
    export setdefault_pseudodir
    export setdefault_server
    export configuredefault_pseudos
    export removedefault_pseudodir
    export removedefault_pseudos
    export setdefault_jobheader
    export getdefault_pseudo
    export getdefault_server
    export getdefault_jobheader
    export getdefault_pseudodir
    export getdefault_pseudodirs

    const dfprintln = println
    const dfprint = print
    include("display/overloads.jl")
    using Requires
    
    function __init__()
        @require Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d" include("display/printing_juno.jl")
        @require Glimpse = "f6e19d58-12a4-5927-8606-ac30a9ce9b69" include("display/glimpse.jl")
		merge!(Unitful.basefactors, localunits)
		Unitful.register(@__MODULE__)
		init_job_registry()
    end

    const pythonpath = Sys.iswindows() ? joinpath(depsdir, "python2", "python") : joinpath(dirname(@__DIR__), "deps", "python2", "bin", "python")
    const cif2cellpath = Sys.iswindows() ? joinpath(depsdir, "python2", "Scripts", "cif2cell") : joinpath(dirname(@__DIR__), "deps", "python2", "bin", "cif2cell")
    
end
