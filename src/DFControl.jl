module DFControl
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
    using FileIO
    import FileIO: save

    using spglib_jll
    const SPGLIB = spglib_jll.libsymspg

	@reexport using Unitful
	import Unitful: Length, @unit, FreeUnits, unit, 𝐋, FreeUnits
	include("units.jl")

	using NearestNeighbors
    using Crayons


    abstract type Package end
    struct Wannier90 <: Package end
    struct QE        <: Package end
    struct Abinit    <: Package end
    struct Elk       <: Package end
    export Wannier90, QE, Abinit, Elk

    const depsdir = joinpath(dirname(@__DIR__), "deps")

    include("types.jl")
    include("atom.jl")
    include("structure.jl")

    include("input.jl")
    include("utils.jl")
    export yesterday, lastweek, lastmonth

    include("job.jl")

    export DFJob, Exec, DFInput, InputData
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
    export @add_default
    export setdefault_input
    export removedefault_input
    export getdefault_pseudo
    export getdefault_server
    export getdefault_jobheader
    export getdefault_pseudodir
    export getdefault_pseudodirs
    export removedefault

    const dfprintln = println
    const dfprint = print
    include("display/overloads.jl")
    using Requires
    function __init__()
        @require Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d" include("display/printing_juno.jl")
        init_defaults(default_file)
		merge!(Unitful.basefactors, localunits)
		Unitful.register(@__MODULE__)
    end

    const pythonpath = Sys.iswindows() ? joinpath(depsdir, "python2", "python") : joinpath(dirname(@__DIR__), "deps", "python2", "bin", "python")
    const cif2cellpath = Sys.iswindows() ? joinpath(depsdir, "python2", "Scripts", "cif2cell") : joinpath(dirname(@__DIR__), "deps", "python2", "bin", "cif2cell")
#     if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
#       include(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
#     else
#       error("DFControl not properly installed. Please run Pkg.build(\"DFControl\")")
# end
    # include("precompile.jl")
    # _precompile_()
end
