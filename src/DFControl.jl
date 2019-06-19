module DFControl
    using LinearAlgebra
    using Statistics
    using Media
    using Dates
    using DelimitedFiles
    import Base.Iterators.flatten

    using RecipesBase
    @reexport using StaticArrays
    using Parameters

	@reexport using Unitful
	Unitful.register(@__MODULE__);
	import Unitful: Length, @unit, FreeUnits, unit

	Base.eltype(::Type{Length{T}}) where T = T
	@unit Ang "Ang" Angstrom              1e-1u"nm"               false
    @unit e₀  "eₒ"  ElementaryCharge      1.602176620898e-19*u"C" false
    @unit kₑ  "kₑ"  CoulombForceConstant  1/(4π)u"ϵ0"             false
    @unit a₀  "a₀"  BohrRadius            1u"ħ^2/(1kₑ*me*e₀^2)"   false
    @unit Eₕ  "Eₕ"  HartreeEnergy         1u"me*e₀^4*kₑ^2/(1ħ^2)" true
    @unit Ry  "Ry"  RydbergEnergy         0.5Eₕ                   true


	@inline function StaticArrays._inv(::StaticArrays.Size{(3,3)}, A::SMatrix{3,3, LT}) where {LT<:Length}

	    @inbounds x0 = SVector{3}(A[1], A[2], A[3])
	    @inbounds x1 = SVector{3}(A[4], A[5], A[6])
	    @inbounds x2 = SVector{3}(A[7], A[8], A[9])

	    y0 = cross(x1,x2)
	    d  = StaticArrays.bilinear_vecdot(x0, y0)
	    x0 = x0 / d
	    y0 = y0 / d
	    y1 = cross(x2,x0)
	    y2 = cross(x0,x1)

	    @inbounds return SMatrix{3, 3}((y0[1], y1[1], y2[1], y0[2], y1[2], y2[2], y0[3], y1[3], y2[3]))
	end

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

    export DFJob, Exec, DFInput
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
