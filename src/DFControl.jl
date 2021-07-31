module DFControl
const CONFIG_DIR = occursin("cache", first(Base.DEPOT_PATH)) ?
                   abspath(Base.DEPOT_PATH[2], "config", "DFControl") :
                   abspath(Base.DEPOT_PATH[1], "config", "DFControl")

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

using Crayons
using REPL: REPL
using REPL.TerminalMenus
using CodeTracking
using Pkg
using LoggingExtras
using StructTypes

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
export list_pseudosets

const dfprintln = println
const dfprint = print
include("display/overloads.jl")

using DaemonMode, Pkg, LoggingExtras, Distributed
include("microservice/Mapper.jl")
include("microservice/Service.jl")
include("microservice/Resource.jl")
include("daemon.jl")

using Requires

function __init__()
    @require Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d" include("display/printing_juno.jl")
    @require Glimpse = "f6e19d58-12a4-5927-8606-ac30a9ce9b69" include("display/glimpse.jl")
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plotting.jl")
    merge!(Unitful.basefactors, localunits)
    Unitful.register(@__MODULE__)
    init_job_registry()
    # if !haskey(ENV, "IS_DAEMON")
    #     init_daemon()
    # else
    #     global_logger(DFControl.daemon_logger())
    # end
    return 
end

const pythonpath = Sys.iswindows() ? joinpath(depsdir, "python2", "python") :
                   joinpath(dirname(@__DIR__), "deps", "python2", "bin", "python")
const cif2cellpath = Sys.iswindows() ? joinpath(depsdir, "python2", "Scripts", "cif2cell") :
                     joinpath(dirname(@__DIR__), "deps", "python2", "bin", "cif2cell")

end
