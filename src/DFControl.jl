module DFControl
const DFC = DFControl
export DFC

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
include("utils.jl")
using .Utils

include("calculation.jl")
include("job.jl")
include("documentation.jl")

include("constants.jl")


include("defaults.jl")

const dfprintln = println
const dfprint = print
include("display/overloads.jl")

using Pkg, LoggingExtras, Distributed
include("Service/Service.jl")
include("Resource/Resource.jl")
include("Client/Client.jl")
export Client

using Requires

function __init__()
    @require Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d" include("display/printing_juno.jl")
    @require Glimpse = "f6e19d58-12a4-5927-8606-ac30a9ce9b69" include("display/glimpse.jl")
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plotting.jl")
    merge!(Unitful.basefactors, localunits)
    Unitful.register(@__MODULE__)
    Service.init_job_registry()
    Client.maybe_create_localhost()
    # Client.maybe_start_server("localhost")
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

# if ccall(:jl_generating_output, Cint, ()) == 1
#     DFJob(joinpath(dirname(dirname(pathof(DFControl))),"test/testassets/reference_job"))
# end
end
