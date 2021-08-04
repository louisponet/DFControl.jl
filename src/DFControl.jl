module DFControl
const DFC = DFControl
export DFC

const CONFIG_DIR = occursin("cache", first(Base.DEPOT_PATH)) ?
                   abspath(Base.DEPOT_PATH[2], "config", "DFControl") :
                   abspath(Base.DEPOT_PATH[1], "config", "DFControl")
const DEPS_DIR = joinpath(dirname(@__DIR__), "deps")

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

const Point{N,F} = SVector{N,F}
const Point3{F}  = SVector{3,F}
const Vec3{F}    = SVector{3,F}
const Mat3{F}    = SMatrix{3,3,F,9}
const Mat4{F}    = SMatrix{4,4,F,16}
const Vec{N,T}   = SVector{N,T}

Point{N,T}(x::T) where {N,T} = Point{N,T}(x, x, x)

function Base.convert(::Type{Point3{T}}, x::Vector{T}) where {T<:AbstractFloat}
    return Point3{T}(x[1], x[2], x[3])
end

export Point, Point3, Vec3, Mat3, Mat4


using Parameters
using JLD2

using Crayons
using REPL: REPL
using REPL.TerminalMenus
using CodeTracking
using Pkg
using LoggingExtras
using StructTypes


include("types.jl")
include("utils.jl")
using .Utils

include("Structures/Structures.jl")
include("Calculations/Calculations.jl")
include("Jobs/Jobs.jl")
include("FileIO/FileIO.jl")


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

const pythonpath = Sys.iswindows() ? joinpath(DEPS_DIR, "python2", "python") :
                   joinpath(dirname(@__DIR__), "deps", "python2", "bin", "python")
const cif2cellpath = Sys.iswindows() ? joinpath(DEPS_DIR, "python2", "Scripts", "cif2cell") :
                     joinpath(dirname(@__DIR__), "deps", "python2", "bin", "cif2cell")

# if ccall(:jl_generating_output, Cint, ()) == 1
#     Job(joinpath(dirname(dirname(pathof(DFControl))),"test/testassets/reference_job"))
# end
end
