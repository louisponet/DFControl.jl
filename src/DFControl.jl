module DFControl
const DFC = DFControl
export DFC

const CONFIG_DIR = occursin("cache", first(Base.DEPOT_PATH)) ?
                   abspath(Base.DEPOT_PATH[2], "config", "DFControl") :
                   abspath(Base.DEPOT_PATH[1], "config", "DFControl")
const DEPS_DIR = joinpath(dirname(@__DIR__), "deps")

const PYTHONPATH = Sys.iswindows() ? joinpath(DEPS_DIR, "python2", "python") :
                   joinpath(dirname(@__DIR__), "deps", "python2", "bin", "python")
                   
const CIF2CELLPATH = Sys.iswindows() ? joinpath(DEPS_DIR, "python2", "Scripts", "cif2cell") :
                     joinpath(dirname(@__DIR__), "deps", "python2", "bin", "cif2cell")

config_path(path...) = joinpath(CONFIG_DIR, path...)

using LinearAlgebra
using Reexport
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

using Parameters, StructTypes, Dates
include("types.jl")
include("utils.jl")
include("Structures/Structures.jl")
include("Calculations/Calculations.jl")
include("Jobs/Jobs.jl")
include("FileIO/FileIO.jl")

include("Servers/Servers.jl")
include("Service/Service.jl")
include("Resource/Resource.jl")
include("Client/Client.jl")
include("Display/Display.jl")
export Client


function __init__()
    Jobs.init_job_registry()
    Servers.maybe_create_localhost()
    # Client.maybe_start_server("localhost")
    # if !haskey(ENV, "IS_DAEMON")
    #     init_daemon()
    # else
    #     global_logger(DFControl.daemon_logger())
    # end
    return 
end

# if ccall(:jl_generating_output, Cint, ()) == 1
#     Job(joinpath(dirname(dirname(pathof(DFControl))),"test/testassets/reference_job"))
# end
end
