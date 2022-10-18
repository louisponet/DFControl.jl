module DFControl
const DFC = DFControl
export DFC

const CONFIG_DIR = occursin("cache", first(Base.DEPOT_PATH)) ?
                   abspath(Base.DEPOT_PATH[2], "config", "DFControl") :
                   abspath(Base.DEPOT_PATH[1], "config", "DFControl")
const DEPS_DIR = joinpath(dirname(@__DIR__), "deps")

const PYTHONPATH = Sys.iswindows() ? joinpath(DEPS_DIR, "python2", "python") :
                   joinpath(dirname(@__DIR__), "deps", "python2", "bin", "python")

const CIF2CELLPATH = Sys.iswindows() ?
                     joinpath(DEPS_DIR, "python2", "Scripts", "cif2cell") :
                     joinpath(dirname(@__DIR__), "deps", "python2", "bin", "cif2cell")

config_path(path...) = joinpath(CONFIG_DIR, gethostname(), path...)

using LinearAlgebra
using Reexport
@reexport using StaticArrays

using Parameters, StructTypes, Dates
include("types.jl")
export Point, Point3, Vec3, Mat3, Mat4

include("utils.jl")

include("Database/Database.jl")
include("Servers/Servers.jl")
@reexport using .Database

include("Structures/Structures.jl")
include("Calculations/Calculations.jl")
include("Jobs/Jobs.jl")
include("FileIO/FileIO.jl")

include("Service/Service.jl")
include("Resource/Resource.jl")
include("Client/Client.jl")

include("Display/Display.jl")

@reexport using .Client

using SnoopPrecompile

@precompile_all_calls begin
    if !exists(Server(name=gethostname()))
        Servers.configure_local(interactive=false)
    end
    s = local_server()
    alive = isalive(s)
    if !exists(Environment("default"))
        save(Environment(name="default"))
    end
end


end
