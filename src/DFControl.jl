module DFControl
const DFC = DFControl
export DFC

const DEPS_DIR = joinpath(dirname(@__DIR__), "deps")

const PYTHONPATH = Sys.iswindows() ? joinpath(DEPS_DIR, "python2", "python") :
                   joinpath(dirname(@__DIR__), "deps", "python2", "bin", "python")

const CIF2CELLPATH = Sys.iswindows() ?
                     joinpath(DEPS_DIR, "python2", "Scripts", "cif2cell") :
                     joinpath(dirname(@__DIR__), "deps", "python2", "bin", "cif2cell")

using LinearAlgebra
using Reexport
@reexport using StaticArrays

using Parameters, StructTypes, Dates
include("types.jl")
export Point, Point3, Vec3, Mat3, Mat4

include("utils.jl")

include("Structures/Structures.jl")
include("Calculations/Calculations.jl")
include("Jobs/Jobs.jl")
include("FileIO/FileIO.jl")
include("Client/Client.jl")
include("Display/Display.jl")

@reexport using .Client

end
