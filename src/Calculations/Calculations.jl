module Calculations
# This module handles all interactions with calculations
using Parameters, StructTypes, LinearAlgebra, JSON3 
using ..DFControl: config_path, DEPS_DIR, Band, Point3, Vec3, Mat3, Mat4, SVector, SArray
using ..Utils
using ..Structures
import ..Database: storage_directory, exists, load, save, Storable, verify

# include("modules.jl")
include("execs.jl")
include("calculation.jl")
include("qe.jl")
include("elk.jl")
include("wannier.jl")
include("abinit.jl")
include("documentation.jl")
include("julia.jl")

export Exec, Calculation, InputData
export Wannier90, QE, Abinit, Elk, Julia
end
