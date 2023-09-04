module Calculations
# This module handles all interactions with calculations
using Parameters, StructTypes, LinearAlgebra, JSON3
using RemoteHPC: Exec, exec
import RemoteHPC
using ..DFControl: DEPS_DIR, Band, Point3, Vec3, Mat3, Mat4, SVector, SArray
using ..Utils
using ..Structures

include("calculation.jl")
include("qe.jl")
include("elk.jl")
include("wannier.jl")
include("abinit.jl")
include("documentation.jl")
include("julia.jl")

export Calculation, InputData
export Wannier90, AbstractQE, QE, QE7_2, Abinit, Elk, Julia, set_flags!
end
