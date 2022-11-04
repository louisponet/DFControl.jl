module Jobs
# This module handles all the Job related functionality
using Pkg
using Parameters, StructTypes, JLD2, Dates, JSON3
using ..Calculations
using ..Calculations: Calculation
using ..Structures
using ..Utils
import RemoteHPC
using RemoteHPC: Server, isalive

include("job.jl")
include("versioning.jl")
export Job, set_flow!
end
