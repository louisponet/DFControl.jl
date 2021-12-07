module Jobs
# This module handles all the Job related functionality
using Parameters, StructTypes, CodeTracking, JLD2, Dates, JSON3
using ..Calculations
using ..Structures
using ..Utils
using ..DFControl

@enum JobState Pending Running Completed Cancelled 

include("environment.jl")
include("job.jl")
include("versioning.jl")
include("registry.jl")
export Job, Environment
end
