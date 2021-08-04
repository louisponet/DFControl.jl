module Jobs
    # This module handles all the Job related functionality
    using Parameters, StructTypes
    using ..Calculations
    using ..Structures
    using ..Utils
    include("job.jl")
    include("versioning.jl")
    include("registry.jl")
    export Job
end
