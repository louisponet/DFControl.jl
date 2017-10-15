using DFControl
using Base.Test




tic()
@testset "File processing tests" begin include("file_processing_tests.jl") end
@testset "Job control tests" begin include("job_control_tests.jl") end
toc()