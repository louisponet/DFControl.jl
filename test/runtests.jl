using DFControl
using Base.Test
import DFControl: Exec, data, add!, execs



tic()
@testset "File processing tests" begin include("file_processing_tests.jl") end
@testset "Job control tests" begin include("job_control_tests.jl") end
#@testset "Plotting tests" begin include("plotting_tests.jl") end
toc()
