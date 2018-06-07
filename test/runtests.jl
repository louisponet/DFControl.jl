using Base.Test

tic()
@testset "constants" begin include("constant_tests.jl") end
@testset "Setting defaults"  begin include("defaults_test.jl") end
@testset "Job from CIF file" begin include("jobfromcif_test.jl") end
@testset "Job control tests" begin include("job_control_tests.jl") end
@testset "Remove defaults"   begin include("rmdefaults_test.jl") end
@testset "Plotting tests" begin include("plotting_tests.jl") end
toc()
