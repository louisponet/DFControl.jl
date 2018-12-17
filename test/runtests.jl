using DFControl, Test, Suppressor

testdir = joinpath(dirname(dirname(pathof(DFControl))), "test")
@time begin
@time @testset "constants"         begin @suppress include("constant_tests.jl")    end
@time @testset "Setting defaults"  begin @suppress include("defaults_test.jl")     end
@time @testset "Job from CIF file" begin @suppress include("jobfromcif_test.jl")   end
@time @testset "Job control tests" begin @suppress include("job_control_tests.jl") end
@time @testset "Remove defaults"   begin @suppress include("rmdefaults_test.jl")   end
# @testset "Plotting tests" begin include("plotting_tests.jl") end
end
