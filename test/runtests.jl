using DFControl, Test, Suppressor

testdir = @__DIR__ 
@time begin
@time @testset "constants"         begin include("constant_tests.jl")      end
@time @testset "documenation"      begin @suppress include("documentation_tests.jl") end
@time @testset "Setting defaults"  begin @suppress include("defaults_tests.jl")      end
@time @testset "Job from CIF file" begin @suppress include("jobfromcif_tests.jl")    end
@time @testset "Job control tests" begin @suppress include("job_control_tests.jl")   end
@time @testset "Remove defaults"   begin @suppress include("rmdefaults_tests.jl")    end
end
