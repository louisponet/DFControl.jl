using DFControl
using Base.Test




tic()
@testset "Type tests" begin
  t_p = Point3(0.2)
  @test t_p + Point3(1.) == Point3(1.2)
  @test t_p - Point3(1.) == Point3(-0.8)
  @test t_p * 2 == 2 * t_p
  @test t_p/t_p == Point3(1.2)/Point3(1.2)
  @test zero(Point3) == Point3(0.)
  @test norm(zero(Point3)) == 0.0
  @test convert(Point3,Point3(1.2)) == Point3(1.2)
end



@testset "File processing tests" begin include("file_processing_tests.jl") end
@testset "Job control tests" begin include("job_control_tests.jl") end
@testset "Plotting tests" begin include("plotting_tests.jl") end
toc()