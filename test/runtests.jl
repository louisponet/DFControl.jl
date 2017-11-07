using DFControl
using Base.Test




tic()
@testset "Type tests" begin
  t_p = Point3D(0.2)
  @test t_p + 1 == Point3D(1.2)
  @test t_p - 1 == Point3D(-0.8)
  @test t_p * 2 == 2 * t_p
  @test t_p/t_p == Point3D(1.2)/Point3D(1.2)
  @test zero(Point3D) == Point3D(0.)
  @test norm(zero(Point3D)) == 0.0
  @test convert(Point3D,Point3D(1.2)) == Point3D(1.2)
  



@testset "File processing tests" begin include("file_processing_tests.jl") end
@testset "Job control tests" begin include("job_control_tests.jl") end
@testset "Plotting tests" begin include("plotting_tests.jl") end
toc()