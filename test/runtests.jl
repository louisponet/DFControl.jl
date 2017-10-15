using DFControl
using Base.Test




tic()
@testset "File processing tests" begin include("file_processing_tests.jl") end
toc()