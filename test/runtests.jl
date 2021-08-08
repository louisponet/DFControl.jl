using DFControl, Test, Suppressor

testdir = @__DIR__
@time begin
    test_server = Server("localhost_test", "ponet", "localhost", 8123, Servers.Bash, "", "julia", homedir() )
    if Servers.isalive(test_server)
        try
            Servers.kill_server(test_server)
        catch
            nothing
        end
    end
    Servers.save(test_server)
    DFControl.Servers.start(test_server)
    @time @testset "constants" begin
        @suppress include("constant_tests.jl")
    end
    @time @testset "documenation" begin
        @suppress include("documentation_tests.jl")
    end
    @time @testset "Setting defaults" begin
        include("defaults_tests.jl")
    end
    @time @testset "Job from CIF file" begin
        include("jobfromcif_tests.jl")
    end
    @time @testset "Job control tests" begin
        include("job_control_tests.jl")
    end
    @time @testset "Remove defaults" begin
        @suppress include("rmdefaults_tests.jl")
    end
    try
        Servers.kill_server(test_server)
    catch
        nothing
    end
    rm(DFControl.config_path("servers/localhost_test.jld2"))
end
