using DFControl, Test
using UUIDs
testdir = @__DIR__
@time begin
    if exists(Server(name="localhost"))
        test_server = Server("localhost")
        if Servers.isalive(test_server)
            try
                Servers.kill(test_server)
            catch
                nothing
            end
        end
    else
        test_server = Server(name="localhost", domain="localhost", julia_exec = joinpath(Sys.BINDIR, "julia"), uuid=string(uuid4()))
        save(test_server) 
    end
        # Servers.save(test_server)
    # DFControl.Servers.start(test_server)
    # if Servers.isalive(test_server)
    #     try
    #         Servers.kill_server(test_server)
    #     catch
    #         nothing
    #     end
    # end
    @async DFControl.Resource.run()
    while !Servers.isalive(Server("localhost"))
        sleep(0.1)
    end
    testserver = Server("localhost")
    @time @testset "constants" begin
        include("constant_tests.jl")
    end
    @time @testset "documenation" begin
        include("documentation_tests.jl")
    end
    @time @testset "Setting defaults" begin
        include("defaults_tests.jl")
    end
    @time @testset "Database" begin
        include("database_tests.jl")
    end
    @time @testset "Job from CIF file" begin
        include("jobfromcif_tests.jl")
    end
    @time @testset "Display tests" begin
        include("display_tests.jl")
    end
    @time @testset "Job control tests" begin
        include("job_control_tests.jl")
    end
    @time @testset "Remove defaults" begin
        include("rmdefaults_tests.jl")
    end
    include("cleanup.jl")
end
