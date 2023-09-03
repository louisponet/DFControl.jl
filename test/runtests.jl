using DFControl, Test
using RemoteHPC: configure_local, julia_main

# tconfdir = "/tmp/remotehpc"
tconfdir = tempname()
if ispath(tconfdir)
    rm(tconfdir; recursive = true)
end
import RemoteHPC: config_path
config_path(p...) = joinpath(tconfdir, p...)

paths = ["jobs",
         "logs/jobs",
         "logs/runtimes",
         "storage/servers",
         "storage/execs",
         "storage/environments"]
         
for p in paths
    mkpath(config_path(p))
end

redirect_stdin(devnull) do
    redirect_stderr(devnull) do
        redirect_stdout(devnull) do
            configure_local(; interactive = false)
            return t = @async julia_main()
        end
    end
end

while !isalive(local_server())
    sleep(0.1)
end

using UUIDs
testdir = @__DIR__
test_server = local_server()
@time begin
    # @testset "constants" begin
    #     include("constant_tests.jl")
    # end
    # @testset "documenation" begin
    #     include("documentation_tests.jl")
    # end
    @testset "Setting defaults" begin
        include("defaults_tests.jl")
    end
    # @testset "Job from CIF file" begin
    #     include("jobfromcif_tests.jl")
    # end
    @testset "Job from CIF file 7.2" begin
        include("jobfromcif_tests_7.2.jl")
    end
    @testset "Display tests" begin
        include("display_tests.jl")
    end
    @testset "Job control tests" begin
        include("job_control_tests.jl")
    end
    @testset "Remove defaults" begin
        include("rmdefaults_tests.jl")
    end
    include("cleanup.jl")
end
