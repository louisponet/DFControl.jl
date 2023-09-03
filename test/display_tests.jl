using DFControl

testjobpath = joinpath(testdir, "testassets", "test_job")

@testset "Basic printing" begin
    job = load(test_server, Job(testjobpath))
    display(job.calculations[1])
    display(job.calculations[1].data[1])
    display(outputdata(job)["scf"][:bands][1])
    display(job.structure)
    display(job.structure.atoms[1])
    display(job.structure.atoms[1].element)
    display(Projection("d"))
    display(DFC.Calculations._ELK_CONTROLBLOCKS()[1])
    display(DFC.Calculations._ELK_CONTROLBLOCKS()[1].flags[1])
end

# @testset "html printing" begin
    # job = Job(testjobpath)
    # display(stdout, MIME("text/html"), job.calculations[1])
# end

# using Plots

# @testset "Plots" begin
    # job = Job(testjobpath)
    # plot(job, -2, 2)
# end
