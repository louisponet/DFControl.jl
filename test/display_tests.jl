using DFControl

testjobpath = joinpath(testdir, "testassets", "test_job")

@testset "Basic printing" begin
    job = load(Job(testjobpath))
    show(job.calculations[1])
    show(job.calculations[1].data[1])
    show(DFC.Calculations.qe_execflag(:nosym))
    show(outputdata(job)["scf"][:bands][1])
    show(job.structure)
    show(job.structure.atoms[1])
    show(job.structure.atoms[1].element)
    show(Projection("d"))
    show(DFC.Calculations._ELK_CONTROLBLOCKS()[1])
    show(DFC.Calculations._ELK_CONTROLBLOCKS()[1].flags[1])
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
