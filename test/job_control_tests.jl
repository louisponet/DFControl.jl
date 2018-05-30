using DFControl, Base.Test

import DFControl:  data
# testjobpath = joinpath(Pkg.dir("DFControl"), "test/testassets/test_job/")
testjobpath = joinpath(@__DIR__, "testassets/test_job/")
job = DFJob(testjobpath);

nscf = DFControl.input(job, "nscf")
nscf2 = DFInput(nscf, "nscf2.in", data=[:testdata => (:testoption, "test"), :k_points => (:blabla, [1,1,1,1,1,1])])

@test data(job, "nscf2", :testdata) == nothing
@test data(nscf2, :testdata).option == :testoption
@test data(nscf2, :testdata).data   == "test"
@test data(nscf2, :testdata).name   == :testdata
@test data(nscf2, :k_points).option == :blabla
@test data(nscf2, :k_points).data   == [1,1,1,1,1,1]

setkpoints!(nscf2, (3,3,3), print=false)
@test data(nscf2, :k_points).data  == kgrid(3, 3, 3, :nscf)


setkpoints!(nscf2, [(3.,3.,3.,1.), (3.,3.,3.,1.)], print=false)
@test data(nscf2, :k_points).option  == :crystal_b
@test data(nscf2, :k_points).data  == [(3.,3.,3.,1.), (3.,3.,3.,1.)]

setkpoints!(nscf2, (3,3,3,0,0,1), print=false)
@test data(nscf2, :k_points).option  == :automatic
@test data(nscf2, :k_points).data  == [3,3,3,0,0,1]

fermi = read_qe_output(outpath(job, "nscf"))[:fermi]
@test fermi == read_fermi_from_qe_output(job.local_dir * "nscf.out")

addbandscalc!(job, [[0.5,0.0,0.5,10.0],[0.0,0.0,0.0,10.0],[0.5,0.5,0.5,1.0]],filename="bands2")
@test flag(job, "bands2", :calculation) == "'bands'"
@test data(job, "bands2", :k_points).data == [[0.5,0.0,0.5,10.0],[0.0,0.0,0.0,10.0],[0.5,0.5,0.5,1.0]]

wanflags = [:write_hr => true, :wannier_plot => true]

addwancalc!(job, nscf,[:Pt => [:s, :p, :d]], Emin=fermi-7.0, Epad=5.0, wanflags=wanflags, print=false)

@test flag(job, "wan.win", :write_hr) == flag(job, "wan.win", :wannier_plot) == true

job.inputs = job.inputs[1:3]

setflags!(job, :nspin => 2, print=false)
@test flag(job, "nscf", :nspin) == 2

addwancalc!(job, nscf,[:Pt => [:s, :p, :d]], Emin=fermi-7.0, Epad=5.0, wanflags=wanflags, print=false)


setflow!(job, [false for i=1:length(job.inputs)])
@test job.inputs[1].run == false
setflow!(job, "nscf" => true, "bands" => true)
@test job.inputs[3].run
setflow!(job, "pw2wan" => true)
@test job.inputs[end].run

print_info(job)
@test inputs(job,["nscf"]) == inputs(job, "nscf")

save(job)

job2 = DFJob(local_dir)

begin
    for (calc, calc2) in zip(job.inputs, job2.inputs)
        @test calc.flags == calc2.flags
        for (b1, b2) in zip(calc.data, calc2.data)
            for name in fieldnames(b1)
                @test getfield(b1, name) == getfield(b2,name)
            end
        end
    end
end

job3 = DFJob(job2, :lspinorb => true)
@test flag(job3, :lspinorb)
rmflags!(job3, :lspinorb, print=false)

begin
    for (calc, calc2) in zip(job.inputs, job3.inputs)
        @test calc.flags == calc2.flags
        for (b1, b2) in zip(calc.data, calc2.data)
            for name in fieldnames(b1)
                @test getfield(b1, name) == getfield(b2,name)
            end
        end
    end
end

testorbs = [:s, :p]
setprojections!(job, :Pt => testorbs)
@test convert.(Symbol, [p.orb for p in projections(job, :Pt)]) == testorbs
setwanenergies!(job, fermi-7.0, read_qe_bands_file(outpath(job, nscf)), Epad=3.0, print=false)

@test flag(job, :dis_froz_max) == 14.2907
@test flag(job, :dis_win_max) == 14.2907 + 3.0

setexecflags!(job, "pw.x", :nk => 230)
@test execs(job, "nscf")[2].flags[:nk] == 230

setexecdir!(job, "pw.x", "~/bin/")
@test execs(job, "nscf")[2].dir == "~/bin/"


rm.(path.(job, job.inputs))
rm(job.local_dir * "job.tt")
