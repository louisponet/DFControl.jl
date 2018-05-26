using DFControl, Base.Test

# testjobpath = joinpath(Pkg.dir("DFControl"), "test/testassets/test_job/")
testjobpath = joinpath(@__DIR__, "testassets/test_job/")
job = DFJob(testjobpath);

nscf = DFControl.input(job, "nscf")
fermi = read_qe_output(outpath(job, "nscf"))[:fermi]
@test fermi == read_fermi_from_qe_output(job.local_dir * "nscf.out")

wanflags = [:write_hr => true, :wannier_plot => true]

addwancalc!(job, nscf,[:Pt => [:s, :p, :d]], Emin=fermi-7.0, Epad=5.0, wanflags=wanflags)

@test length(job.calculations) == 5
@test flag(job, "wan.win", :write_hr) == flag(job, "wan.win", :wannier_plot) == true

job.calculations = job.calculations[1:3]

setflags!(job, :nspin => 2, print=false)
@test flag(job, "nscf", :nspin) == 2

addwancalc!(job, nscf,[:Pt => [:s, :p, :d]], Emin=fermi-7.0, Epad=5.0, wanflags=wanflags)

@test length(job.calculations) == 7

setflow!(job, [false for i=1:length(job.calculations)])
@test job.calculations[1].run == false
setflow!(job, "nscf" => true, "bands" => true)
@test job.calculations[3].run
setflow!(job, "pw2wan" => true)
@test job.calculations[end].run

print_info(job)
@test inputs(job,["nscf"]) == inputs(job, "nscf")

save(job)

job2 = DFJob(local_dir)

@test length(job2.calculations) == length(job.calculations) == 7
begin
    for (calc, calc2) in zip(job.calculations, job2.calculations)
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
setwanenergies!(job, fermi-7.0, read_qe_bands_file(outpath(job, nscf)), Epad=3.0)

@test flag(job, :dis_froz_max) == 17.3737
@test flag(job, :dis_win_max) == 17.3737 + 3.0

rm.(path.(job, job.calculations))
rm(job.local_dir * "job.tt")
