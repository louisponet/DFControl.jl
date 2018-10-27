using DFControl, Test

testjobpath = joinpath(testdir, "testassets", "test_job")
job = DFJob(testjobpath);

#output stuff
out = outputdata(job;print=false,onlynew=false);
@test haskey(out, "nscf")
@test haskey(out["nscf"], :fermi)

nscf = DFControl.input(job, "nscf")
nscf2 = DFInput(nscf, "nscf2", data=[:testdata => (:testoption, "test"), :k_points => (:blabla, [1,1,1,1,1,1])])

@test data(nscf2, :testdata).option == :testoption
@test data(nscf2, :testdata).data   == "test"
@test data(nscf2, :testdata).name   == :testdata
@test data(nscf2, :k_points).option == :blabla
@test data(nscf2, :k_points).data   == [1,1,1,1,1,1]
@test job["scf"][:k_points].option == :automatic
@test_throws ErrorException job["bladbsflk"]
@test_throws ErrorException job["nscf"][:bladkfj]


setkpoints!(nscf2, (3,3,3), print=false)
@test data(nscf2, :k_points).data  == kgrid(3, 3, 3, :nscf)


setkpoints!(nscf2, [(3.,3.,3.,1.), (3.,3.,3.,1.)], print=false)
@test data(nscf2, :k_points).option  == :crystal_b
@test data(nscf2, :k_points).data  == [(3.,3.,3.,1.), (3.,3.,3.,1.)]

setkpoints!(nscf2, (3,3,3,0,0,1), print=false)
@test data(nscf2, :k_points).option  == :automatic
@test data(nscf2, :k_points).data  == [3,3,3,0,0,1]

fermi = read_qe_output(outpath(job, "nscf"))[:fermi]
@test fermi == read_fermi_from_qe_output(joinpath(job.local_dir, "nscf.out"))

addcalc!(job, [(0.5,0.0,0.5,10.0),(0.0,0.0,0.0,10.0),(0.5,0.5,0.5,1.0)], name="bands2")
@test flag(job, "bands2", :calculation) == "'bands'"
@test data(job, "bands2", :k_points).data == [(0.5,0.0,0.5,10.0),(0.0,0.0,0.0,10.0),(0.5,0.5,0.5,1.0)]
addcalc!(job, (10,10,10), name="nscf2")
@test flag(job, "nscf2", :calculation) == "'nscf'"
addcalc!(job, (5,5,5,1,1,1), name="1scf2")
@test flag(job, "1scf2", :calculation) == "'scf'"
@test job["nscf2"][:calculation] == flag(job, "nscf2", :calculation)


wanflags = [:write_hr => true, :wannier_plot => true]

addwancalc!(job, "nscf",:Pt => [:s, :p, :d], Emin=fermi-7.0, Epad=5.0, wanflags=wanflags, print=false)
@test flag(job, "wan", :write_hr) == flag(job, "wan", :wannier_plot) == true


job.inputs = job.inputs[1:3]

setflags!(job, :nspin => 2, print=false)
@test flag(job, "nscf", :nspin) == 2
job["nscf"][:nspin] = 3
@test job["nscf"][:nspin] == 3

job[:nspin] = 10
@test job["nscf"][:nspin] == 10

job["nscf"][:Hubbard_U] = [1.0, 2.0]
@test job["nscf"][:Hubbard_U] == [1.0, 2.0]


addwancalc!(job, nscf,:Pt => [:s, :p, :d], Emin=fermi-7.0, Epad=5.0, wanflags=wanflags, print=false)

setflow!(job, ""=>false)
@test job.inputs[1].run == false
setflow!(job, "nscf" => true, "bands" => true)
@test job.inputs[3].run
setflow!(job, "pw2wan" => true)
@test job.inputs[end].run

@test inputs(job,["nscf"]) == inputs(job, "nscf")

save(job)

job2 = DFJob(local_dir)

begin
    for (calc, calc2) in zip(job.inputs, job2.inputs)
        @test calc.flags == calc2.flags
        for (b1, b2) in zip(calc.data, calc2.data)
            for n in fieldnames(typeof(b1))
                @test getfield(b1, n) == getfield(b2, n)
            end
        end
    end
end

job3 = DFJob(job2, :lspinorb => true)
@test all(atoms(job3).==atoms(job2))
@test flag(job3, :lspinorb)
rmflags!(job3, :lspinorb, print=false)

begin
    for (calc, calc2) in zip(job.inputs, job3.inputs)
        @test calc.flags == calc2.flags
        for (b1, b2) in zip(calc.data, calc2.data)
            for n in fieldnames(typeof(b1))
                @test getfield(b1, n) == getfield(b2,n)
            end
        end
    end
end

testorbs = [:s, :p]
setprojections!(job, :Pt => testorbs)
@test convert.(Symbol, [p.orb for p in projections(job, :Pt)]) == testorbs
setwanenergies!(job, fermi-7.0, read_qe_bands_file(outpath(nscf)), Epad=3.0, print=false)

@test flag(job, :dis_froz_max) == 14.285699999999999
@test flag(job, :dis_win_max) == 14.285699999999999 + 3.0

setexecflags!(job, "pw.x", :nk => 230)
@test execs(job, "nscf")[2].flags[:nk] == 230
rmexecflags!(job, "pw.x", :nk)
@test !haskey(execs(job, "nscf")[2].flags,:nk)

setexecdir!(job, "pw.x", joinpath(homedir(), "bin"))
@test execs(job, "nscf")[2].dir == joinpath(homedir(), "bin")

setname!(job, "nscf", "test")
@test inpath(job, "test") == joinpath(job.local_dir, "test.in")
setname!(job, "test", "nscf")

setserverdir!(job, "localhost")
@test job.server_dir == "localhost"

setheaderword!(job, "defpart", "frontend", print=false)
@test any(occursin.("frontend",job.header))

undo!(job)
@test any(occursin.("defpart", job.header))


setdataoption!(job, "nscf",:k_points, :blabla, print=false)
@test data(job, "nscf", :k_points).option == :blabla
setdataoption!(job, :k_points, :test, print=false)
@test data(job, "nscf", :k_points).option == :test

job = undo(job)
@test data(job, "nscf", :k_points).option == :blabla

report = progressreport(job; onlynew=false, print=false)
@test report[:fermi] == 17.4572
@test length(report[:accuracy]) == 9

rm.(inpath.(job.inputs))
rm(joinpath(job.local_dir, "job.tt"))
