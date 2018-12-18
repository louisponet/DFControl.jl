using DFControl, Test

testjobpath = joinpath(testdir, "testassets", "test_job")
job = DFJob(testjobpath);

#output stuff
out = outputdata(job;print=false,onlynew=false);
@test haskey(out, "nscf")
@test haskey(out["nscf"], :fermi)

@test isconverged(job["scf"])
nscf = DFControl.input(job, "nscf")
show(nscf)
nscf2 = DFInput(nscf, "nscf2", data=[:testdata => (:testoption, "test"), :k_points => (:blabla, [1,1,1,1,1,1])])

@test data(nscf2, :testdata).option == :testoption
@test data(nscf2, :testdata).data   == "test"
@test data(nscf2, :testdata).name   == :testdata
@test data(nscf2, :k_points).option == :blabla
@test data(nscf2, :k_points).data   == [1,1,1,1,1,1]
@test data(job["scf"], :k_points).option == :automatic
@test_throws ErrorException job["bladbsflk"]
@test_throws ErrorException nscf2[:bladkfj]


setkpoints!(nscf2, (3,3,3), print=false)
@test data(nscf2, :k_points).data  == DFControl.kgrid(3, 3, 3, :nscf)


setkpoints!(nscf2, [(3.,3.,3.,1.), (3.,3.,3.,1.)], print=false)
@test data(nscf2, :k_points).option  == :crystal_b
@test data(nscf2, :k_points).data  == [(3.,3.,3.,1.), (3.,3.,3.,1.)]

setkpoints!(nscf2, (3,3,3,0,0,1), print=false)
@test data(nscf2, :k_points).option  == :automatic
@test data(nscf2, :k_points).data  == [3,3,3,0,0,1]

fermi = outputdata(job, "nscf")[:fermi]
@test fermi == 17.4572

push!(job, gencalc_bands(job["bands"], [(0.5,0.0,0.5,10.0),(0.0,0.0,0.0,10.0),(0.5,0.5,0.5,1.0)], name="bands2"))
@test job["bands2"][:calculation] == "bands"
@test data(job, "bands2", :k_points).data == [(0.5,0.0,0.5,10.0),(0.0,0.0,0.0,10.0),(0.5,0.5,0.5,1.0)]
push!(job, gencalc_nscf(job["nscf"], (10,10,10), name="nscf2"))
@test job["nscf2"][:calculation] == "nscf"
push!(job, gencalc_scf(job["scf"], (5,5,5,1,1,1), name="1scf2"))
@test job["1scf2"][:calculation] == "scf"


wanflags = [:write_hr => true, :wannier_plot => true]

setprojections!(job, :Pt => [:s, :p, :d])
wancalc = gencalc_wan(structure(job), searchinput(job, "nscf"), fermi-7.0; Epad=5.0, wanflags=wanflags)
push!.((job,), wancalc)
@test job["wan"][:write_hr] == job["wan"][:wannier_plot] == true

wanout = outputdata(job, "wan")
@test length(wanout[:final_state]) == 20
@test length(wanout[:wannierise]) == 203
job.inputs = job.inputs[1:4]

setflags!(job, :nspin => 2, print=false)
@test job["nscf"][:nspin] == 2
job["nscf"][:nspin] = 3
@test job["nscf"][:nspin] == 3

job[:nspin] = 10
@test job["nscf"][:nspin] == 10

job["nscf"][:Hubbard_U] = [1.0]
@test job["nscf"][:Hubbard_U] == [1.0]

job["nscf"][:starting_magnetization] = [[1.0, 2.0, 1.0]]
@test job["nscf"][:starting_magnetization] == [[1.0, 2.0, 1.0]]

job["nscf"][:Hubbard_J] = [0 1 2]
@test job["nscf"][:Hubbard_J] == [0 1 2]
push!.((job,), gencalc_wan(structure(job), nscf, fermi-7.0, Epad=5.0, wanflags=wanflags))

#TODO: add test with the new wancalc from projwfc

setflow!(job, ""=>false)
@test job.inputs[1].run == false
setflow!(job, "nscf" => true, "bands" => true)
@test job.inputs[3].run


save(job)
job2 = DFJob(local_dir)


begin
    for (calc, calc2) in zip(job.inputs, job2.inputs)

        for (f, v) in calc.flags
            if f in (:Hubbard_J, :pseudo_dir)
                continue
            else
                @test v == calc2.flags[f]
            end
        end
        for (b1, b2) in zip(calc.data, calc2.data)
            for n in fieldnames(typeof(b1))
                @test getfield(b1, n) == getfield(b2, n)
            end
        end
    end
end

job3 = DFJob(job2, :lspinorb => true)
@test all(atoms(job3).==atoms(job2))
@test length(job3[:lspinorb]) == length(searchinputs(job3, QE))
rmflags!(job3, :lspinorb, print=false)

begin
    for (calc, calc2) in zip(job.inputs, job3.inputs)
        for (f, v) in calc.flags
            if f in (:Hubbard_J, :pseudo_dir)
                continue
            else
                @test v == calc2.flags[f]
            end
        end
        for (b1, b2) in zip(calc.data, calc2.data)
            for n in fieldnames(typeof(b1))
                @test getfield(b1, n) == getfield(b2,n)
            end
        end
    end
end

setcutoffs!(job)
@test job["scf"][:ecutwfc] == 32.0

setpseudos!(job, "pseudos", :Pt => "Pt.UPF")
@test job.structure.atoms[1].pseudo == "Pt.UPF"
@test job["nscf"][:pseudo_dir] == "pseudos"
testorbs = [:s, :p]
setprojections!(job, :Pt => testorbs)
@test convert.(Symbol, [p.orb for p in projections(job, :Pt)]) == testorbs
setwanenergies!(job, nscf, fermi-7.0, Epad=3.0)

@test job["wanup"][:dis_froz_max] == 13.2921
@test job["wanup"][:dis_win_max] == 16.292099999999998

setexecflags!(job, "pw.x", :nk => 230)
@test DFControl.getfirst(x->x.symbol==:nk, execs(job, "nscf")[2].flags).value == 230
rmexecflags!(job, "pw.x", :nk)
@test isempty(filter(x->x.symbol == :nk, execs(job, "nscf")[2].flags))

setexecdir!(job, "pw.x", joinpath(homedir(), "bin"))
@test execs(job, "nscf")[2].dir == joinpath(homedir(), "bin")

setname!(job, "nscf", "test")
@test DFControl.inpath(job, "test") == joinpath(job.local_dir, "test.in")
setname!(job, "test", "nscf")

setserverdir!(job, "localhost")
@test job.server_dir == "localhost"

setheaderword!(job, "defpart", "frontend", print=false)
@test any(occursin.("frontend",job.header))

setdataoption!(job, "nscf",:k_points, :blabla, print=false)
@test data(job, "nscf", :k_points).option == :blabla
setdataoption!(job, :k_points, :test, print=false)
@test data(job, "nscf", :k_points).option == :test


report = progressreport(job; onlynew=false, print=false)
@test report[:fermi] == 17.4572
@test length(report[:accuracy]) == 9
newatompos = outputdata(job, "vc_relax", onlynew=false)[:final_structure]
job.structure = newatompos


rm.(DFControl.inpath.(job.inputs))

rm(joinpath(splitdir(DFControl.inpath(job.inputs[1]))[1], "pw2wan_wanup.in"))
rm(joinpath(splitdir(DFControl.inpath(job.inputs[1]))[1], "pw2wan_wandn.in"))
rm(joinpath(job.local_dir, "job.tt"))
