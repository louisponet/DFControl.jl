using DFControl, Test

testjobpath = joinpath(testdir, "testassets", "test_job")

job = DFJob(testjobpath);
job4 = DFJob(testjobpath)

#output stuff
out = outputdata(job;print=false,onlynew=false);
@test haskey(out, "nscf")
@test haskey(out["nscf"], :fermi)

@test isconverged(job["scf"])
nscf = DFControl.input(job, "nscf")
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


tcalc = gencalc_bands(job["bands"], [(0.5,0.0,0.5,10.0),(0.0,0.0,0.0,10.0),(0.5,0.5,0.5,1.0)], name="bands2")

@test gencalc_bands(job, [(0.5,0.0,0.5,10.0),(0.0,0.0,0.0,10.0),(0.5,0.5,0.5,1.0)], name = "bands2", template_name="bands") ==
    tcalc

push!(job, tcalc)
@test job["bands2"][:calculation] == "bands"
@test data(job, "bands2", :k_points).data == [(0.5,0.0,0.5,10.0),(0.0,0.0,0.0,10.0),(0.5,0.5,0.5,1.0)]

tcalc = gencalc_nscf(job["nscf"], (10,10,10), name="nscf2")
@test gencalc_nscf(job, (10,10,10), name="nscf2", template_name="nscf") == tcalc
push!(job, tcalc)
@test job["nscf2"][:calculation] == "nscf"
push!(job, gencalc_scf(job["scf"], (5,5,5,1,1,1), name="1scf2"))
@test job["1scf2"][:calculation] == "scf"


wanflags = [:write_hr => true, :wannier_plot => true]

setprojections!(job, :Pt => [:s, :p, :d])
wancalc = gencalc_wan(job["nscf"], structure(job), fermi-7.0, wanflags...; Epad=5.0)
@test wancalc == gencalc_wan(job, fermi-7.0, wanflags...; Epad=5.0)
@show DFControl.iscolin(job.structure)
append!(job, wancalc)
@test job["wan"][:write_hr] == job["wan"][:wannier_plot] == true

wanout = outputdata(job, "wan")
@test length(wanout[:final_state]) == 20
@test length(wanout[:wannierise]) == 203
job.inputs = job.inputs[1:4]

job[:nbnd] = 300
@test job["scf"][:nbnd] == 300
@test job["nscf"][:nbnd] == 300
rmflags!(job, :nbnd)
@test !haskey(job["scf"].flags, :nbnd)

job["scf"][:nbnd] = 300
@test job["scf"][:nbnd] == 300
rmflags!(job, :nbnd)

job["nscf"][:Hubbard_U] = [1.0]
@test job["nscf"][:Hubbard_U] == [1.0]

job["nscf"][:starting_magnetization] = [[1.0, 2.0, 1.0]]
@test job["nscf"][:starting_magnetization] == [[1.0, 2.0, 1.0]]

job["nscf"][:Hubbard_J] = [0 1 2]
@test job["nscf"][:Hubbard_J] == [0 1 2]
rmflags!(job, :Hubbard_J, :starting_magnetization, :Hubbard_U)
set_magnetization!(job, :Pt => [0.,0.,1.0])

push!.((job,), gencalc_wan(nscf, structure(job), fermi-7.0, wanflags..., Epad=5.0))
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

setpseudos!(job, :Pt => Pseudo("Pt.UPF", joinpath(testdir, "testassets/pseudos")))
@test job.structure.atoms[1].pseudo == Pseudo("Pt.UPF", joinpath(testdir, "testassets/pseudos"))
setpseudos!(job, :Pt, :test)
@test job.structure.atoms[1].pseudo == Pseudo("Pt.UPF", joinpath(testdir, "testassets/pseudos"))

setpseudos!(job, :test)
@test job.structure.atoms[1].pseudo == Pseudo("Pt.UPF", joinpath(testdir, "testassets/pseudos"))


testorbs = [:s, :p]
setprojections!(job, :Pt => testorbs)
@test convert.(Symbol, [p.orb for p in projections(job, :Pt)]) == testorbs
setwanenergies!(job, nscf, fermi-7.0, Epad=3.0)

@test job["wanup"][:dis_froz_max] == 13.2921
@test job["wanup"][:dis_win_max] == 16.292099999999998

setexecflags!(job, "pw.x", :nk => 230, :ndiag => 3)
@test DFControl.getfirst(x->x.symbol==:nk, execs(job, "nscf")[2].flags).value == 230
@test DFControl.getfirst(x->x.symbol==:ndiag, execs(job, "nscf")[2].flags).value == 3
rmexecflags!(job, "pw.x", :nk)
@test isempty(filter(x->x.symbol == :nk, execs(job, "nscf")[2].flags))

setexecdir!(job, "pw.x", joinpath(homedir(), "bin"))
@test execs(job, "nscf")[2].dir == joinpath(homedir(), "bin")

setname!(job, "nscf", "test")
@test DFControl.inpath(job, "test") == joinpath(job.local_dir, "nscf.in")
setname!(job, "test", "nscf")

setserverdir!(job, "localhost")
@test job.server_dir == "localhost"

setheaderword!(job, "defpart" => "frontend", print=false)
@test any(occursin.("frontend",job.header))

setdataoption!(job, "nscf",:k_points, :blabla, print=false)
@test data(job, "nscf", :k_points).option == :blabla
setdataoption!(job, :k_points, :test, print=false)
@test data(job, "nscf", :k_points).option == :test

report = progressreport(job; onlynew=false, print=false)
@test report[:fermi] == 17.4572
@test length(report[:accuracy]) == 9
oldat = atoms(job)[1]
newatompos = outputdata(job, "vc_relax", onlynew=false)[:final_structure]
job.structure = newatompos
push!(job.structure.atoms, oldat)

cp(joinpath(testdir, "testassets", "pseudos"), joinpath(testdir, "testassets", "pseudos_copy"), force=true)
setpseudos!(job, :Si => Pseudo("Si.UPF", joinpath(testdir, "testassets", "pseudos_copy")))
save(job)
@test ispath(joinpath(job.local_dir, "Si.UPF"))
@test atom(job, :Si).pseudo == Pseudo("Si.UPF", job.local_dir)
rm(joinpath(testdir, "testassets", "pseudos_copy"), recursive=true)

job["nscf"][:occupations] = "smearing"
job["nscf"][:degauss] = 2.0
job["nscf"][:smearing] = "mp"
projwfc = gencalc_projwfc(job["nscf"], -20, 10, 0.05)
@test projwfc[:Emin] == -20
@test projwfc[:degauss] == job["nscf"][:degauss]
@test projwfc[:ngauss] == 1

set_Hubbard_U!(job, :Si => 1.0)
@test atoms(job, :Si)[1].dftu.U == 1.0


prev_a = cell(job)[1, :]
prev_b = cell(job)[2, :]
prev_c = cell(job)[3, :]
prev_pos = position_cart.(atoms(job))
scale_cell!(job, 2)
@test prev_a .* 2 == cell(job)[1, :]
@test prev_b .* 2 == cell(job)[2, :]
@test prev_c .* 2 == cell(job)[3, :]
for (p, at) in zip(prev_pos, atoms(job))
	@test round.(DFControl.ustrip(p * 2), digits=5) == round.(DFControl.ustrip(position_cart(at)), digits=5)
end

set_magnetization!(job, :Pt => [1.0, 0.0, 0.0])
@test magnetization(atom(job, :Pt)) == DFControl.Vec3(1.0, 0.0, 0.0)


at = atoms(job)[1]
c = cell(job)
orig_pos_cart  = position_cart(at)
orig_pos_cryst = position_cryst(at)

set_position!(at, orig_pos_cart .+ c * [0.1, 0.1, 0.1], c)
@test isapprox(position_cryst(at), orig_pos_cryst .+ [0.1, 0.1, 0.1])

set_position!(at, orig_pos_cryst .+ [0.1, 0.1, 0.1], c)
@test position_cart(at) == c * position_cryst(at)

at2 = atoms(job)[2]
p1, p2 = position_cart(at), position_cart(at2)
mid = (p1+p2)/2
bondlength = distance(at, at2)
scale_bondlength!(at, at2, 0.5, c)

@test isapprox((position_cart(at) + position_cart(at2))/2, mid)
@test isapprox(distance(at, at2), bondlength/2) 

@test isapprox(bandgap(job), 2.6701999999999995)

t = job["nscf"]
curlen = length(job.inputs)
id = findfirst(x->x.name == "nscf", job.inputs)
n = pop!(job, "nscf")
@test n == t
@test length(job.inputs) == curlen - 1

nscf = job["scf"]
rmflags!(job, :nspin, :lda_plus_u, :noncolin)
set_magnetization!(job, :Pt => [0.2, 1.0, 0.2])
DFControl.sanitizeflags!(job)
set_Hubbard_U!(job, :Pt => 2.3)
@test job["scf"][:lda_plus_u]
@test job["scf"][:noncolin]
@test job["scf"][:lda_plus_u_kind] == 1

rmflags!(job, :nspin, :lda_plus_u, :noncolin)
set_magnetization!(job, :Pt => [0.0, 0.0, 0.5])
DFControl.sanitizeflags!(job)
@test job["scf"][:nspin] == 2

using LinearAlgebra
new_str = create_supercell(structure(job), 1, 0, 0)
prevcell = cell(job)
@test norm(cell(new_str)[:,1]) == norm(prevcell[:, 1]) * 2
prevlen = length(atoms(job))
@test length(atoms(new_str)) == 2*prevlen
prevlen_Pt = length(atoms(job, :Pt))

set_magnetization!(job, :Pt => [0, 0, 1])
orig_projs = projections(atoms(job, :Pt)[1])

job.structure = create_supercell(structure(job), 1, 0, 0, make_afm=true)

DFControl.sanitize_magnetization!(job)
@test length(atoms(job, :Pt)) == prevlen_Pt
@test length(atoms(job, :Pt1)) == prevlen_Pt

@test magnetization(atoms(job, :Pt)[1]) == [0, 0, 1]
@test magnetization(atoms(job, :Pt1)[1]) == [0, 0, -1]

DFControl.sanitize_projections!(job)
@test projections(atoms(job, :Pt)[1]) != projections(atoms(job, :Pt1)[1])

job4.server_dir = "/tmp"
save(job4)
job3 = DFJob(job4.local_dir)

@test job3.server_dir == "/tmp"

rm.(DFControl.inpath.(job.inputs))

# rm(joinpath(splitdir(DFControl.inpath(job.inputs[1]))[1], "pw2wan_wanup.in"))
# rm(joinpath(splitdir(DFControl.inpath(job.inputs[1]))[1], "pw2wan_wandn.in"))
rm(joinpath(job.local_dir, "job.tt"))
rm.(joinpath.((job.local_dir,), filter(x -> occursin("UPF", x), readdir(job.local_dir))))


