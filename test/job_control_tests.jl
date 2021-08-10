using DFControl, Test

testassetspath = joinpath(testdir, "testassets")
testjobpath = joinpath(testassetspath, "test_job")

@testset "Structure manipulation" begin
    job = Job(testjobpath, "localhost_test")

    @test length(Structures.symmetry_operators(job.structure)[1]) == 48
    
    job.structure = Structures.create_supercell(job.structure, 1, 0, 0, make_afm = true)
    @test length(Structures.symmetry_operators(job.structure)[1]) == 12
    @test length(job.structure) == 2 
    @test isapprox(job.structure[:Ni1][1].magnetization, [0,0,-0.1])
    @test isapprox(job.structure[:Ni][1].magnetization, [0,0,0.1])
    @test isapprox(map(x->x.magnetization, job.structure[Structures.element(:Ni)]), [[0,0,0.1], [0,0,-0.1]])
    ngl = Structures.niggli_reduce(job.structure)
    @test Structures.volume(job.structure) == Structures.volume(ngl)
    @test job.structure.cell != ngl.cell
    prev = job.structure[Structures.element(:Ni)][1].position_cryst[1]
    prev_vol = Structures.volume(job.structure)
    Structures.scale_cell!(job.structure, diagm(0=>[2,1,1]))
    @test job.structure[Structures.element(:Ni)][1].position_cryst[1] == 0.5prev
    @test 2*prev_vol == Structures.volume(job.structure)

    out_str = outputdata(job)["scf"][:initial_structure]
    Structures.update_geometry!(job.structure, out_str)
    @test out_str == job.structure
end

@testset "supercell" begin
    job = Job(testjobpath, "localhost_test")
    struct2 = Structures.create_supercell(job.structure, 1, 2, 1)
    newpositions = [at.position_cart for at in struct2.atoms]
    oldposition = job.structure.atoms[1].position_cart
    cell_ = job.structure.cell

    @test job.structure.atoms[1].position_cart ==
          job.structure.cell' * job.structure.atoms[1].position_cryst
    @test oldposition == newpositions[1]
    @test oldposition + cell_[:, 1] ∈ newpositions
    @test oldposition + 2 * cell_[:, 2] ∈ newpositions
    @test oldposition + 1 * cell_[:, 3] ∈ newpositions
    @test oldposition + (1 * cell_[:, 3] + 1 * cell_[:, 2]) ∈ newpositions
end

@testset "registry" begin
    job = Job(testjobpath,  "localhost_test")
    rm(testjobpath, recursive=true)
    prevlen = length(Jobs.registered_jobs())
    save(job)
    @test length(Jobs.registered_jobs()) == prevlen + 1
    @test Jobs.registered_jobs()[end][1] == job.dir
end

@testset "execs" begin
    job = Job(testjobpath,  "localhost_test")
    for c in job.calculations
        for e in c.execs
            if e.exec == "pw.x"
                Calculations.set_flags!(e, :nk => 4)
                e.dir = "test/test"
            elseif e.exec == "mpirun"
                Calculations.set_flags!(e, :np => 4)
                @test_throws ErrorException Calculations.set_flags!(e, :nasdf => 4)
            end
        end
    end
    ex = job.calculations[1].execs[2]
    mpi = job.calculations[1].execs[1]
    @test ex.dir == "test/test"
    @test ex.flags[1].symbol == :nk
    @test ex.flags[1].value == 4
    @test mpi.flags[1].symbol == :np
    @test mpi.flags[1].value == 4

    save(job)
    job = Job(testjobpath, "localhost_test")
    ex = job.calculations[1].execs[2]
    mpi = job.calculations[1].execs[1]
    @test ex.flags[1].symbol == :nk
    @test ex.flags[1].value == 4
    @test mpi.flags[1].symbol == :np
    @test mpi.flags[1].value == 4

    for c in job.calculations
        for e in c.execs
            if e.exec == "pw.x"
                Calculations.rm_flags!(e, :nk)
                e.dir = "test/test"
            elseif e.exec == "mpirun"
                Calculations.rm_flags!(e, :np)
            end
        end
    end
   
    @test isempty(ex.flags)
    @test isempty(mpi.flags)
    save(job)
end

@testset "calculation management" begin
    job = Job(testjobpath, "localhost_test")
    ncalcs = length(job.calculations)
    t = pop!(job, "scf")
    @test length(job.calculations) == ncalcs - 1
    @test job[2].name == "bands"
    insert!(job, 2, t)
    @test job.calculations[2].name == "scf"

    Jobs.set_flow!(job, "" => false, "scf" => true)
    @test job["scf"].run == true
    @test job["nscf"].run == true
    @test job["bands"].run == false
    push!(job.header, "module load testmod")
    Jobs.set_headerword!(job, "testmod" => "testmod2")
    @test job.header[end] == "module load testmod2"
    
end

@testset "versioning" begin
    job = Job(testjobpath, "localhost_test")
    job[:nbnd] = 30
    curver = Client.last_version(job)
    save(job)
    @test job.version == curver + 1
    save(job)
    @test ispath(joinpath(job, Jobs.VERSION_DIR_NAME))
    @test ispath(joinpath(job, Jobs.VERSION_DIR_NAME, "$(curver+1)"))
    job[:nbnd] = 40
    save(job)
    @test job.version == curver + 3
    @test job["scf"][:nbnd] == 40
    switch_version!(job, curver + 1)
    @test job.version == curver + 1
    # @test DFControl.last_version(job) == 2
    @test job["scf"][:nbnd] == 30
    rm_version!(job, curver + 1)
    @test !ispath(joinpath(Jobs.main_job_dir(job), Jobs.VERSION_DIR_NAME,
                           "$(curver + 1)"))
    @test ispath(joinpath(Jobs.main_job_dir(job), Jobs.VERSION_DIR_NAME,
                          "$(curver + 2)"))
end


rm(testjobpath, recursive=true)
# function copy_outfiles()
#     for f in readdir(joinpath(testassetspath, "outputs"))
#         cp(joinpath(testassetspath, "outputs", f), joinpath(testjobpath, f); force = true)
#     end
# end

# copy_outfiles()

# job = Job(testjobpath);
# job4 = Job(testjobpath)

# #output stuff
# out = outputdata(job; print = false, onlynew = false);
# @test haskey(out, "nscf")
# @test haskey(out["nscf"], :fermi)

# @test isconverged(job["scf"])
# nscf = DFControl.calculation(job, "nscf")
# nscf2 = Calculation(nscf, "nscf2";
#                       data = [:testdata => (:testoption, "test"),
#                               :k_points => (:blabla, [1, 1, 1, 1, 1, 1])])

# @test data(nscf2, :testdata).option == :testoption
# @test data(nscf2, :testdata).data == "test"
# @test data(nscf2, :testdata).name == :testdata
# @test data(nscf2, :k_points).option == :blabla
# @test data(nscf2, :k_points).data == [1, 1, 1, 1, 1, 1]
# @test data(job["scf"], :k_points).option == :automatic
# @test_throws ErrorException job["bladbsflk"]
# @test_throws ErrorException nscf2[:bladkfj]

# set_kpoints!(nscf2, (3, 3, 3); print = false)
# @test data(nscf2, :k_points).data == DFControl.kgrid(3, 3, 3, :nscf)

# set_kpoints!(nscf2, [(3.0, 3.0, 3.0, 1.0), (3.0, 3.0, 3.0, 1.0)]; print = false)
# @test data(nscf2, :k_points).option == :crystal_b
# @test data(nscf2, :k_points).data == [(3.0, 3.0, 3.0, 1.0), (3.0, 3.0, 3.0, 1.0)]

# set_kpoints!(nscf2, (3, 3, 3, 0, 0, 1); print = false)
# @test data(nscf2, :k_points).option == :automatic
# @test data(nscf2, :k_points).data == [3, 3, 3, 0, 0, 1]

# fermi = outputdata(job, "nscf")[:fermi]
# @test fermi == 17.4572

# tcalc = gencalc_bands(job["bands"],
#                       [(0.5, 0.0, 0.5, 10.0), (0.0, 0.0, 0.0, 10.0), (0.5, 0.5, 0.5, 1.0)];
#                       name = "bands2")

# push!(job, tcalc)
# @test job["bands2"][:calculation] == "bands"
# @test data(job["bands2"], :k_points).data ==
#       [(0.5, 0.0, 0.5, 10.0), (0.0, 0.0, 0.0, 10.0), (0.5, 0.5, 0.5, 1.0)]

# tcalc = gencalc_nscf(job["nscf"], (10, 10, 10); name = "nscf2")
# push!(job, tcalc)
# @test job["nscf2"][:calculation] == "nscf"
# push!(job, gencalc_scf(job["scf"], (5, 5, 5, 1, 1, 1); name = "1scf2"))
# @test job["1scf2"][:calculation] == "scf"

# wanflags = [:write_hr => true, :wannier_plot => true]

# set_projections!(job, :Pt => ["s", "p", "d"])
# wancalc = gencalc_wan(job["nscf"], structure(job), fermi - 7.0, wanflags...; Epad = 5.0)
# append!(job, wancalc)
# @test job["wan"][:write_hr] == job["wan"][:wannier_plot] == true

# copy_outfiles()
# wanout = outputdata(job, "wan")
# @test length(wanout[:final_state]) == 20
# @test length(wanout[:wannierise]) == 203
# job.calculations = job.calculations[1:4]

# job[:nbnd] = 300
# @test job["scf"][:nbnd] == 300
# @test job["nscf"][:nbnd] == 300
# rm_flags!(job, :nbnd)
# @test !haskey(job["scf"].flags, :nbnd)

# job["scf"][:nbnd] = 300
# @test job["scf"][:nbnd] == 300
# rm_flags!(job, :nbnd)

# job["nscf"][:Hubbard_U] = [1.0]
# @test job["nscf"][:Hubbard_U] == [1.0]

# job["nscf"][:starting_magnetization] = [[1.0, 2.0, 1.0]]
# @test job["nscf"][:starting_magnetization] == [[1.0, 2.0, 1.0]]

# job["nscf"][:Hubbard_J] = [0 1 2]
# @test job["nscf"][:Hubbard_J] == [0 1 2]
# rm_flags!(job, :Hubbard_J, :starting_magnetization, :Hubbard_U)
# set_magnetization!(job, :Pt => [0.0, 0.0, 1.0])

# push!.((job,), gencalc_wan(nscf, structure(job), fermi - 7.0, wanflags...; Epad = 5.0))
# #TODO: add test with the new wancalc from projwfc

# set_flow!(job, "" => false)
# @test job.calculations[1].run == false
# set_flow!(job, "nscf" => true, "bands" => true)
# @test job.calculations[3].run

# save(job)
# job2 = Job(dir)
# begin
#     for (calc, calc2) in zip(job.calculations, job2.calculations)
#         for (f, v) in calc.flags
#             if f in (:Hubbard_J, :pseudo_dir, :wannier_plot)
#                 continue
#             else
#                 @test v == calc2.flags[f]
#             end
#         end
#         for (b1, b2) in zip(calc.data, calc2.data)
#             for n in fieldnames(typeof(b1))
#                 @test getfield(b1, n) == getfield(b2, n)
#             end
#         end
#     end
# end
# copy_outfiles()
# set_cutoffs!(job)
# @test job["scf"][:ecutwfc] == 32.0

# set_pseudos!(job, :Pt => Pseudo("Pt.UPF", joinpath(testdir, "testassets/pseudos")))
# @test job.structure.atoms[1].pseudo ==
#       Pseudo("Pt.UPF", joinpath(testdir, "testassets/pseudos"))
# set_pseudos!(job, :Pt, :test)
# @test job.structure.atoms[1].pseudo ==
#       Pseudo("Pt.UPF", joinpath(testdir, "testassets/pseudos"))

# set_pseudos!(job, :test)
# @test job.structure.atoms[1].pseudo ==
#       Pseudo("Pt.UPF", joinpath(testdir, "testassets/pseudos"))

# testorbs = ["s", "p"]
# set_projections!(job, :Pt => testorbs)
# @test convert.(String, [p.orb for p in Iterators.flatten(projections.(job.structure.atoms :Pt)))]) ==
#       testorbs
# set_wanenergies!(job, nscf, fermi - 7.0; Epad = 3.0)

# @test job["wanup"][:dis_froz_max] == 13.2921
# @test job["wanup"][:dis_win_max] == 16.292099999999998

# set_execflags!(job, "pw.x", :nk => 230, :ndiag => 3)
# @test DFControl.getfirst(x -> x.symbol == :nk, execs(job["nscf"])[2].flags).value == 230
# @test DFControl.getfirst(x -> x.symbol == :ndiag, execs(job["nscf"])[2].flags).value == 3
# rm_execflags!(job, "pw.x", :nk)
# @test isempty(filter(x -> x.symbol == :nk, execs(job["nscf"])[2].flags))

# set_execdir!(job, "pw.x", joinpath(homedir(), "bin"))
# @test execs(job["nscf"])[2].dir == joinpath(homedir(), "bin")

# set_name!(job["nscf"], "test")
# @test DFControl.inpath(job, "test") == joinpath(job.dir, "test.in")
# set_name!(job["test"], "nscf")

# set_serverdir!(job, "localhost")
# @test job.server_dir == "localhost"

# set_headerword!(job, "defpart" => "frontend"; print = false)
# @test any(occursin.("frontend", job.header))

# report = progressreport(job; onlynew = false, print = false)
# @test report[:fermi] == 17.4572
# @test length(report[:accuracy]) == 9
# oldat = job.structure.atoms[1]
# copy_outfiles()
# newatompos = outputdata(job, "vc_relax"; onlynew = false)[:final_structure]
# job.structure = newatompos
# push!(job.structure.atoms, oldat)

# cp(joinpath(testdir, "testassets", "pseudos"),
#    joinpath(testdir, "testassets", "pseudos_copy"); force = true)
# set_pseudos!(job, :Si => Pseudo("Si.UPF", joinpath(testdir, "testassets", "pseudos_copy")))
# save(job)
# @test ispath(joinpath(job.dir, "Si.UPF"))
# @test job.structure.atoms :Si)[1].pseudo == Pseudo("Si.UPF", job.dir)
# rm(joinpath(testdir, "testassets", "pseudos_copy"); recursive = true)

# job["nscf"][:occupations] = "smearing"
# job["nscf"][:degauss] = 2.0
# job["nscf"][:smearing] = "mp"
# projwfc = gencalc_projwfc(job["nscf"], -20, 10, 0.05)
# @test projwfc[:Emin] == -20
# @test projwfc[:degauss] == job["nscf"][:degauss]
# @test projwfc[:ngauss] == 1

# set_Hubbard_U!(job, :Si => 1.0)
# @test job.structure.atoms :Si)[1].dftu.U == 1.0

# prev_a = job.structure.cell[1, :]
# prev_b = job.structure.cell[2, :]
# prev_c = job.structure.cell[3, :]
# prev_pos = position_cart.(job.structure.atoms)
# scale_cell!(job, [2 0 0; 0 2 0; 0 0 2])
# @test prev_a .* 2 == job.structure.cell[1, :]
# @test prev_b .* 2 == job.structure.cell[2, :]
# @test prev_c .* 2 == job.structure.cell[3, :]
# for (p, at) in zip(prev_pos, job.structure.atoms)
#     @test round.(DFControl.ustrip.(p * 2), digits = 3) ==
#           round.(DFControl.ustrip.(at.position_cart), digits = 3)
# end

# set_magnetization!(job, :Pt => [1.0, 0.0, 0.0])
# @test magnetization(job.structure.atoms :Pt)[1]) == DFControl.Vec3(1.0, 0.0, 0.0)

# at             = job.structure.atoms[1]
# c              = job.structure.cell
# orig_pos_cart  = at.position_cart
# orig_pos_cryst = at.position_cryst

# set_position!(at, orig_pos_cart .+ c * [0.1, 0.1, 0.1], c)
# @test isapprox(at.position_cryst, orig_pos_cryst .+ [0.1, 0.1, 0.1])

# set_position!(at, orig_pos_cryst .+ [0.1, 0.1, 0.1], c)
# @test at.position_cart == c * at.position_cryst

# at2 = job.structure.atoms[2]
# p1, p2 = at.position_cart, at.position_cart)
# mid = (p1 + p2) / 2
# bondlength = distance(at, at2)
# scale_bondlength!(at, at2, 0.5, c)

# @test isapprox((at.position_cart + at.position_cart)) / 2, mid)
# @test isapprox(distance(at, at2), bondlength / 2)
# copy_outfiles()
# @test isapprox(bandgap(job), 2.6701999999999995)

# t = job["nscf"]
# curlen = length(job.calculations)
# id = findfirst(x -> x.name == "nscf", job.calculations)
# n = pop!(job, "nscf")
# @test n == t
# @test length(job.calculations) == curlen - 1

# nscf = job["scf"]
# rm_flags!(job, :nspin, :lda_plus_u, :noncolin)
# set_magnetization!(job, :Pt => [0.2, 1.0, 0.2])
# DFControl.sanitize_flags!(job)
# set_Hubbard_U!(job, :Pt => 2.3)
# DFControl.sanitize_magnetization!(job)
# DFControl.set_hubbard_flags!.(filter(x -> eltype(x) == QE,
#                                      job.calculations), (job.structure,))
# DFControl.set_starting_magnetization_flags!.(filter(x -> eltype(x) == QE,
#                                                     job.calculations),
#                                              (job.structure,))
# @test job["scf"][:lda_plus_u]
# @test job["scf"][:noncolin]
# @test job["scf"][:lda_plus_u_kind] == 1

# rm_flags!(job, :nspin, :lda_plus_u, :noncolin)
# set_magnetization!(job, :Pt => [0.0, 0.0, 0.5])
# DFControl.sanitize_flags!(job)
# DFControl.sanitize_magnetization!(job)
# DFControl.set_hubbard_flags!.(filter(x -> eltype(x) == QE,
#                                      job.calculations), (job.structure,))
# DFControl.set_starting_magnetization_flags!.(filter(x -> eltype(x) == QE,
#                                                     job.calculations),
#                                              (job.structure,))
# @test job["scf"][:nspin] == 2

# using LinearAlgebra
# new_str = create_supercell(structure(job), 1, 0, 0)
# prevcell = job.structure.cell
# @test norm(new_str.cell[:, 1]) == norm(prevcell[:, 1]) * 2
# prevlen = length(job.structure.atoms)
# @test length(new_str.atoms) == 2 * prevlen
# prevlen_Pt = length(job.structure.atoms :Pt))

# set_magnetization!(job, :Pt => [0, 0, 1])
# orig_projs = projections(job.structure.atoms :Pt)[1])

# job.structure = create_supercell(structure(job), 1, 0, 0; make_afm = true)
# DFControl.sanitize_magnetization!(job)
# DFControl.set_starting_magnetization_flags!.(filter(x -> eltype(x) == QE,
#                                                     job.calculations),
#                                              (job.structure,))
# @test length(job.structure.atoms :Pt)) == prevlen_Pt
# @test length(job.structure.atoms :Pt1)) == prevlen_Pt

# @test magnetization(job.structure.atoms :Pt)[1]) == [0, 0, 1]
# @test magnetization(job.structure.atoms :Pt1)[1]) == [0, 0, -1]

# DFControl.sanitize_projections!(job)
# @test projections(job.structure.atoms :Pt)[1]) != projections(job.structure.atoms :Pt1)[1])

# # job4.server_dir = "/tmp"
# # save(job4)
# # job3 = Job(job4.dir)

# # @test job3.server_dir == "/tmp"
# set_data_option!(job["scf"], :k_points, :blabla; print = false)
# @test data(job["scf"], :k_points).option == :blabla
# set_data_option!(job, :k_points, :test; print = false)
# @test data(job["scf"], :k_points).option == :test

# rm.(DFControl.inpath.(job.calculations))
# job.calculations = [job.calculations[2]]
# set_kpoints!(job["scf"], (6, 6, 6, 1, 1, 1))
# rm(joinpath(job, DFControl.VERSION_DIR_NAME); recursive = true)

# rm(joinpath(DFControl.main_job_dir(job), DFControl.VERSION_DIR_NAME); recursive = true)
# set_dir!(job, DFControl.main_job_dir(job))
# rm.(DFControl.inpath.(job.calculations))

# rm(joinpath(job, "job.tt"))
# rm(joinpath(job, ".metadata.jld2"))
# rm(joinpath(job, "pw2wan_wandn.in"))
# rm(joinpath(job, "pw2wan_wanup.in"))
# rm.(joinpath.((job.dir,), filter(x -> occursin("UPF", x), readdir(job.dir))))
