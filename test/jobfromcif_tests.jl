using DFControl, Test, LinearAlgebra

testdir = joinpath(dirname(dirname(pathof(DFControl))), "test")
testjobpath = joinpath(testdir, "testassets", "test_job")
@testset "initial creation" begin
    if ispath(testjobpath)
        rm(testjobpath, recursive=true)
    end

    name = "Ni"
    dir = testjobpath
    bin_dir = joinpath(homedir(), "Software/qe/bin")
    pw_excs = [Exec("mpirun", "", :np => 4), Exec("pw.x", bin_dir, :nk => 4)]

    pseudoset = :test

    header = ["#SBATCH -N 1"]

    str = Structure(joinpath(testdir, "testassets/Ni.cif"); name = "Ni")

    calculations = [DFCalculation{QE}("vcrelax", :calculation => "vc-relax", :verbosity => "high", :ion_dynamics => "bfgs", :cell_dynamics => "bfgs";
                                      execs = pw_excs,
                                      data = [InputData(:k_points, :automatic,
                                                        [6, 6, 6, 1, 1, 1])]),
                    DFCalculation{QE}(; name = "scf", execs = pw_excs,
                                      flags = Dict(:calculation => "scf", :verbosity => "high"),
                                      data = [InputData(:k_points, :automatic,
                                                        [4, 4, 4, 1, 1, 1])])]
    job = DFJob(name, str, calculations, :ecutwfc => 40.0, :occupations => "smearing", :degauss=>0.01, :conv_thr => 1e-6, :nbnd => 18;
                #kwargs
                header = header, dir = dir)

    set_pseudos!(job, :test)


    set_kpoints!(job["scf"], (6, 6, 6, 1, 1, 1))

    set_magnetization!(atoms(job, element(:Ni))[1], [0,0, 0.1])
    set_Hubbard_U!(job, element(:Ni) => 4.0)

    push!(job, gencalc_bands(job["scf"], high_symmetry_kpath(job, 20)))
    push!(job, gencalc_nscf(job["scf"], (5,5,5)))

    push!(job, gencalc_projwfc(job["nscf"], 2.0, 20.0, 0.1))

    save(job)
    @test job.version == 1
    @test length(job) == 5
    @test data(job["scf"], :k_points).data == [6,6,6,1,1,1]
    @test job["nscf"].execs == pw_excs
    @test job["projwfc"].execs == [pw_excs[1], Exec("projwfc.x", pw_excs[2].dir)]
    @test show(job) == nothing

    job2 = DFJob(job.dir)
    for (c1, c2) in zip(job2.calculations, job.calculations)
        @test c2 == c1
    end
    @test job2.structure == job.structure
    @test all(values(job[:ecutwfc]) .== 40.0)
    @test DFControl.find_cutoffs(job) == (41.0, 236.0)
end

refjobpath =joinpath(testdir, "testassets", "reference_job")

@testset "reference comparison" begin
    job = DFJob(testjobpath)
    orig_job = deepcopy(job)
    job.structure = create_supercell(job, 1, 0, 0, make_afm = true)
    
    job2 = DFJob(refjobpath)
    @test job2.structure == job.structure
    
    for f in DFControl.searchdir(job2, ".out")
        cp(f, joinpath(job, splitdir(f)[2]), force=true)
    end
    for f in DFControl.searchdir(job2, "dos")
        cp(f, joinpath(job, splitdir(f)[2]), force=true)
    end

    set_projections!(job, element(:Ni) => ["s", "p", "d"])
    wanexec = Exec("wannier90.x", joinpath(homedir(), "Software/wannier90"))
    append!(job, gencalc_wan(job, 0.000011, wanexec = wanexec))
    for (c1, c2) in zip(job2.calculations, job.calculations)
        @test c2 == c1
    end
    save(job)
    @test !ispath(joinpath(job, "scf.out"))
    job = DFJob(testjobpath)
    
    for (c1, c2) in zip(job2.calculations, job.calculations)
        @test c2 == c1
    end
    save(orig_job)
    for f in DFControl.searchdir(job2, ".out")
        cp(f, joinpath(job, splitdir(f)[2]), force=true)
    end
    for f in DFControl.searchdir(job2, "dos")
        cp(f, joinpath(job, splitdir(f)[2]), force=true)
    end
end
