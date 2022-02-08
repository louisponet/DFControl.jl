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
    pw_exec = Exec("pw", "pw.x", bin_dir, :nk => 4)
    @test !isempty(pw_exec.flags)
    pseudoset = :test

    str = Structure(joinpath(testdir, "testassets/Ni.cif"))

    calculations = [Calculation("vcrelax", :calculation => "vc-relax", :verbosity => "high", :ion_dynamics => "bfgs", :cell_dynamics => "bfgs";
                                      exec = pw_exec,
                                      data = [InputData(:k_points, :automatic,
                                                        [6, 6, 6, 1, 1, 1])]),
                    Calculation(; name = "scf", exec = pw_exec,
                                      flags = Dict(:calculation => "scf", :verbosity => "high"),
                                      data = [InputData(:k_points, :automatic,
                                                        [4, 4, 4, 1, 1, 1])])]
    job = Job(name, str, calculations, :ecutwfc => 40.0, :occupations => "smearing", :degauss=>0.01, :conv_thr => 1e-6, :nbnd => 18;
                #kwargs
                dir = dir, server="localhost")


    set_pseudos!(job, :test)
    job.environment = "test_default"

    set_kpoints!(job["scf"], (6, 6, 6, 1, 1, 1))

    job.structure[element(:Ni)][1].magnetization = Vec3(0,0,0.1)
    job.structure[element(:Ni)][1].dftu.U = 4.0

    push!(job, Calculations.gencalc_bands(job["scf"], Structures.high_symmetry_kpath(job.structure, 20)))
    push!(job, Calculations.gencalc_nscf(job["scf"], (5,5,5)))

    push!(job, Calculations.gencalc_projwfc(job["nscf"], 2.0, 20.0, 0.1))

    @test all(values(job[:ecutwfc]) .== 40.0)
    for c in job.calculations
        pop!(c, :ecutwfc, nothing)
    end

    
    save(job)
    @test all(values(job[:ecutwfc]) .== 41.0)
    @test all(values(job[:ecutrho]) .== 236.0)
    @test job.version == 1
    @test length(job) == 5
    @test data(job["scf"], :k_points).data == [6,6,6,1,1,1]
    @test job["nscf"].exec == pw_exec
    @test job["projwfc"].exec.exec == "projwfc.x"
    @test job["projwfc"].exec.dir == pw_exec.dir
    @test show(job) == nothing
    job[:ecutwfc] = 40.0
    for c in job.calculations
        pop!(c, :ecutrho, nothing)
    end
    save(job)

    job2 = load(test_server, Job(abspath(job)))
    for (c1, c2) in zip(job2.calculations, job.calculations)
        @test c1 == c2
    end
    @test job2.structure == job.structure
end

refjobpath =joinpath(testdir, "testassets", "reference_job")

@testset "reference comparison" begin
    job = load(test_server, Job(testjobpath))
    
    orig_job = deepcopy(job)
    job.structure = Structures.create_supercell(job.structure, 1, 0, 0, make_afm = true)
    
    job2 = load(test_server, Job(refjobpath))
    @test job2.structure == job.structure
    
    for f in DFControl.Utils.searchdir(job2, ".out")
        cp(f, joinpath(job, splitdir(f)[2]), force=true)
    end
    for f in DFControl.Utils.searchdir(job2, "dos")
        cp(f, joinpath(job, splitdir(f)[2]), force=true)
    end

    for a in job.structure.atoms
        a.projections = [Projection("s"), Projection("p"), Projection("d")]
    end
    wanexec = Exec("wan","wannier90.x", joinpath(homedir(), "Software/wannier90"))
    append!(job, Calculations.gencalc_wan(job, 0.000011, wanexec = wanexec))
    for (c1, c2) in zip(job2.calculations, job.calculations)
        @test c1 == c2
    end
    save(job)
    job = load(test_server, Job(testjobpath))
    
    for (c1, c2) in zip(job2.calculations, job.calculations)
        @test c1 == c2
    end
    save(orig_job)
    for f in DFControl.Utils.searchdir(job2, ".out")
        cp(f, joinpath(job, splitdir(f)[2]), force=true)
    end
    for f in DFControl.Utils.searchdir(job2, "dos")
        cp(f, joinpath(job, splitdir(f)[2]), force=true)
    end
end
