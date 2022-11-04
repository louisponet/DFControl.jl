using DFControl, Test

testassetspath = joinpath(testdir, "testassets")
testjobpath = joinpath(testassetspath, "test_job")

@testset "Structure manipulation" begin
    job = load(test_server, Job(testjobpath))

    @test length(Structures.symmetry_operators(job.structure)[1]) == 48
    
    job.structure = Structures.create_supercell(job.structure, 1, 0, 0, make_afm = true)
    @test length(Structures.symmetry_operators(job.structure)[1]) == 12
    @test length(job.structure) == 2 
    @test isapprox(job.structure[:Ni1][1].magnetization, [0,0,-0.1])
    @test isapprox(job.structure[:Ni][1].magnetization, [0,0,0.1])
    @test isapprox(map(x->x.magnetization, job.structure[Structures.element(:Ni)]), [[0,0,0.1], [0,0,-0.1]])
    ngl = Structures.niggli_reduce(job.structure)
    @test isapprox(Structures.volume(job.structure),  Structures.volume(ngl))
    @test job.structure.cell != ngl.cell
    prev = job.structure[Structures.element(:Ni)][1].position_cryst[1]
    prev_vol = Structures.volume(job.structure)
    Structures.scale_cell!(job.structure, diagm(0=>[2,1,1]))
    @test isapprox(job.structure[Structures.element(:Ni)][1].position_cryst[1],0.5prev)
    @test isapprox(2*prev_vol,  Structures.volume(job.structure))

    out_str = outputdata(job)["scf"][:initial_structure]
    Structures.update_geometry!(job.structure, out_str)
    for a in job.structure[element(:Ni)]
        a.dftu.l = 2
    end
    @test out_str == job.structure
end

@testset "supercell" begin
    job = load(test_server, Job(testjobpath))
    struct2 = Structures.create_supercell(job.structure, 1, 2, 1)
    newpositions = [at.position_cart for at in struct2.atoms]
    oldposition = job.structure.atoms[1].position_cart
    cell_ = job.structure.cell

    @test isapprox(job.structure.atoms[1].position_cart,
          job.structure.cell' * job.structure.atoms[1].position_cryst)
    @test isapprox(oldposition, newpositions[1])
    @test oldposition + cell_[:, 1] ∈ newpositions
    @test oldposition + 2 * cell_[:, 2] ∈ newpositions
    @test oldposition + 1 * cell_[:, 3] ∈ newpositions
    @test oldposition + (1 * cell_[:, 3] + 1 * cell_[:, 2]) ∈ newpositions
end

@testset "calculation management" begin
    job = load(test_server, Job(testjobpath))
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
    
end

@testset "versioning" begin
    job = load(test_server, Job(testjobpath))
    job[:nbnd] = 30
    curver = job.version
    save(job, fillexecs=false)
    @test job.version == curver + 1
    save(job, fillexecs=false)
    @test ispath(joinpath(job, Jobs.VERSION_DIR_NAME))
    @test ispath(joinpath(job, Jobs.VERSION_DIR_NAME, "$(curver+1)"))
    job[:nbnd] = 40
    save(job, fillexecs=false)
    @test job.version == curver + 3
    c = job["scf"]
    @test c[:nbnd] == 40
    switch_version!(job, curver + 1)
    @test job.version == curver + 1
    # @test DFControl.last_version(job) == 2
    c = job["scf"]
    @test c[:nbnd] == 30
    rm_version!(job, curver + 1)
    @test !ispath(joinpath(Jobs.main_job_dir(job), Jobs.VERSION_DIR_NAME,
                           "$(curver + 1)"))
    @test ispath(joinpath(Jobs.main_job_dir(job), Jobs.VERSION_DIR_NAME,
                          "$(curver + 2)"))
end
