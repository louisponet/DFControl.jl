using DFControl, Test

testdir = joinpath(dirname(dirname(pathof(DFControl))), "test")
testjobpath = joinpath(testdir, "testassets", "test_job")

name = "Pt"
local_dir = testjobpath
server_dir = ""
bin_dir = joinpath(homedir(), "bin")
excs = [Exec("mpirun", bin_dir, :np => 24), Exec("pw.x", bin_dir, :nk => 2)]

pseudoset = :test

header = ["#SBATCH -N 1", "#SBATCH --ntasks-per-node=24",
          "#SBATCH --time=24:00:00", "#SBATCH -p defpart",
          "module load open-mpi/gcc/1.10.4-hfi", "module load mkl/2016.1.056"
         ]

calculations = [DFCalculation{QE}(name="vc_relax", execs = excs,flags = Dict(:calculation => "vc-relax", :verbosity => "low"), data = [InputData(:k_points, :automatic, [6,6,6,1,1,1])]),
                DFCalculation{QE}(name="scf", execs = excs,   flags = Dict(:calculation => "scf", :verbosity => "low"), data = [InputData(:k_points, :automatic, [6,6,6,1,1,1])]),
                DFCalculation{QE}(name="bands", execs = excs, flags = Dict(:calculation => "bands", :verbosity => "high", :nbnd=>8), data = [InputData(:k_points,:crystal_b, [[0.5, 0.5, 0.5, 100.],
                                                                                                                                    [0.0, 0.0, 0.0, 100.],
                                                                                                                                    [0.0, 0.5, 0.0, 1.]])]),
                DFCalculation{QE}(name="nscf", execs = excs, flags = Dict(:calculation => "nscf", :verbosity => "low"))]
str = Structure(joinpath(testjobpath, "Pt.cif"), name="Pt")
set_pseudos!(str, :test)
job = DFJob(name, str, calculations,
      :prefix       => "$name",
      :restart_mode => "from_scratch",
      :mixing_mode  => "plain",
      :mixing_beta  => 0.7,
      :conv_thr     => 1.0e-8,
      #kwargs
      header      = header,
      local_dir = local_dir
     )
set_kpoints!(job["nscf"], (10,10,10))
save(job)
show(job)

@test data(job["scf"], :k_points).data == [6, 6, 6, 1, 1, 1]
@test data(job, "nscf", :k_points).data == DFControl.kgrid(10, 10, 10, :nscf)
@test all(values(job[:ecutwfc]) .== 32.0)
@test job["scf"][:prefix] == job["nscf"][:prefix] == "$name"
@test job["bands"][:verbosity] == "high"

set_flags!(job, :prefix => "blabla", print=false)
@test job["scf"][:prefix] == job["nscf"][:prefix] == "blabla"
set_flags!(job, :Hubbard_U => [4], print=false)
set_flags!(job, :Hubbard_J => [4 4 5], print=false)
@test job["scf"][:Hubbard_U] == job["nscf"][:Hubbard_U] == [4.0]
@test job["scf"][:Hubbard_J] == job["nscf"][:Hubbard_J] == [4.0  4.0  5.0]


struct2 = DFControl.create_supercell(job.structure, 1, 2, 1)
newpositions = [DFControl.position_cart(at) for at in atoms(struct2)]
oldposition = atoms(job.structure)[1].position_cart
cell_ = DFControl.cell(job.structure)

@test atoms(job.structure)[1].position_cart == job.structure.cell' * atoms(job.structure)[1].position_cryst
@test oldposition == newpositions[1]
@test oldposition + cell_[:, 1] ∈ newpositions
@test oldposition + 2*cell_[:, 2] ∈ newpositions
@test oldposition + 1*cell_[:, 3] ∈ newpositions
@test oldposition + (1*cell_[:, 3] + 1*cell_[:, 2]) ∈ newpositions
