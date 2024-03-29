# # Basic Tutorial
# !!! note
#     Make sure to first go through the [`Configuration`](@ref Configuration) steps.
# Since DFControl is aimed at improving the day to day quality of life of a material's scientist/anyone running DFT codes,
# we will look at a simple demonstration of how by creating and submitting some
# [Quantum-Espresso](https://www.quantum-espresso.org/) calculations on Si
# starting from a cif file specifying the structure.

using DFControl

# First we download the cif file, extract the `Structure` and assign the right pseudos to it.
# In this case Si (F d -3 m :1) from http://www.crystallography.net/cod/9011998.cif
using Downloads
cif_file = Downloads.download("http://www.crystallography.net/cod/9011998.cif", "Si.cif")

structure = Structure(cif_file)
if false#hide
set_pseudos!(structure, load(local_server(), PseudoSet("pbesol")))
end#hide

# This assumes that the `"pbesol"` pseudopotential set was installed during the
#md # [configuration step](@ref Configuration).

# Next we specify the executables with which to run the QE calculations.
# We assume here that `mpirun` is installed and in the user's PATH, and that
# QE is installed, and to be found in the `/opt/qe/bin`, change this according to your own
# setup. The first argument to the constructor can be used as a label to later retrieve the executable after it was saved.

pw_exec = Exec(name="pw", path="/opt/qe/bin/pw.x",flags=Dict("-nk" => 4))

# Then we create the first calculation for our job, we name it scf, which will be used to reference it later.
# We also pass the executables to be used and additional flags to be set to the constructor.
# Afterwards we set the kpoints to be used in the scf calculation.
scf_calculation = Calculation("scf", :calculation => "scf"; exec = pw_exec)
set_kpoints!(scf_calculation, (6, 6, 6, 1, 1, 1))

# The code recognizes internally that this 6-Tuple corresponds to a
# `K_POINTS (automatic)` block in QE. Alternatively (leading to an identical final result):

scf_calculation = Calculation("scf", :calculation => "scf"; exec = pw_exec,
                                  data = [InputData(:k_points, :automatic,
                                                    (6, 6, 6, 1, 1, 1))])
# We can now define our job:  job = Job("Si", structure, [scf_calculation], :ecutwfc => 20, :conv_thr => 1e-6; dir="job")
# Additional calculations would be be added to the list `[scf_calculation]`.
# The flag => value pairs will set the specified flags to that value for all calculations in the job
# that allow recognize that flag, so it's ideal for things like cutoffs and smearing etc.

if false#hide
job = Job("Si", structure, [scf_calculation], :ecutwfc => 40.0, :occupations => "smearing", :degauss=>0.01, :conv_thr => 1e-6, :nbnd => 18;
            dir = dir, server=gethostname(), environment="default")
end#hide

# We are now ready to submit the job, which will run in the current working directory
if false #hide
submit(job)
else #hide
    global job = load(local_server(), Job(joinpath(splitdir(pathof(DFControl))[1], "..", "docs","src","assets", "job")))#hide
    pop!(job) #hide
end #hide

# This will generate and save all calculation files, and the corresponding job script (`job.sh`),
# to the server specified in `job.server`, and then the job will be submitted on the server.

# After the job finishes the outputs can be parsed through
outputdata(job)
# this returns a dictionary with the input names as keys and a Dict with parsed outputs as values.
# i.e.
outputdata(job)["scf"]
# takes the outputdata of the input that we previously named "scf" when creating it.

# Now that the scf calculation finished succesfully, the next step is usually to
# have a look at the bandstructure. For this we generate a bands calculation,
# using the scf calculation as a template and a generated high symmetry k-point path
# with 20 kpoints per segment.

bands_calc = Calculations.gencalc_bands(job["scf"], Structures.high_symmetry_kpath(job.structure, 20))

# Observe the `:calculation => "bands"`, and automatic setting of the `:verbosity => "high"` flags.
# We now `push!` this calculation to the job queue.
push!(job, bands_calc)

# However, since we know the scf succeeded there is no need to rerun it.
# To un-schedule it we do
job["scf"].run = false

# Printing the job will now highlight the scheduled calculations differently from the non-scheduled ones
job

# Seeing that all is right we submit the job again
job.dir = "job"; #hide
if false #hide
submit(job)
else #hide
    global job = load(local_server(),Job(joinpath(splitdir(pathof(DFControl))[1], "..", "docs","src","assets", "job")));#hide
end #hide

# We can access the bands through
bands = readbands(job);
# or
bands = outputdata(job)["bands"][:bands];

# They can be plotted too
fermi = readfermi(job) # = outputdata(job)["scf"][:fermi]
using Plots
plot(bands; fermi = fermi)

# Since more info (such as the structure) is available in the job,
# plotting the job leads to a richer plot
plot(job, -10, 1)

# As can be seen in the [Advanced Usage](@ref), additional information
# generated by additional calculations will be picked up by DFControl
# in order to create richer plots.
