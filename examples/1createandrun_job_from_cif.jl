using DFControl

#This example goes through how one would create a job from scratch, using a .cif file for the Structure.

#First go to your favourite Crystal structure database and download the .cif you want to use.
#e.g. Si (F d -3 m :1) : http://www.crystallography.net/cod/9011998.cif
cif_file = Downloads.download("http://www.crystallography.net/cod/9011998.cif", "Si.cif")

# Construct your structure from it, and set the pseudos to one of the sets we defined in the previous example
str = Structure(cif_file; name = "Si")
set_pseudos!(str, :pbesol)

# Since we are going to run QuantumEspresso, we define the executables to be used
pw_execs = [Exec("mpirun", "", :np => 4), Exec("pw.x", "/opt/qe/bin/", :nk => 4)]

# Then we generate the first calculation for our job, we name it scf, pass the executables to be used
# and set some specific flags.
# Afterwards the kpoints can also be set
scf_calculation = Calculation{QE}("scf", pw_execs, :calculation => "scf")
set_kpoints!(scf_calculation, (6, 6, 6, 1, 1, 1))
# Or

scf_calculation = Calculation{QE}("scf", pw_execs, :calculation => "scf";
                                    data = [InputData(:k_points, :automatic,
                                                      (6, 6, 6, 1, 1, 1))])

# Using these we can now define our job, if we would have more calculations they would be added
# to the list [scf_calculation].
# The flag => value pairs will set the specified flags to that value for all calculations in the job
# that allow that flag to be set, so it's ideal for things like cutoffs and smearing etc.
# The header can contain any lines that will be pasted in the job script in front of
# the calculation lines.
job = DFJob("Si", str, [scf_calculation], :ecutwfc => 20, #these flags will be set_ on all calculations that are passed to the job
            :verbosity => "high", :conv_thr => 1e-6; header = ["export OMP_NUM_THREADS=1"])

# Now the job can be submitted to be ran.
submit(job)
# this first saves the job and it's calculation files to the `job.local_dir` then pushes the `job.tt` file and the calculations to the `job.server_dir` on `job.server`, and tries to first submit it through slurm, otherwise just runs `bash job.tt`
# You can check the job.local_dir to see the calculation files and `job.tt` script.

# If the job was submitted via slurm, one can check if it's running through
slurm_isrunning(job)

# Hopefully everything went according to plan and we can retrieve the outputs
out = outputdata(job)
# Or, to retrieve the output of a specific calculation,
out = outputdata(job["scf"])

# Indeed specific calculations of the job can be accessed through their name  
# If for some reason we discover that smearing was needed, the flags related to that
# can be changed either job-wide (the flag will be set for every calculation that allows for that flag):
job[:occupations] = "smearing"
job[:degauss] = 0.001
# Or on a specific calculation
job["scf"][:occupations] = "smearing"
job["scf"][:degauss] = "smearing"

# Next step would be to run a bands calculation:
bands_in = gencalc_bands(job["scf"], high_symmetry_kpath(job.structure))
# Here we generate a bands calculation using the scf as the template, and the k-points
# generated from the structure symmetries.
# k-points can also be passed as e.g. :
bands_in = gencalc_bands(job["scf"], [(0.5, 0.5, 0.5, 10.0), (0.0, 0.0, 0.0, 1.0)])
# mimicking the usual structure in the bands pw.x calculation.

# Then we add the bands calculation to the job,
push!(job, bands_in)
# and since we already ran the scf we can do
set_flow!(job, "scf" => false, "bands" => true)

# The String => Bool pairs supplied to setflow will sequentially fuzzy match any name
# in which the string occurs. i.e.
set_flow!(job, "" => false, "bands" => true)
# would first set all calculations to not run, followed by setting the bands calculation to run
set_flow!(job, "" => false, "scf" => true)
# would run calculations with both scf and nscf names.

# Anyway now we submit the bands calculation
set_flow!(job, "" => false, "bands" => true)
submit(job)

#now the bandstructure can be plotted
bands = outputdata(job["bands"])[:bands]
#or
bands = readbands(job["bands"])
#or
bands = readbands(job)
# The last will search for a bands calculation in the job, otherwise it will use the output of an nscf calculation
# to generate the bands

# similarly 
fermi = outputdata(job["scf"])[:fermi]
fermi = readfermi(job["scf"])
fermi = readfermi(job)

using Plots

# Now the bands are plotted
plot(bands; fermi = fermi)

# To access the dos, we can first generate an nscf to generate a uniform k-grid calculation
# followed by a projwfc

push!(job, gencalc_nscf(job["scf"], (6, 6, 6)))
push!(job, gencalc_projwfc(job["nscf"], fermi - 10, fermi + 10, 0.1))

setflow!(job, "" => false, "prowjfc" => true, "nscf" => true)
submit(job)

# Now we can do the following:
plot(job, -5, 5)
# which will plot the band structure colored by the dos between fermi-5 and fermi+5

# This concludes this example.
