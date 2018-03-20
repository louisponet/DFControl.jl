using DFControl

#This example goes through how one would create a job from scratch, using a .cif file for the Structure.

#First go to your favourite Crystal structure database and download the .cif you want to use.
#e.g. Si (F d -3 m :1) : http://www.crystallography.net/cod/9011998.cif

#We want to run an 'scf' and 'bands' calculation using QuantumEspresso.


#Variables that will be passed to the `DFJob` constructor.

name = "Si" #this name will also be given to the Structure inside the DFJob
local_dir = "/home/ponet/Documents/Si"
server_dir = "/home/ponet/Si"
bin_dir = "/usr/local/bin" #this is defaulted to the users bin dir = "~/bin/", it is the directory where pw.x etc will be called from
run_command = Exec("mpirun", bin_dir, Dict{Symbol, Any}(:np => 24)) #this is the run command before the executable of the calculation and it's flags

pseudo_set = :pbesol #nonrelativistic calculation ( assumes you set up the pseudos, as demonstrated in the README)
pseudo_specifier = "paw" #this selects the correct pseudo if multiple belong to the pseudo_set. If you don't specify this, the first one in the set will be used.

#The header holds all the other information inside a job scriptfile that is not recognized as input and output.
header = ["#SBATCH -N 1", "#SBATCH --ntasks-per-node=24",
          "#SBATCH --time=24:00:00", "#SBATCH -p defpart",
          "module load open-mpi/gcc/1.10.4-hfi", "module load mkl/2016.1.056"
         ]

#The various calculations we want to run and the flags and data to pass to them are defined in two ways:
#   - calculation specific flags and data are associated with the calculation they belong to
#   - common flags can be defined as Pair{Symbol, Any} varargs at the end of the constructor call.
scf_data = Dict(:k_points => (6, 6, 6, 1, 1, 1), :flags => [:verbosity => "'low'"])
bands_data = Dict(:k_points => [(0.5, 0.5, 0.5, 100.),
                                (0.0, 0.0, 0.0, 100.),
                                (0.0, 0.5, 0.0, 1.)],
                  :flags => [:verbosity => "'high'", :nbnd => 8])

calculations = [:scf => scf_data, :bands => bands_data] #the order here is the order in which the calculations will run!

#Now we load the cif file and create a `DFJob` from it.

job = DFJob(name, local_dir, "/home/ponet/Downloads/9011998.cif", calculations,
      :prefix       => "'silicon'",
      :restart_mode => "'from_scratch'",
      :ecutwfc      => 18.0,
      :mixing_mode  => "'plain'",
      :mixing_beta  => 0.7,
      :conv_thr     => 1.0e-8,
      #kwargs
      server_dir  = server_dir,
      bin_dir     = bin_dir,
      run_command = run_command,
      header      = header,
      pseudo_set  = :pbesol,
      pseudo_specifier = pseudo_specifier
     )
# An additional kwarg is `server=get_default_server()`, which is set to the server you have defined while following the setup in README.
# This can ofcourse be changed to a different server.

#Now the job can be submitted to the server to run.
submit_job(job)
#this first saves the job and it's input files to the `job.local_dir` then pushes the `job.tt` file and the inputs to the `job.server_dir` on `job.server`, and runs the `qsub job.tt` command.
#You can check the job.local_dir to see the input files and `job.tt` script.

#you can observe the slurm qstat by doing
qstat()
#or watch it by
watch_qstat()
#these default to run the commands on the default server

#hopefully everything went according to plan and we can watch our outputs
outputs = pull_outputs(job)

#now the bandstructure can be plotted
bands = read_qe_output(outputs[2])[:bands]
#alt:
bands = read_qe_bands_file(outputs[2])

fermi = read_qe_output(outputs[1])[:fermi]
#alt:
fermi = read_fermi_from_qe_output(outputs[1])
using Plots

plot(bands, fermi=fermi)