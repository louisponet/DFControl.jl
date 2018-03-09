#This example assumes that you ran through the first one.
# It will show how to load a job from a server, and change flags.

using DFControl

server_dir = "Si/"
local_dir  = "/home/ponet/Documents/Si"
job = load_server_job(server_dir, local_dir)
#shorthand:
job = ldsj(server_dir, local_dir) #for other shorthands see DFControl/src/shorthands.jl

#change flags
change_flags!(job, :ecutwfc => 25.0)  #more `Pair{Symbol, Any}`'s can be given as varargs.
#shorthand: chfls!
#change_flags! will only go through already set flags to change them,
#if you want to set new flags simply do (it can also be used to change already set flags):

set_flags!(job, :ecut_rho => 80.) #won't work, it goes through the QE documentation to find the allowed flags for the input files in the job, and also tries to convert the given value to what it should be.

set_flags!(job, :ecutrho => 80., :diagonalization => "'david'")
#shorthand: stfls!

#now we might not want to run all the calculations again.
#one can change the job "flow":
change_flow!(job, "scf" => false) #again, more `Pair{String, Bool}`'s can be given as varargs
#This looks through all the calculations inside the job, if the `calculation.name` (derived from it's filename) contains a given `String`, it will set `calculation.run` to the corresponding `Bool`.
#If a job has calculations/inputs that don't need to run, their line will be commented out in `job.tt` script.

#now we can resubmit the job with the changed flags
sbmj(job)
