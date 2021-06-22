#This example assumes that you ran through the first one.
# It will show how to load a job from a server, and set flags.

using DFControl

server_dir = "Si/"
local_dir  = "/home/ponet/Documents/Si"
job        = DFJob(server_dir, local_dir)

#set flags
setflags!(job, :ecutwfc => 25.0)  #more `Pair{Symbol, Any}`'s can be given as varargs.
#setflags! will only go through already set flags to set them,
#if you want to set new flags simply do (it can also be used to set already set flags):

setflags!(job, :ecut_rho => 80.0) #won't work, it goes through the QE documentation to find the allowed flags for the calculation files in the job, and also tries to convert the given value to what it should be.

setflags!(job, :ecutrho => 80.0, :diagonalization => "david")

#now we might not want to run all the calculations again.
#one can set the job "flow":
setflow!(job, "scf" => false) #again, more `Pair{String, Bool}`'s can be given as varargs
#This looks through all the calculations inside the job, if the `calculation.name` (derived from it's filename) contains a given `String`, it will set `calculation.run` to the corresponding `Bool`.
#If a job has calculations/calculations that don't need to run, their line will be commented out in `job.tt` script.

#now we can resubmit the job with the setd flags
submit(job)
