# # Advanced Usage
# In this tutorial we will continue from the job created in the [Basic Tutorial](@ref),
# and demonstrate some more advanced functionality that DFControl offers.

# To load a previously saved job we here provide a valid job directory with a job.sh script
# in it.
using DFControl

tjob = load(Job(joinpath(splitdir(pathof(DFControl))[1], "..", "docs","src","assets", "job")))#hide
tjob2 = load(Job(joinpath(splitdir(pathof(DFControl))[1], "..", "docs","src","assets", "Job2")))#hide
if false#hide
job = load(Servers.local_server(), Job("job"))
else#hide
    global job = deepcopy(tjob)#hide
    job.dir= "job" #hide
    job#hide
end#hide

# Since the job created in the [Basic Tutorial](@ref) was saved in the "job" directory this will work, see the section on [Jobs](@ref) for
# further details and options on how to load previously saved jobs.

# The next thing we may want to do is to change the directory where the job is running.
if false#hide
job.dir = "Job2"
else#hide
    global job = deepcopy(tjob2)#hide
    pop!(job)#hide
    pop!(job)#hide
    job#hide
end#hide
# With the `copy=true` flag we let DFControl know that not only to create and set the
# new directory, but also to copy the previous results and temporary files to the
# new directory so we don't have to rerun the scf calculation.

# Next we would like to plot the projected density of states.
# For that we create both an nscf calculation to get a uniform k-grid, and projwfc calculation.
push!(job, Calculations.gencalc_nscf(job["scf"], (6, 6, 6)))

# The second argument of gencalc_nscf is the kgrid. When passing a 3-Tuple,
# the code will assume that an explicit k-grid is requested, which can be verified by
data(job["nscf"], :k_points)

# Next we generate a projwfc calculation using the fermi level as a guide for the
# energy window.
# The arguments are structured as (template, Emin, Emax, deltaE) respectively.
fermi = readfermi(job)
push!(job, Calculations.gencalc_projwfc(job["nscf"], fermi - 10, fermi + 1, 0.1))
# Next we disable the bands calculation and run the new ones.
job["bands"].run = false
if false#hide
    submit(job)
else#hide
    global job = deepcopy(tjob2)#hide
end#hide

# ## [Plot Results](@id results_plotting)
using Plots
plot(job, -10, 1)
# As we can see, again DFControl identifies the additional information that is now present in the job, and uses it
# to display in the plot.

# In the demonstrated case we see that everything went according to plan, however, often things need to be changed
# in a trial and error way until the desired results are found.

# On common occurence is that calculation flags have to be set, or changed. This can be done in two ways
job[:ecutwfc] = 40.0
# will go through all the calculations of the job and set the flag if it is allowed, i.e. the flag will not
# be set in the projwfc calculation since it makes no sense.
job["bands"][:nbnd] = 30
# This will set a flag for one specific calculation, again checking whether the flag is valid, and the type
# will be converted to the correct one.

# In order to quickly specify what calculations to schedule and which not, one can use
Jobs.set_flow!(job, "" => false, "scf" => true)
# As we can see, only the scf and nscf calculations are scheduled to run now,
# this is because for each of the pairs in the arguments of `set_flow!`, every calculation inside
# the job for which the string occurs in the name will be set to run or not depending on the Bool.
