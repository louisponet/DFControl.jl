using DFControl

#This example shows how to setup and interact with the user_defaults.jl file which will get loaded each time you do `using DFControl`
setdefault_server("blabla@123.312.124")
getdefault_pseudodirs()
setdefault_pseudodir(:sssp, "~/pseudos/sssp")
setdefault_pseudodir(:pbesol, "~/pseudos/pbesol/PSEUDOPOTENTIALS")
setdefault_pseudodir(:pbesolrel, "~/pseudos/pbesolrel/PSEUDOPOTENTIALS")
configuredefault_pseudos()
#look at the user_defaults.jl file inside ~/.julia/config/DFControl

# stored default pseudos can be accessed through e.g.
getdefault_pseudo(:Si, :sssp)
