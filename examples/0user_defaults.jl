using DFControl

#This example shows how to setup and interact with the user_defaults.jl file which will get loaded each time you do `using DFControl`
setdefault_server("blabla@123.312.124")
setdefault_server("ponet@10.255.9.115")
setdefault_pseudodir(:sssp, "~/pseudos/sssp")
setdefault_pseudodir(:pbesol, "~/pseudos/pbesol/PSEUDOPOTENTIALS")
setdefault_pseudodir(:pbesolrel, "~/pseudos/pbesolrel/PSEUDOPOTENTIALS")
configure_defaultpseudos()
#look at the user_defaults file inside ~/.julia/v0.x/DFControl/user_defaults/

#We also allow for other defaults to be defined upon loading the package, this can be anything in the form of `var = expr`.

@add_default phd_dir = "/Users/user/Documents/physics/PhD/"
