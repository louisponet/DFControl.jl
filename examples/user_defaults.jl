using DFControl

#This example shows how to setup and interact with the user_defaults.jl file which will get loaded each time you do `using DFControl`
set_default_server("blabla@123.312.124")
set_default_pseudo_dir(:sssp, "~/pseudos/sssp")
set_default_pseudo_dir(:pbesol, "~/pseudos/pbesol/PSEUDOPOTENTIALS")
set_default_pseudo_dir(:pbesolrel, "~/pseudos/pbesolrel/PSEUDOPOTENTIALS")
configure_default_pseudos()
#look at the user_defaults file inside ~/.julia/v0.x/DFControl/user_defaults/

#We also allow for other defaults to be defined upon loading the package, this can be anything in the form of `var = expr`.

@add_default phd_dir = "/Users/ponet/Documents/Fysica/PhD/"
