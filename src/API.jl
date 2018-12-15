include("inputAPI.jl")
#Basic Interaction with DFInputs
export flag, setflags!, rmflags!, data, setdata!, setdataoption!, exec, execs,
       setexecflags!, rmexecflags!, setexecdir!, runcommand, outputdata

#Extended Interaction with DFInputs
export setkpoints!, readbands

#generating new DFInputs
export gencalc_scf, gencalc_nscf, gencalc_bands

include("jobAPI.jl")
#Basic Job Control Functionality
export save, submit, abort, setflow!, setheaderword!, isrunning, progressreport,
       setserverdir!, setlocaldir!

#Basic Interaction with DFInputs inside DFJob
export searchinput, searchinputs, setcutoffs!, setname!

#Extended Interaction with DFInputs inside DFJob
export  addwancalc!, setwanenergies!

#Interacting with the Structure inside DFJob
export atom, atoms, setatoms!, setpseudos!, projections, setprojections!
