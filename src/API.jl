export Structure, Atom, Pseudo, DFTU
#Units exported
export Ang, e₀, kₑ, a₀, Eₕ, Ry


include("inputAPI.jl")
export DFInput
#Basic Interaction with DFInputs
export flag, set_flags!, rm_flags!, data, set_data!, set_dataoption!, exec, execs,
       set_execflags!, rmexecflags!, set_execdir!, runcommand, outputdata, set_name!

#Extended Interaction with DFInputs
export set_kpoints!, readbands, readfermi, set_wanenergies!, isconverged

#generating new DFInputs
export gencalc_scf, gencalc_vcrelax, gencalc_nscf, gencalc_bands,
       gencalc_projwfc, gencalc_wan

include("jobAPI.jl")
#Basic Job Control Functionality
export save, submit, abort, set_flow!, set_headerword!, isrunning, progressreport,
       set_serverdir!, set_localdir!, structure, scale_cell!, volume,
       switch_version, version, versions, registered_jobs, rm_version!, rm_versions!, rm_tmp_dirs!,
       cleanup

#Basic Interaction with DFInputs inside DFJob
export searchinput, searchinputs, set_cutoffs!

#Interacting with the Structure inside DFJob
export atom, atoms, set_atoms!, set_pseudos!, projections, set_projections!, cell, a, b, c,
	   set_magnetization!, symmetry_operators, international, niggli_reduce, update_geometry!,
	   high_symmetry_kpath, high_symmetry_kpoints

export create_supercell

# Atom interface functions
export name, position_cart, position_cryst, element, pseudo, projections, magnetization, dftu
export distance, set_position!, scale_bondlength!, polyhedron

#Bands related functionality
export bandgap

#other postprocessing
export pdos

include("serverAPI.jl")
export qstat, watch_qstat, qdel

#Slurm interactions
export slurm_history_jobdir, slurm_jobid, slurm_isrunning, slurm_mostrecent, slurm_isqueued

include("documentationAPI.jl")
export documentation
