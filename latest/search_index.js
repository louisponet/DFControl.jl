var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#DFControl-1",
    "page": "Home",
    "title": "DFControl",
    "category": "section",
    "text": ""
},

{
    "location": "job.html#",
    "page": "Job",
    "title": "Job",
    "category": "page",
    "text": ""
},

{
    "location": "job.html#DFControl.DFJob",
    "page": "Job",
    "title": "DFControl.DFJob",
    "category": "type",
    "text": "Represents a full DFT job with multiple input files and calculations.\n\n\n\n"
},

{
    "location": "job.html#DFControl.DFJob",
    "page": "Job",
    "title": "DFControl.DFJob",
    "category": "type",
    "text": "DFJob(job_dir::String, T=Float64; job_fuzzy = \"job\", new_job_name=nothing, new_local_dir=nothing, server=getdefault_server(),server_dir=\"\")\n\nLoads and returns a local DFJob. If local_dir is not specified the job directory will be registered as the local one.\n\n\n\n"
},

{
    "location": "job.html#DFControl.DFJob",
    "page": "Job",
    "title": "DFControl.DFJob",
    "category": "type",
    "text": "DFJob(server_dir::String, local_dir::String, server=getdefault_server(); job_fuzzy=\"*job*\", new_job_name=\"\")\n\nPulls a server job to local directory and then loads it. A fuzzy search for the job file will be performed and the found input files will be pulled.\n\n\n\n"
},

{
    "location": "job.html#DFControl.abort-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.abort",
    "category": "method",
    "text": "abort(job::DFJob)\n\nWill look for the job id inside it\'s metadata and try to remove it from the server queue.\n\n\n\n"
},

{
    "location": "job.html#DFControl.addcalc!",
    "page": "Job",
    "title": "DFControl.addcalc!",
    "category": "function",
    "text": "Creates a new DFInput from the template with the new flags and new data, then adds it to the inputs of the job at the specified index.\n\n\n\n"
},

{
    "location": "job.html#DFControl.addcalc!-Tuple{DFControl.DFJob,Array{#s48,1} where #s48<:(NTuple{4,T} where T),Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.addcalc!",
    "category": "method",
    "text": "addcalc!(job::DFJob, kpoints::Vector{NTuple{4}}, newflags...; name=\"bands\", run=true, template=\"scf\")\n\nSearches for the given template and creates a bands calculation from it.\n\n\n\n"
},

{
    "location": "job.html#DFControl.addcalc!-Tuple{DFControl.DFJob,NTuple{6,T} where T,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.addcalc!",
    "category": "method",
    "text": "addcalc!(job::DFJob, kpoints::NTuple{6}, newflags...; name=\"scf\", run=true, template=\"nscf\")\n\nSearches for the given template and creates a bands calculation from it.\n\n\n\n"
},

{
    "location": "job.html#DFControl.addcalc!-Tuple{DFControl.DFJob,Tuple{T,T,T} where T,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.addcalc!",
    "category": "method",
    "text": "addcalc!(job::DFJob, kpoints::NTuple{3}, newflags...; name=\"nscf\", run=true, template=\"scf\")\n\nSearches for the given template and creates a bands calculation from it.\n\n\n\n"
},

{
    "location": "job.html#DFControl.addwancalc!-Tuple{DFControl.DFJob,DFControl.DFInput{DFControl.QE},Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.addwancalc!",
    "category": "method",
    "text": "addwancalc!(job::DFJob, nscf::DFInput{QE}, projections;\n                 Emin=-5.0,\n                 Epad=5.0,\n                 wanflags=SymAnyDict(),\n                 pw2wanexec=Exec(\"pw2wannier90.x\", nscf.execs[2].dir, nscf.execs[2].flags),\n                 wanexec=Exec(\"wannier90.x\", nscf.execs[2].dir),\n                 bands=readbands(nscf))\n\nAdds a wannier calculation to a job. For now only works with QE.\n\n\n\n"
},

{
    "location": "job.html#DFControl.atom",
    "page": "Job",
    "title": "DFControl.atom",
    "category": "function",
    "text": "Returns the ith atom with id atsym.\n\n\n\n"
},

{
    "location": "job.html#DFControl.atoms-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.atoms",
    "category": "method",
    "text": "atoms(job::DFJob)\n\nReturns a list the atoms in the structure.\n\n\n\n"
},

{
    "location": "job.html#DFControl.data-Tuple{DFControl.DFJob,Any,Symbol}",
    "page": "Job",
    "title": "DFControl.data",
    "category": "method",
    "text": "data(job::DFJob, calc_filenames, name::Symbol)\n\nLooks through the calculation filenames and returns the data with the specified symbol.\n\n\n\n"
},

{
    "location": "job.html#DFControl.execs-Tuple{DFControl.DFJob,Any}",
    "page": "Job",
    "title": "DFControl.execs",
    "category": "method",
    "text": "Returns the executables attached to a given input.\n\n\n\n"
},

{
    "location": "job.html#DFControl.flag-Tuple{DFControl.DFJob,Any,Symbol}",
    "page": "Job",
    "title": "DFControl.flag",
    "category": "method",
    "text": "flag(job::DFJob, names, flag_name::Symbol)\n\nLooks through the input names and returns the value of the specified flag.\n\n\n\n"
},

{
    "location": "job.html#DFControl.flag-Tuple{DFControl.DFJob,Symbol}",
    "page": "Job",
    "title": "DFControl.flag",
    "category": "method",
    "text": "flag(job::DFJob, flag_name::Symbol)\n\nLooks through all the calculations and returns the value of the specified flag.\n\n\n\n"
},

{
    "location": "job.html#DFControl.inputs-Tuple{DFControl.DFJob,Array{T,1} where T}",
    "page": "Job",
    "title": "DFControl.inputs",
    "category": "method",
    "text": "inputs(job::DFJob, names::Vector)\n\nReturns an array of the inputs that match the names.\n\n\n\n"
},

{
    "location": "job.html#DFControl.outpath-Tuple{DFControl.DFJob,String}",
    "page": "Job",
    "title": "DFControl.outpath",
    "category": "method",
    "text": "Finds the input corresponding to the name and returns the full output path.\n\n\n\n"
},

{
    "location": "job.html#DFControl.outputdata-Tuple{DFControl.DFJob,Array{DFControl.DFInput,1}}",
    "page": "Job",
    "title": "DFControl.outputdata",
    "category": "method",
    "text": "Finds the output files for each of the inputs of a job, and groups all found data into a dictionary.\n\n\n\n"
},

{
    "location": "job.html#DFControl.projections",
    "page": "Job",
    "title": "DFControl.projections",
    "category": "function",
    "text": "Returns the projections inside the job for the specified ith atom in the job with id atsym.\n\n\n\n"
},

{
    "location": "job.html#DFControl.projections-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.projections",
    "category": "method",
    "text": "Returns all the projections inside the job.\n\n\n\n"
},

{
    "location": "job.html#DFControl.rmflags!-Tuple{DFControl.DFJob,Array{#s48,1} where #s48<:AbstractString,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.rmflags!",
    "category": "method",
    "text": "rmflags!(job::DFJob, names, flags...)\n\nLooks through the input names and removes the specified flags.\n\n\n\n"
},

{
    "location": "job.html#DFControl.save",
    "page": "Job",
    "title": "DFControl.save",
    "category": "function",
    "text": "save(job::DFJob)\n\nSaves a DFJob, it\'s job file and all it\'s input files.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setatoms!-Tuple{DFControl.DFJob,Dict{Symbol,#s47} where #s47<:(Array{#s46,1} where #s46<:(GeometryTypes.Point{3,T} where T))}",
    "page": "Job",
    "title": "DFControl.setatoms!",
    "category": "method",
    "text": "setatoms!(job::DFJob, atoms::Dict{Symbol,<:Array{<:Point3,1}}, pseudo_setname=nothing, pseudospecifier=nothing, option=:angstrom)\n\nSets the data data with atomic positions to the new one. This is done for all calculations in the job that have that data. If default pseudopotentials are defined, a set can be specified, together with a fuzzy that distinguishes between the possible multiple pseudo strings in the pseudo set. These pseudospotentials are then set in all the calculations that need it. All flags which specify the number of atoms inside the calculation also gets set to the correct value.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setcell!-Tuple{DFControl.DFJob,StaticArrays.SArray{Tuple{3,3},T,2,9} where T}",
    "page": "Job",
    "title": "DFControl.setcell!",
    "category": "method",
    "text": "setcell_parameters!(job::DFJob, cell_::Mat3)\n\nsets the cell parameters of the structure in the job.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setdata!",
    "page": "Job",
    "title": "DFControl.setdata!",
    "category": "function",
    "text": "setdata!(job::DFJob, filenames, name, data, option=:none)\n\nAdds a block to the specified filenames.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setdata!-Tuple{DFControl.DFJob,Any,DFControl.InputData}",
    "page": "Job",
    "title": "DFControl.setdata!",
    "category": "method",
    "text": "setdata!(job::DFJob, names, data::InputData)\n\nAdds a block to the specified names.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setdata!-Tuple{DFControl.DFJob,Any,Symbol,Any}",
    "page": "Job",
    "title": "DFControl.setdata!",
    "category": "method",
    "text": "setdata!(job::DFJob, calc_filenames, data_block_name::Symbol, new_block_data; option=nothing)\n\nLooks through the calculation filenames and sets the data of the datablock with data_block_name to new_block_data. if option is specified it will set the block option to it.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setdataoption!-Tuple{DFControl.DFJob,Array{String,1},Symbol,Symbol}",
    "page": "Job",
    "title": "DFControl.setdataoption!",
    "category": "method",
    "text": "setdataoption!(job::DFJob, names::Vector{String}, dataname::Symbol, option::Symbol)\n\nsets the option of specified data in the specified inputs.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setdataoption!-Tuple{DFControl.DFJob,Symbol,Symbol}",
    "page": "Job",
    "title": "DFControl.setdataoption!",
    "category": "method",
    "text": "setdataoption!(job::DFJob, name::Symbol, option::Symbol)\n\nsets the option of specified data block in all calculations that have the block.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setexecflags!-Tuple{DFControl.DFJob,Any,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.setexecflags!",
    "category": "method",
    "text": "setexecflags!(job::DFJob, exec, flags...)\n\nGoes through the calculations of the job and if the name contains any of the inputnames it sets the exec flags to the specified ones.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setflags!-Tuple{DFControl.DFJob,Array{String,1},Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.setflags!",
    "category": "method",
    "text": "setflags!(job::DFJob, names::Vector{String}, flags...; print=true)\n\nSets the flags in the names to the flags specified. This only happens if the specified flags are valid for the names. If necessary the correct control block will be added to the calculation (e.g. for QEInputs).\n\nThe values that are supplied will be checked whether they are valid.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setflow!-Tuple{DFControl.DFJob,Array{String,1},Any}",
    "page": "Job",
    "title": "DFControl.setflow!",
    "category": "method",
    "text": "setflow!(job::DFJob, filenames::Array{String,1}, should_run)\n\nGoes throug the calculation filenames and sets whether it should run or not.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setflow!-Tuple{DFControl.DFJob,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.setflow!",
    "category": "method",
    "text": "setflow!(job::DFJob, should_runs...)\n\nSets whether or not calculations should be run. Calculations are specified using their indices.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setheaderword!-Tuple{DFControl.DFJob,String,String}",
    "page": "Job",
    "title": "DFControl.setheaderword!",
    "category": "method",
    "text": "setheaderword!(job::DFJob, word::String, new_word::String)\n\nReplaces the specified word in the header with the new word.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setkpoints!-Tuple{DFControl.DFJob,Any,Any}",
    "page": "Job",
    "title": "DFControl.setkpoints!",
    "category": "method",
    "text": "setkpoints!(job::DFJob, n, k_points)\n\nsets the data in the k point DataBlock inside the specified inputs.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setlocaldir!-Tuple{Any,Any}",
    "page": "Job",
    "title": "DFControl.setlocaldir!",
    "category": "method",
    "text": "Sets the local dir of the job.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setprojections!-Tuple{DFControl.DFJob,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.setprojections!",
    "category": "method",
    "text": "sets the projections of the specified atoms inside the job structure.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setpseudos!",
    "page": "Job",
    "title": "DFControl.setpseudos!",
    "category": "function",
    "text": "sets the pseudopotentials to the specified one in the default pseudoset.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setserverdir!-Tuple{Any,Any}",
    "page": "Job",
    "title": "DFControl.setserverdir!",
    "category": "method",
    "text": "Sets the server dir of the job.\n\n\n\n"
},

{
    "location": "job.html#DFControl.setwanenergies!-Tuple{DFControl.DFJob,AbstractFloat,Any}",
    "page": "Job",
    "title": "DFControl.setwanenergies!",
    "category": "method",
    "text": "Automatically calculates and sets the wannier energies. This uses the projections, Emin and the bands to infer the other limits. Epad allows one to specify the padding around the inner and outer energy windows\n\n\n\n"
},

{
    "location": "job.html#DFControl.submit-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.submit",
    "category": "method",
    "text": "submit(job::DFJob; server=job.server, server_dir=job.server_dir)\n\nSaves the job locally, and then either runs it locally using qsub (when job.server == \"localhost\") or sends it to the specified job.server in job.server_dir, and submits it using qsub on the server.\n\n\n\n"
},

{
    "location": "job.html#DFControl.undo!-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.undo!",
    "category": "method",
    "text": "undo!(job::DFJob)\n\nUndos the last set to the calculations of the job.\n\n\n\n"
},

{
    "location": "job.html#DFControl.undo-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.undo",
    "category": "method",
    "text": "undo(job::DFJob)\n\nReturns the previous state of the job.\n\n\n\n"
},

{
    "location": "job.html#Job-1",
    "page": "Job",
    "title": "Job",
    "category": "section",
    "text": "Modules = [DFControl]\nPages   = [\"job.jl\"]\nPrivate = false\nOrder   = [:type, :function]"
},

{
    "location": "input.html#",
    "page": "Inputs",
    "title": "Inputs",
    "category": "page",
    "text": ""
},

{
    "location": "input.html#Inputs-1",
    "page": "Inputs",
    "title": "Inputs",
    "category": "section",
    "text": ""
},

{
    "location": "input.html#DFControl.DFInput-Tuple{DFControl.DFInput,Any,Vararg{Any,N} where N}",
    "page": "Inputs",
    "title": "DFControl.DFInput",
    "category": "method",
    "text": "DFInput(template::DFInput, name, newflags...; runcommand=template.runcommand, run=true)\n\nCreates a new DFInput from a template DFInput, setting the newflags in the new one.\n\n\n\n"
},

{
    "location": "input.html#DFControl.outputdata-Tuple{DFControl.DFInput}",
    "page": "Inputs",
    "title": "DFControl.outputdata",
    "category": "method",
    "text": "Returns the outputdata for the input.\n\n\n\n"
},

{
    "location": "input.html#DFControl.rmflags!-Tuple{DFControl.DFInput,Vararg{Any,N} where N}",
    "page": "Inputs",
    "title": "DFControl.rmflags!",
    "category": "method",
    "text": "rmflags!(input::DFInput, flags...)\n\nRemove the specified flags.\n\n\n\n"
},

{
    "location": "input.html#DFControl.setdata!-Tuple{DFControl.DFInput,DFControl.InputData}",
    "page": "Inputs",
    "title": "DFControl.setdata!",
    "category": "method",
    "text": "setdata!(input::DFInput, data::InputData)\n\nAdds the given data to the input. Should put it in the correct arrays.\n\n\n\n"
},

{
    "location": "input.html#DFControl.setdata!-Tuple{DFControl.DFInput,Symbol,Any}",
    "page": "Inputs",
    "title": "DFControl.setdata!",
    "category": "method",
    "text": "setdata!(input::DFInput, block_name::Symbol, new_block_data; option=nothing, print=true)\n\nsets the data of the specified \'InputData\' to the new data. Optionally also sets the \'InputData\' option.\n\n\n\n"
},

{
    "location": "input.html#DFControl.setdataoption!-Tuple{DFControl.DFInput,Symbol,Symbol}",
    "page": "Inputs",
    "title": "DFControl.setdataoption!",
    "category": "method",
    "text": "setdataoption!(input::DFInput, name::Symbol, option::Symbol; print=true)\n\nSets the option of specified data.\n\n\n\n"
},

{
    "location": "input.html#DFControl.setflags!-Union{Tuple{DFControl.DFInput{T},Vararg{Any,N} where N}, Tuple{T}} where T",
    "page": "Inputs",
    "title": "DFControl.setflags!",
    "category": "method",
    "text": "setflags!(input::DFInput, flags...; print=true)\n\nSets the specified flags in the input.\n\n\n\n"
},

{
    "location": "input.html#DFControl.setkpoints!-Tuple{DFControl.DFInput{DFControl.Wannier90},Tuple{Int64,Int64,Int64}}",
    "page": "Inputs",
    "title": "DFControl.setkpoints!",
    "category": "method",
    "text": "setkpoints!(input::DFInput, k_grid)\n\nSets the kpoints of the input. Will automatically generate the kgrid values if necessary.\n\n\n\n"
},

{
    "location": "input.html#DFInput-1",
    "page": "Inputs",
    "title": "DFInput",
    "category": "section",
    "text": "Modules = [DFControl]\nPages   = [\"input.jl\"]\nPrivate = false\nOrder   = [:type, :function]"
},

{
    "location": "structure.html#",
    "page": "Structure",
    "title": "Structure",
    "category": "page",
    "text": ""
},

{
    "location": "structure.html#DFControl.atoms-Tuple{DFControl.AbstractStructure,Symbol}",
    "page": "Structure",
    "title": "DFControl.atoms",
    "category": "method",
    "text": "Returns all the atoms inside the structure with the specified symbol\n\n\n\n"
},

{
    "location": "structure.html#DFControl.setprojections!-Tuple{DFControl.Structure,Vararg{Any,N} where N}",
    "page": "Structure",
    "title": "DFControl.setprojections!",
    "category": "method",
    "text": "sets the projections of the specified atoms.\n\n\n\n"
},

{
    "location": "structure.html#Structure-1",
    "page": "Structure",
    "title": "Structure",
    "category": "section",
    "text": "Modules = [DFControl]\nPages   = [\"structure.jl\"]\nPrivate = false\nOrder   = [:type, :function]"
},

{
    "location": "atom.html#",
    "page": "Atom",
    "title": "Atom",
    "category": "page",
    "text": ""
},

{
    "location": "atom.html#Atom-1",
    "page": "Atom",
    "title": "Atom",
    "category": "section",
    "text": "Modules = [DFControl]\nPages   = [\"atom.jl\"]\nPrivate = false\nOrder   = [:type, :function]"
},

{
    "location": "fileio.html#",
    "page": "FileIO",
    "title": "FileIO",
    "category": "page",
    "text": ""
},

{
    "location": "fileio.html#DFControl.read_fermi_from_qe_output",
    "page": "FileIO",
    "title": "DFControl.read_fermi_from_qe_output",
    "category": "function",
    "text": "read_fermi_from_qe_output(filename::String,T=Float64)\n\nReads the Fermi level from a Quantum Espresso scf calculation output file (if there is one).\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_ks_from_qe_output",
    "page": "FileIO",
    "title": "DFControl.read_ks_from_qe_output",
    "category": "function",
    "text": "read_ks_from_qe_output(filename::String, T=Float64)\n\nRead k-points from a Quantum Espresso bands output file in cartesian (2pi/alat in Angstrom^-1!) and crystalline coordinates. Returns (k_points_cart,k_points_cryst).\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_qe_bands_file",
    "page": "FileIO",
    "title": "DFControl.read_qe_bands_file",
    "category": "function",
    "text": "read_qe_bands(filename::String, T=Float64)\n\nReads the output file of a \'bands\' calculation in Quantum Espresso. Returns an array of DFBands each with the same k_points and their respective energies.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_qe_input",
    "page": "FileIO",
    "title": "DFControl.read_qe_input",
    "category": "function",
    "text": "read_qe_input(filename, T=Float64; exec=\"pw.x\",  runcommand=\"\", run=true, structure_name=\"NoName\")\n\nReads a Quantum Espresso input file. The exec get\'s used to find which flags are allowed in this input file, and convert the read values to the correct Types. Returns a DFInput{QE} and the Structure that is found in the input.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_qe_kpdos",
    "page": "FileIO",
    "title": "DFControl.read_qe_kpdos",
    "category": "function",
    "text": "read_qe_kpdos(filename::String,column=1;fermi=0)\n\nReads the k_resolved partial density of states from a Quantum Espresso projwfc output file. Only use this if the flag kresolveddos=true in the projwfc input file!! The returned matrix can be readily plotted using heatmap() from Plots.jl! Optional input: column = 1 (column of the output, 1 = first column after ik and E) fermi  = 0 (possible fermi offset of the read energy values) Return:         Array{Float64,2}(length(k_points),length(energies)) , (ytickvals,yticks)\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_qe_output",
    "page": "FileIO",
    "title": "DFControl.read_qe_output",
    "category": "function",
    "text": "read_qe_output(filename::String, T=Float64)\n\nReads a generic quantum espresso input, returns a dictionary with all found data in the file. Possible keys:\n\n:fermi\n:polarization\n:pol_mod\n:k_cryst\n:k_cart\n:alat\n:cell_parameters\n:pos_option\n:atomic_positions\n:total_force\n:colin_mag_moments\n:bands\n:accuracy\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_qe_pdos",
    "page": "FileIO",
    "title": "DFControl.read_qe_pdos",
    "category": "function",
    "text": "read_qe_pdos(filename::String, column=1; fermi=0)\n\nReads partial dos file. One can specify the column of values to read.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_qe_polarization",
    "page": "FileIO",
    "title": "DFControl.read_qe_polarization",
    "category": "function",
    "text": "read_qe_polarization(filename::String, T=Float64)\n\nReturns the polarization and modulus.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_wannier_input",
    "page": "FileIO",
    "title": "DFControl.read_wannier_input",
    "category": "function",
    "text": "read_wannier_input(filename::String, T=Float64; runcommand= Exec(\"\"), run=true, exec=Exec(\"wannier90.x\"), structure_name=\"NoName\")\n\nReads a DFInput{Wannier90} and the included Structure from a WANNIER90 input file.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.save",
    "page": "FileIO",
    "title": "DFControl.save",
    "category": "function",
    "text": "save(input::DFInput{Wannier90}, structure, filename::String=inpath(input))\n\nWrites the DFInput{Wannier90} and structure to a file, that can be interpreted by WANNIER90. The atoms in the structure must have projections defined.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.save",
    "page": "FileIO",
    "title": "DFControl.save",
    "category": "function",
    "text": "save(input::DFInput{QE}, structure, filename::String=inpath(input))\n\nWrites a Quantum Espresso input file.\n\n\n\n"
},

{
    "location": "fileio.html#FileIO-1",
    "page": "FileIO",
    "title": "FileIO",
    "category": "section",
    "text": "Modules = [DFControl]\nPages   = [\"fileio.jl\"]\nPrivate = false\nOrder   = [:type, :function]"
},

{
    "location": "defaults.html#",
    "page": "Defaults",
    "title": "Defaults",
    "category": "page",
    "text": ""
},

{
    "location": "defaults.html#DFControl.configuredefault_pseudos-Tuple{}",
    "page": "Defaults",
    "title": "DFControl.configuredefault_pseudos",
    "category": "method",
    "text": "configuredefault_pseudos(server = getdefault_server(), pseudo_dirs=getdefault_pseudodirs())\n\nReads the specified default_pseudo_dirs on the default_server and sets up the default_pseudos variable, and also adds all the entries to the user_defaults.jl file.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.getdefault_pseudo",
    "page": "Defaults",
    "title": "DFControl.getdefault_pseudo",
    "category": "function",
    "text": "getdefault_pseudo(atom::Symbol, pseudo_setname=:default; pseudospecifier=nothing)\n\nReturns the pseudo potential string linked to the atom.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.getdefault_pseudodirs-Tuple{}",
    "page": "Defaults",
    "title": "DFControl.getdefault_pseudodirs",
    "category": "method",
    "text": "getdefault_pseudodirs()\n\nReturns the default pseudo dirs if it\'s defined. If it is not defined return nothing.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.getdefault_server-Tuple{}",
    "page": "Defaults",
    "title": "DFControl.getdefault_server",
    "category": "method",
    "text": "getdefault_server()\n\nReturns the default server if it\'s defined. If it is not defined return \"\".\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.removedefault_input-Tuple{Symbol}",
    "page": "Defaults",
    "title": "DFControl.removedefault_input",
    "category": "method",
    "text": "removedefault_input(input::Symbol)\n\nRemove input from the default_inputs variable. Also removes the stored input file.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.removedefault_pseudodir-Tuple{Symbol}",
    "page": "Defaults",
    "title": "DFControl.removedefault_pseudodir",
    "category": "method",
    "text": "removedefault_pseudodir(pseudo_symbol::Symbol)\n\nRemoves entry with flag pseudo_symbol from the default_pseudodirs and user_defaults.jl file.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.removedefault_pseudos-Tuple{Symbol}",
    "page": "Defaults",
    "title": "DFControl.removedefault_pseudos",
    "category": "method",
    "text": "removedefault_pseudos(pseudo_symbol::Symbol)\n\nRemoves all pseudo entries with flag pseudo_symbol from the default_pseudos.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.setdefault_input-Union{Tuple{DFControl.DFInput{T},Any,Symbol}, Tuple{T}} where T",
    "page": "Defaults",
    "title": "DFControl.setdefault_input",
    "category": "method",
    "text": "setdefault_input(input::dfinput, structure, calculation::Symbol)\n\nAdds the input to the default_inputs variable, and writes it to a file in user_defaults folder to be read every time on load.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.setdefault_jobheader-Tuple{Any}",
    "page": "Defaults",
    "title": "DFControl.setdefault_jobheader",
    "category": "method",
    "text": "setdefault_jobheader(lines)\n\nSets the header that will get added to each job.tt file, if no other header was specified.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.setdefault_pseudodir-Tuple{Symbol,String}",
    "page": "Defaults",
    "title": "DFControl.setdefault_pseudodir",
    "category": "method",
    "text": "setdefault_pseudodir(pseudo_symbol::Symbol, dir::String)\n\nAdds an entry inside the default_pseudodirs with flag pseudo_symbol, and adds it to the user_defaults.jl file.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.setdefault_server-Tuple{String}",
    "page": "Defaults",
    "title": "DFControl.setdefault_server",
    "category": "method",
    "text": "setdefault_server(server::String)\n\nSets the default server variable, and also adds it to the user_defaults.jl file.\n\n\n\n"
},

{
    "location": "defaults.html#Defaults-1",
    "page": "Defaults",
    "title": "Defaults",
    "category": "section",
    "text": "Modules = [DFControl]\nPages   = [\"defaults.jl\"]\nPrivate = false\nOrder   = [:type, :function]"
},

{
    "location": "utils.html#",
    "page": "Utils",
    "title": "Utils",
    "category": "page",
    "text": ""
},

{
    "location": "utils.html#DFControl.kgrid-Tuple{Any,Any,Any,Type{DFControl.QE}}",
    "page": "Utils",
    "title": "DFControl.kgrid",
    "category": "method",
    "text": "kgrid(na, nb, nc, input)\n\nReturns an array of k-grid points that are equally spaced, input can be either :wan or :nscf, the returned grids are appropriate as inputs for wannier90 or an nscf calculation respectively.\n\n\n\n"
},

{
    "location": "utils.html#Utils-1",
    "page": "Utils",
    "title": "Utils",
    "category": "section",
    "text": "Modules = [DFControl]\nPages   = [\"utils.jl\"]\nPrivate = false\nOrder   = [:type, :function]"
},

]}
