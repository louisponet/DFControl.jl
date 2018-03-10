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
    "location": "job.html#DFControl.add_bands_calculation!-Union{Tuple{DFControl.DFJob,Array{Array{T,1},1}}, Tuple{T}} where T<:AbstractFloat",
    "page": "Job",
    "title": "DFControl.add_bands_calculation!",
    "category": "method",
    "text": "add_bands_calculation!(job::DFJob, k_path::Array{Array{<:AbstractFloat,1},1})\n\nChecks if there is an scf calculation in the job and takes it\'s inputs to generate a bands calculation along the given k-path.\n\n\n\n"
},

{
    "location": "job.html#DFControl.add_block!-Tuple{DFControl.DFJob,Any,DFControl.Block}",
    "page": "Job",
    "title": "DFControl.add_block!",
    "category": "method",
    "text": "add_block!(job::DFJob, filenames, block::Block)\n\nAdds a block to the specified filenames.\n\n\n\n"
},

{
    "location": "job.html#DFControl.add_calculation!",
    "page": "Job",
    "title": "DFControl.add_calculation!",
    "category": "function",
    "text": "add_calculation!(job::DFJob, input::DFInput, index::Int=length(job.calculations)+1; run_command=input.run_command, filename=input.filename)\n\nAdds a calculation to the job, at the specified index.\n\n\n\n"
},

{
    "location": "job.html#DFControl.add_data!",
    "page": "Job",
    "title": "DFControl.add_data!",
    "category": "function",
    "text": "add_data!(job::DFJob, filenames, block_symbol, data, option=:none)\n\nAdds a block to the specified filenames.\n\n\n\n"
},

{
    "location": "job.html#DFControl.add_wan_calc!-Tuple{DFControl.DFJob,Any}",
    "page": "Job",
    "title": "DFControl.add_wan_calc!",
    "category": "method",
    "text": "add_wan_calc!(job::DFJob, k_grid;\n                   nscf_file          = \"nscf.in\",\n                   wan_file           = \"wan.win\",\n                   pw2wan_file        = \"pw2wan.in\",\n                   wan_run_command    = \"~/bin/wannier90.x \",\n                   pw2wan_run_command = \"mpirun -np 24 ~/bin/pw2wannier90.x\",\n                   wan_flags          = nothing,\n                   pw2wan_flags       = nothing)\n\nAdds a wannier calculation to a job. For now only works with QE.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_atoms!-Tuple{DFControl.DFJob,Dict{Symbol,#s35} where #s35<:(Array{#s34,1} where #s34<:(GeometryTypes.Point{3,T} where T))}",
    "page": "Job",
    "title": "DFControl.change_atoms!",
    "category": "method",
    "text": "change_atoms!(job::DFJob, atoms::Dict{Symbol,<:Array{<:Point3,1}}, pseudo_set_name=:default, pseudo_specifier=nothing, option=:angstrom)\n\nSets the data blocks with atomic positions to the new one. This is done for all calculations in the job that have that data. If default pseudopotentials are defined, a set can be specified, together with a fuzzy that distinguishes between the possible multiple pseudo strings in the pseudo set. These pseudospotentials are then set in all the calculations that need it. All flags which specify the number of atoms inside the calculation also gets set to the correct value.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_cell!-Tuple{DFControl.DFJob,Array{T,2} where T}",
    "page": "Job",
    "title": "DFControl.change_cell!",
    "category": "method",
    "text": "change_cell_parameters!(job::DFJob, cell_param::Matrix)\n\nChanges the cell parameters of the structure in the job.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_data!-Tuple{DFControl.DFJob,Any,Symbol,Any}",
    "page": "Job",
    "title": "DFControl.change_data!",
    "category": "method",
    "text": "change_data!(job::DFJob, calc_filenames, data_block_name::Symbol, new_block_data; option=nothing)\n\nLooks through the calculation filenames and changes the data of the datablock with data_block_name to new_block_data. if option is specified it will set the block option to it.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_data_option!-Tuple{DFControl.DFJob,Array{String,1},Symbol,Symbol}",
    "page": "Job",
    "title": "DFControl.change_data_option!",
    "category": "method",
    "text": "change_data_option!(job::DFJob, filenames::Array{String,1}, block_symbol::Symbol, option::Symbol)\n\nChanges the option of specified data block in the specified calculations.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_data_option!-Tuple{DFControl.DFJob,Symbol,Symbol}",
    "page": "Job",
    "title": "DFControl.change_data_option!",
    "category": "method",
    "text": "change_data_option!(job::DFJob, block_symbol::Symbol, option::Symbol)\n\nChanges the option of specified data block in all calculations that have the block.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_filename!-Tuple{DFControl.DFJob,String,String}",
    "page": "Job",
    "title": "DFControl.change_filename!",
    "category": "method",
    "text": "change_filename!(job::DFJob, old_filename::String, new_filename::String)\n\nChanges the filename from the old to the new one.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_flags!-Tuple{DFControl.DFJob,Array{String,1},Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.change_flags!",
    "category": "method",
    "text": "change_flags!(job::DFJob, calc_filenames::Vector{String}, new_flag_data...)\n\nLooks through the given calculations for the specified flags. If any that match and have the same types are found, they will get replaced by the new ones.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_flags!-Tuple{DFControl.DFJob,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.change_flags!",
    "category": "method",
    "text": "change_flags!(job::DFJob, new_flag_data...)\n\nLooks through all the calculations for the specified flags. If any that match and have the same types are found, they will get replaced by the new ones.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_flow!-Tuple{DFControl.DFJob,Array{String,1},Any}",
    "page": "Job",
    "title": "DFControl.change_flow!",
    "category": "method",
    "text": "change_flow!(job::DFJob, filenames::Array{String,1}, should_run)\n\nGoes throug the calculation filenames and sets whether it should run or not.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_flow!-Tuple{DFControl.DFJob,Union{Array{Tuple{String,Bool},1}, Dict{String,Bool}}}",
    "page": "Job",
    "title": "DFControl.change_flow!",
    "category": "method",
    "text": "change_flow!(job::DFJob, should_runs::Union{Dict{String,Bool},Array{Tuple{String,Bool}}})\n\nRuns through the calculation filenames and sets whether it should run or not.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_flow!-Tuple{DFControl.DFJob,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.change_flow!",
    "category": "method",
    "text": "change_flow!(job::DFJob, should_runs...)\n\nSets whether or not calculations should be run. Calculations are specified using their indices.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_header_word!-Tuple{DFControl.DFJob,String,String}",
    "page": "Job",
    "title": "DFControl.change_header_word!",
    "category": "method",
    "text": "replace_header_word!(job::DFJob, word::String, new_word::String)\n\nReplaces the specified word in the header with the new word.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_kpoints!-Tuple{DFControl.DFJob,Any,Any}",
    "page": "Job",
    "title": "DFControl.change_kpoints!",
    "category": "method",
    "text": "change_kpoints!(job::DFJob, calc_filename, k_points)\n\nChanges the data in the k point DataBlock inside the specified calculation.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_projections!-Tuple{DFControl.DFJob,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.change_projections!",
    "category": "method",
    "text": "Changes the projections of the specified atoms inside the job structure.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_pseudo_set!",
    "page": "Job",
    "title": "DFControl.change_pseudo_set!",
    "category": "function",
    "text": "Changes the pseudopotentials to the specified one in the default pseudo_set.\n\n\n\n"
},

{
    "location": "job.html#DFControl.change_run_command!-Tuple{DFControl.DFJob,Any,Any}",
    "page": "Job",
    "title": "DFControl.change_run_command!",
    "category": "method",
    "text": "change_run_command!(job::DFJob, filenames, run_command)\n\nGoes through the calculation filenames and sets the run command of the calculation.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_atoms-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.get_atoms",
    "category": "method",
    "text": "get_atoms(job::DFJob)\n\nReturns a list the atoms in the structure.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_block-Tuple{DFControl.DFJob,Any,Symbol}",
    "page": "Job",
    "title": "DFControl.get_block",
    "category": "method",
    "text": "get_block(job::DFJob, calc_filenames, block_symbol::Symbol)\n\nLooks through the calculation filenames and returns the block with the specified symbol.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_data-Tuple{DFControl.DFJob,Any,Symbol}",
    "page": "Job",
    "title": "DFControl.get_data",
    "category": "method",
    "text": "get_data(job::DFJob, calc_filenames, block_symbol::Symbol)\n\nLooks through the calculation filenames and returns the data with the specified symbol.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_flag-Tuple{DFControl.DFJob,Any,Symbol}",
    "page": "Job",
    "title": "DFControl.get_flag",
    "category": "method",
    "text": "get_flag(job::DFJob, calc_filenames, flag_name::Symbol)\n\nLooks through the calculation filenames and returns the value of the specified flag.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_flag-Tuple{DFControl.DFJob,Symbol}",
    "page": "Job",
    "title": "DFControl.get_flag",
    "category": "method",
    "text": "get_flag(job::DFJob, flag_name::Symbol)\n\nLooks through all the calculations and returns the value of the specified flag.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_input-Tuple{DFControl.DFJob,Array{String,1}}",
    "page": "Job",
    "title": "DFControl.get_input",
    "category": "method",
    "text": "get_input(job::DFJob, filenames::Array)\n\nReturns an array of the inputs that match one of the filenames.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_input-Tuple{DFControl.DFJob,String}",
    "page": "Job",
    "title": "DFControl.get_input",
    "category": "method",
    "text": "get_input(job::DFJob, filename::String)\n\nReturns the input that matches the filename.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_inputs-Tuple{DFControl.DFJob,Array{T,1} where T}",
    "page": "Job",
    "title": "DFControl.get_inputs",
    "category": "method",
    "text": "get_inputs(job::DFJob, filenames::Array)\n\nReturns an array of the inputs that match one of the filenames.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_inputs-Tuple{DFControl.DFJob,String}",
    "page": "Job",
    "title": "DFControl.get_inputs",
    "category": "method",
    "text": "get_inputs(job::DFJob, filename::String)\n\nReturns an array of the input that matches the filename.\n\n\n\n"
},

{
    "location": "job.html#DFControl.get_run_command-Tuple{DFControl.DFJob,Any}",
    "page": "Job",
    "title": "DFControl.get_run_command",
    "category": "method",
    "text": "get_run_command(job::DFJob, filename)\n\nReturns the run command for the specified calculation.\n\n\n\n"
},

{
    "location": "job.html#DFControl.load_job",
    "page": "Job",
    "title": "DFControl.load_job",
    "category": "function",
    "text": "load_job(job_dir::String, T=Float64; job_fuzzy = \"job\", new_job_name=nothing, new_homedir=nothing, server=get_default_server(),server_dir=\"\")\n\nLoads and returns a DFJob. If local_dir is not specified the job directory will be registered as the local one.\n\n\n\n"
},

{
    "location": "job.html#DFControl.load_server_job-Tuple{String,String}",
    "page": "Job",
    "title": "DFControl.load_server_job",
    "category": "method",
    "text": "load_server_job(server_dir::String, local_dir::String; server=get_default_server(), job_fuzzy=\"*job*\", new_job_name=\"\")\n\nPulls a server job to local directory and then loads it. A fuzzy search for the job file will be performed and the found input files will be pulled.\n\n\n\n"
},

{
    "location": "job.html#DFControl.remove_flags!-Tuple{DFControl.DFJob,Array{#s51,1} where #s51<:AbstractString,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.remove_flags!",
    "category": "method",
    "text": "remove_flags!(job::DFJob, calc_filenames, flags...)\n\nLooks through the calculation filenames and removes the specified flags.\n\n\n\n"
},

{
    "location": "job.html#DFControl.remove_flags!-Tuple{DFControl.DFJob,Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.remove_flags!",
    "category": "method",
    "text": "remove_flags!(job::DFJob, flags...)\n\nLooks through all the calculations and removes the flags.\n\n\n\n"
},

{
    "location": "job.html#DFControl.save_job-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.save_job",
    "category": "method",
    "text": "save_job(job::DFJob)\n\nSaves a DFJob, it\'s job file and all it\'s input files.\n\n\n\n"
},

{
    "location": "job.html#DFControl.set_flags!-Tuple{DFControl.DFJob,Array{String,1},Vararg{Any,N} where N}",
    "page": "Job",
    "title": "DFControl.set_flags!",
    "category": "method",
    "text": "set_flags!(job::DFJob, calculations::Vector{String}, flags...; print=true)\n\nSets the flags in the calculations to the flags specified. This only happens if the specified flags are valid for the calculations. If necessary the correct control block will be added to the calculation (e.g. for QEInputs).\n\nThe values that are supplied will be checked whether they are valid.\n\n\n\n"
},

{
    "location": "job.html#DFControl.set_flow!-Tuple{DFControl.DFJob,Array{Bool,1}}",
    "page": "Job",
    "title": "DFControl.set_flow!",
    "category": "method",
    "text": "set_flow!(job::DFJob, should_runs::Vector{Bool})\n\nSets whether calculations should be ran or not. should_runs should have the same length as the amount of calculations in the job.\n\n\n\n"
},

{
    "location": "job.html#DFControl.set_server_dir!-Tuple{Any,Any}",
    "page": "Job",
    "title": "DFControl.set_server_dir!",
    "category": "method",
    "text": "Sets the server dir of the job.\n\n\n\n"
},

{
    "location": "job.html#DFControl.submit_job-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.submit_job",
    "category": "method",
    "text": "submit_job(job::DFJob; server=nothing, server_dir=nothing)\n\nSubmit a DFJob. First saves it locally, pushes it to the server then runs the job file on the server.\n\n\n\n"
},

{
    "location": "job.html#DFControl.undo!-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.undo!",
    "category": "method",
    "text": "undo!(job::DFJob)\n\nUndos the last change to the calculations of the job.\n\n\n\n"
},

{
    "location": "job.html#DFControl.undo-Tuple{DFControl.DFJob}",
    "page": "Job",
    "title": "DFControl.undo",
    "category": "method",
    "text": "undo(job::DFJob)\n\nUndos the last change to the calculations of the job and returns as a new one.\n\n\n\n"
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
    "location": "input.html#DFControl.QEInput-Tuple{DFControl.QEInput,Any,Vararg{Any,N} where N}",
    "page": "Inputs",
    "title": "DFControl.QEInput",
    "category": "method",
    "text": "QEInput(template::QEInput, filename, newflags...; run_command=template.run_command, run=true, new_data...)\n\nCreates a new QEInput from a template QEInput, setting the newflags in the new one.\n\n\n\n"
},

{
    "location": "input.html#DFControl.change_flags!-Tuple{DFControl.QEInput,Vararg{Any,N} where N}",
    "page": "Inputs",
    "title": "DFControl.change_flags!",
    "category": "method",
    "text": "change_flags!(input::QEInput, new_flag_data...)\n\nChanges the flags inside the input to the new ones if they are already defined and if the new ones have the same type.\n\n\n\n"
},

{
    "location": "input.html#DFControl.change_kpoints!-Tuple{DFControl.QEInput,Array{NTuple{4,#s35} where #s35<:AbstractFloat,1}}",
    "page": "Inputs",
    "title": "DFControl.change_kpoints!",
    "category": "method",
    "text": "change_kpoints!(input::QEInput, k_grid::Vector{NTuple{4, <:AbstractFloat}};\nk_option=:crystal_b)\n\nChanges the data in the k point DataBlock inside the specified calculation. The format is [(ka, kb, kc, nk),...]. This format is to be used with a \'bands\' calculation.\n\n\n\n"
},

{
    "location": "input.html#DFControl.change_kpoints!-Tuple{DFControl.QEInput,Union{NTuple{6,Int64}, Tuple{Int64,Int64,Int64}}}",
    "page": "Inputs",
    "title": "DFControl.change_kpoints!",
    "category": "method",
    "text": "change_kpoints!(input::QEInput, k_grid::Union{NTuple{3, Int}, NTuple{6, Int}})\n\nChanges the data in the k point DataBlock inside the specified calculation. If the specified calculation is \'nscf\' the accepted format is (nka, nkb, nkc), and the k_grid will be generated. If the calculation is \'scf\' the format is (nka, nkb, nkc, sta, stb, stc).\n\n\n\n"
},

{
    "location": "input.html#DFControl.get_block-Tuple{DFControl.QEInput,Symbol}",
    "page": "Inputs",
    "title": "DFControl.get_block",
    "category": "method",
    "text": "get_block(input::QEInput, block_symbol::Symbol)\n\nReturns the block with name block_symbol.\n\n\n\n"
},

{
    "location": "input.html#DFControl.get_flag-Tuple{DFControl.QEInput,Symbol}",
    "page": "Inputs",
    "title": "DFControl.get_flag",
    "category": "method",
    "text": "get_flag(input::QEInput, flag::Symbol)\n\nReturns the value of the flag.\n\n\n\n"
},

{
    "location": "input.html#DFControl.remove_flags!-Tuple{DFControl.QEInput,Vararg{Any,N} where N}",
    "page": "Inputs",
    "title": "DFControl.remove_flags!",
    "category": "method",
    "text": "remove_flags!(input::QEInput, flags...)\n\nRemove the specified flags.\n\n\n\n"
},

{
    "location": "input.html#DFControl.set_flags!-Tuple{DFControl.QEInput,Vararg{Any,N} where N}",
    "page": "Inputs",
    "title": "DFControl.set_flags!",
    "category": "method",
    "text": "set_flags!(input::QEInput, flags...)\n\nSets the specified flags in the input, if they are allowed. The flag values will be converted to the correct type according to the Documentation provided by QE. A ControlBlock will be added to the input if necessary.\n\n\n\n"
},

{
    "location": "input.html#QEInput-1",
    "page": "Inputs",
    "title": "QEInput",
    "category": "section",
    "text": "Modules = [DFControl]\nPages   = [\"qe/input.jl\"]\nPrivate = false\nOrder   = [:type, :function]"
},

{
    "location": "input.html#DFControl.add_block!-Tuple{DFControl.DFInput,DFControl.Block}",
    "page": "Inputs",
    "title": "DFControl.add_block!",
    "category": "method",
    "text": "add_block(input::DFInput, block::Block)\n\nAdds the given block to the input. Should put it in the correct arrays.\n\n\n\n"
},

{
    "location": "input.html#DFControl.add_data!",
    "page": "Inputs",
    "title": "DFControl.add_data!",
    "category": "function",
    "text": "add_data!(input::DFInput, block_name::Symbol, block_data, block_option=:none)\n\nAdds a block with the given name data and option to the calculation.\n\n\n\n"
},

{
    "location": "input.html#DFControl.change_data!-Tuple{DFControl.DFInput,Symbol,Any}",
    "page": "Inputs",
    "title": "DFControl.change_data!",
    "category": "method",
    "text": "change_data!(input::DFInput, block_name::Symbol, new_block_data; option=nothing, print=true)\n\nChanges the data of the specified \'DataBlock\' to the new data. Optionally also changes the \'DataBlock\' option.\n\n\n\n"
},

{
    "location": "input.html#DFControl.change_data_option!-Tuple{DFControl.DFInput,Symbol,Symbol}",
    "page": "Inputs",
    "title": "DFControl.change_data_option!",
    "category": "method",
    "text": "change_data_option!(input::DFInput, block_symbol::Symbol, option::Symbol;; print=true)\n\nChanges the option of specified data block.\n\n\n\n"
},

{
    "location": "input.html#DFControl.change_flags!-Tuple{Union{DFControl.AbinitInput, DFControl.WannierInput},Vararg{Any,N} where N}",
    "page": "Inputs",
    "title": "DFControl.change_flags!",
    "category": "method",
    "text": "change_flags!(input::DFInput, new_flag_data...)\n\nChanges the flags inside the input to the new ones if they are already defined and if the new ones have the same type.\n\n\n\n"
},

{
    "location": "input.html#DFControl.change_kpoints!-Tuple{DFControl.WannierInput,Tuple{Int64,Int64,Int64}}",
    "page": "Inputs",
    "title": "DFControl.change_kpoints!",
    "category": "method",
    "text": "change_kpoints!(input::WannierInput, k_grid::NTuple{3, Int}; print=true)\n\nChanges the data in the k point DataBlock inside the specified calculation.\n\n\n\n"
},

{
    "location": "input.html#DFControl.get_block-Tuple{Union{DFControl.AbinitInput, DFControl.WannierInput},Symbol}",
    "page": "Inputs",
    "title": "DFControl.get_block",
    "category": "method",
    "text": "get_block(input::DFInput, block_symbol::Symbol)\n\nReturns the block with name block_symbol.\n\n\n\n"
},

{
    "location": "input.html#DFControl.get_data-Tuple{DFControl.DFInput,Symbol}",
    "page": "Inputs",
    "title": "DFControl.get_data",
    "category": "method",
    "text": "get_data(input::DFInput, block_symbol::Symbol)\n\nReturns the specified \'DataBlock\'.\n\n\n\n"
},

{
    "location": "input.html#DFControl.get_flag-Tuple{Union{DFControl.AbinitInput, DFControl.WannierInput},Symbol}",
    "page": "Inputs",
    "title": "DFControl.get_flag",
    "category": "method",
    "text": "get_flag(input::DFInput, flag::Symbol)\n\nReturns the value of the flag.\n\n\n\n"
},

{
    "location": "input.html#DFControl.remove_flags!-Tuple{Union{DFControl.AbinitInput, DFControl.WannierInput},Vararg{Any,N} where N}",
    "page": "Inputs",
    "title": "DFControl.remove_flags!",
    "category": "method",
    "text": "remove_flags!(input::DFInput, flags...)\n\nRemove the specified flags.\n\n\n\n"
},

{
    "location": "input.html#DFControl.set_flags!-Tuple{Union{DFControl.AbinitInput, DFControl.WannierInput},Vararg{Any,N} where N}",
    "page": "Inputs",
    "title": "DFControl.set_flags!",
    "category": "method",
    "text": "set_flags!(input::DFInput, flags...)\n\nSets the specified flags in the input. A controlblock will be added if necessary.\n\n\n\n"
},

{
    "location": "input.html#WannierInput-and-AbinitInput-1",
    "page": "Inputs",
    "title": "WannierInput & AbinitInput",
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
    "location": "structure.html#DFControl.change_projections!-Tuple{DFControl.Structure,Vararg{Any,N} where N}",
    "page": "Structure",
    "title": "DFControl.change_projections!",
    "category": "method",
    "text": "Changes the projections of the specified atoms.\n\n\n\n"
},

{
    "location": "structure.html#DFControl.get_atoms-Tuple{DFControl.AbstractStructure,Symbol}",
    "page": "Structure",
    "title": "DFControl.get_atoms",
    "category": "method",
    "text": "Returns all the atoms inside the structure with the specified symbol\n\n\n\n"
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
    "location": "fileio.html#DFControl.read_abi_ebands",
    "page": "FileIO",
    "title": "DFControl.read_abi_ebands",
    "category": "function",
    "text": "Reads an abinit EBANDS.agr output file and returns the found DFBands. K-points only given in crystalline coordinates.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_abi_eig",
    "page": "FileIO",
    "title": "DFControl.read_abi_eig",
    "category": "function",
    "text": "Reads and abinit _EIG output file and returns the found DFBands. K-points are only given in crystalline coordinates.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_abi_fatbands",
    "page": "FileIO",
    "title": "DFControl.read_abi_fatbands",
    "category": "function",
    "text": "Reads an abinit FATBANDS output file and returns the found DFBands, with pdos values in the data field. K-points are not given (they aren\'t present in the output file).\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.read_abi_input",
    "page": "FileIO",
    "title": "DFControl.read_abi_input",
    "category": "function",
    "text": "read_abi_input(filename::String, T=Float64)\n\nReturns an ABINIT input. We assume that jdtset is on a seperate line. If # DATASET is supplied as first line, it will make sure that amount of datasets are read no matter what ndtset and jdtset are. ndtset and jdtset will be taken into account to decide which calculations will be marked as \'should run\'. Either all structures of the input are the same or all are different, otherwise there will be an error.\n\n\n\n"
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
    "text": "read_qe_input(filename, T=Float64; exec=\"pw.x\",  run_command=\"\", run=true, structure_name=\"NoName\")\n\nReads a Quantum Espresso input file. The exec get\'s used to find which flags are allowed in this input file, and convert the read values to the correct Types. Returns a QEInput and the Structure that is found in the input.\n\n\n\n"
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
    "text": "read_qe_output(filename::String, T=Float64)\n\nReads a generic quantum espresso input, returns a dictionary with all found data in the file. Possible keys:\n\n:fermi\n:polarization\n:pol_mod\n:k_cryst\n:k_cart\n:alat\n:cell_parameters\n:pos_option\n:atomic_positions\n:total_force\n:colin_mag_moments\n:bands\n\n\n\n"
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
    "text": "read_wannier_input(filename::String, T=Float64; run_command=\"\", run=true, exec=\"wannier90.x\", structure_name=\"NoName\")\n\nReads a WannierInput and the included Structure from a WANNIER90 input file.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.write_input",
    "page": "FileIO",
    "title": "DFControl.write_input",
    "category": "function",
    "text": "write_input(input::QEInput, structure, filename::String=input.filename)\n\nWrites a Quantum Espresso input file.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.write_input",
    "page": "FileIO",
    "title": "DFControl.write_input",
    "category": "function",
    "text": "write_input(input::WannierInput, structure, filename::String=input.filename)\n\nWrites the WannierInput and structure to a file, that can be interpreted by WANNIER90. The atoms in the structure must have projections defined.\n\n\n\n"
},

{
    "location": "fileio.html#DFControl.write_job_files-Tuple{DFControl.DFJob}",
    "page": "FileIO",
    "title": "DFControl.write_job_files",
    "category": "method",
    "text": "write_job_files(job::DFJob)\n\nWrites all the input files and job file that are linked to a DFJob.\n\n\n\n"
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
    "location": "defaults.html#DFControl.configure_default_pseudos",
    "page": "Defaults",
    "title": "DFControl.configure_default_pseudos",
    "category": "function",
    "text": "configure_default_pseudos(server = get_default_server(), pseudo_dirs=get_default_pseudo_dirs())\n\nReads the specified default_pseudo_dirs on the default_server and sets up the default_pseudos variable, and also adds all the entries to the user_defaults.jl file.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.get_default_pseudo",
    "page": "Defaults",
    "title": "DFControl.get_default_pseudo",
    "category": "function",
    "text": "get_default_pseudo(atom::Symbol, pseudo_set_name=:default; pseudo_specifier=nothing)\n\nReturns the pseudo potential string linked to the atom.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.remove_default_input-Tuple{Symbol}",
    "page": "Defaults",
    "title": "DFControl.remove_default_input",
    "category": "method",
    "text": "remove_default_input(input::Symbol)\n\nRemove input from the default_inputs variable. Also removes the stored input file.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.remove_default_pseudo_dir-Tuple{Symbol}",
    "page": "Defaults",
    "title": "DFControl.remove_default_pseudo_dir",
    "category": "method",
    "text": "remove_default_pseudo_dir(pseudo_symbol::Symbol)\n\nRemoves entry with flag pseudo_symbol from the default_pseudo_dirs and user_defaults.jl file.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.set_default_input-Tuple{DFControl.DFInput,Any,Symbol}",
    "page": "Defaults",
    "title": "DFControl.set_default_input",
    "category": "method",
    "text": "set_default_input(input::dfinput, structure, calculation::Symbol)\n\nAdds the input to the default_inputs variable, and writes it to a file in user_defaults folder to be read every time on load.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.set_default_job_header-Tuple{Any}",
    "page": "Defaults",
    "title": "DFControl.set_default_job_header",
    "category": "method",
    "text": "set_default_job_header(lines)\n\nSets the header that will get added to each job.tt file, if no other header was specified.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.set_default_pseudo_dir-Tuple{Symbol,String}",
    "page": "Defaults",
    "title": "DFControl.set_default_pseudo_dir",
    "category": "method",
    "text": "set_default_pseudo_dir(pseudo_symbol::Symbol, dir::String)\n\nAdds an entry inside the default_pseudo_dirs with flag pseudo_symbol, and adds it to the user_defaults.jl file.\n\n\n\n"
},

{
    "location": "defaults.html#DFControl.set_default_server-Tuple{String}",
    "page": "Defaults",
    "title": "DFControl.set_default_server",
    "category": "method",
    "text": "set_default_server(server::String)\n\nSets the default server variable, and also adds it to the user_defaults.jl file.\n\n\n\n"
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
    "location": "utils.html#DFControl.getfirst-Tuple{Function,Any}",
    "page": "Utils",
    "title": "DFControl.getfirst",
    "category": "method",
    "text": "It\'s like filter()[1].\n\n\n\n"
},

{
    "location": "utils.html#Utils-1",
    "page": "Utils",
    "title": "Utils",
    "category": "section",
    "text": "Modules = [DFControl]\nPages   = [\"utils.jl\"]\nPrivate = false\nOrder   = [:type, :function]"
},

]}
