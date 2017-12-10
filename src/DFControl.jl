module DFControl
  # using Reexport
  using RecipesBase
  using Reactive

  const print_s = Signal("hello")
  if Pkg.installed("GtkReactive") != nothing
    # using GLVisualize, GLWindow, GLAbstraction
    include("dashboard.jl")
    # preserve(map(println,print_s))
  else
    preserve(map(println,print_s))
  end
  dfprintln(s::String) = push!(print_s,s)
  Base.display(x) = dfprintln("$(print(x))")
  export dfprintln 

  include("types.jl")
  export Point3D
  export Band
  export DFBand
  export ControlBlock
  export QEControlBlock
  export DataBlock
  export QEDataBlock
  export DFInput
  export QEInput
  export WannierInput
  export DFJob
  export ELEMENTS
  
  include("utils.jl")
  export print_qe_flags
  export print_qe_namelists
  export apply_fermi_level
  export apply_fermi_level!
  export gen_k_grid
  
  #@Cleanup this should all be just one thing without qe
  #@Cleanup what do we actually want to have as frontend?
  include("file_processing.jl")
  export read_qe_bands_file
  export read_ks_from_qe_bands_file
  export read_fermi_from_qe_file
  export read_qe_kpdos
  export read_qe_polarization
  export read_qe_input
  export read_qe_output
  export read_wannier_input
  export write_input
  export write_job_files
  export expr2file
  
  include("input_control.jl")
  include("job_control.jl")
  export create_job
  export load_job
  export pull_job
  export load_server_job
  export save_job
  export push_job
  export submit_job
  
  export change_flags!
  export get_flag
  export change_data!
  export get_data
  export get_block
  export add_flags!
  export remove_flags!
  export change_flow!
  export set_flow!
  export add_calculation!
  export get_run_command
  export change_run_command!
  export get_inputs
  export get_input
  export print_run_command
  export print_flow
  export print_block
  export print_blocks
  export print_data
  export print_filename
  export print_info 
  export print_flags
  export print_flag
  export change_atoms!
  export change_cell_parameters!
  export change_k_points!
  export change_data_option!
  export change_header_word!
  
  include("plotting.jl")
  export plot_qe_bands
  export plot_qe_kpdos
  
  include("server_comm.jl")
  export read_errors
  export pull_outputs
  export pull_file
  export qstat
  export watch_qstat
  
  include("defaults.jl")
  export set_default_pseudo_dir
  export set_default_server
  export configure_default_pseudos
  export remove_default_pseudo_dir
  export set_default_job_header
  export @add_default
  export load_defaults
  export set_default_input
  export remove_default_input

  init_defaults(default_file)
  qe_input_files = search_dir(joinpath(@__DIR__,"../assets/inputs/qe/"),"INPUT")
  const QEControlFlags = vcat([read_qe_flags(joinpath(@__DIR__,"../assets/inputs/qe/") * file) for file in qe_input_files]...)
  export QEControlFlags
end
