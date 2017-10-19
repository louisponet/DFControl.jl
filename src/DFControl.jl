module DFControl
  # using Reexport
  using RecipesBase
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
  export apply_fermi_level
  export apply_fermi_level!

  #@Cleanup this should all be just one thing without qe
  #@Cleanup what do we actually want to have as frontend?
  include("file_processing.jl")
  export read_qe_bands_file
  export read_ks_from_qe_bands_file
  export read_fermi_from_qe_file
  export read_qe_kpdos
  export read_qe_input
  export write_qe_input
  export read_wannier_input
  export write_wannier_input
  export write_df_input
  export write_job_files

  include("job_control.jl")
  export load_job
  export load_qe_job
  export load_qe_server_job
  export pull_job
  export load_server_job
  export save_job
  export push_job
  export submit_job
  export change_job_control_flags!
  export change_input_control_flags!
  export change_input_data!
  export change_job_data!
  export set_input_control_flags!
  export set_job_control_flags!
  export remove_input_control_flag!
  export remove_job_control_flag!
  export remove_job_control_flags!
  export set_should_run!

  include("plotting.jl")
  export plot_qe_bands
  export plot_qe_kpdos

  include("server_comm.jl")
  export read_errors
  export pull_job_outputs
end
