module DFControl
  using Reexport
  @reexport using Plots
  include("types.jl")
  export Point3D
  export Band
  export DFBand
  export DFInput
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
  export write_df_input
  export write_job_files

  include("job_control.jl")
  export load_qe_job
  export load_qe_server_job
  export pull_job
  export save_job
  export push_job
  export submit_job
  export change_job_data!
  export set_job_data!

  include("plotting.jl")
  export plot_qe_bands
  export plot_qe_kpdos

end