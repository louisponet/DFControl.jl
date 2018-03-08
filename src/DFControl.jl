__precompile__(true)
module DFControl
# using Reexport
    using RecipesBase
    using StaticArrays
    using GeometryTypes
    include("types.jl")
    export element
    export AbstractAtom
    export Atom
    export AbstractStructure
    export Structure

    export Band
    export DFBand

    include("utils.jl")

    include("input.jl")
    export QEDataBlock
    export WannierDataBlock
    export AbinitDataBlock
    export QEInput
    export AbinitInput
    export WannierInput

    include("job.jl")
    export DFJob
    export load_job
    export load_server_job
    export save_job
    export submit_job

    export change_flags!
    export get_flag
    export change_data!
    export get_data
    export get_block
    export add_block!
    export add_data!
    export set_flags!
    export remove_flags!
    export change_flow!
    export set_flow!
    export add_calculation!
    export get_run_command
    export change_run_command!
    export get_inputs
    export get_input
    export change_filename!
    export change_kpoints!
    export change_data_option!
    export change_header_word!
    export change_pseudo_set!
    export add_bands_calculation!
    export get_path
    export undo!
    export undo
    export set_server_dir!
    export get_atoms
    export change_atoms!
    export get_cell
    export change_cell!
    export get_path
    export add_wan_calc!
    export change_local_dir!
    export change_server_dir!
    export change_projections!

    include("constants.jl")

    include("fileio.jl")
    export read_abi_input
    export read_abi_output
    export read_abi_fatbands
    export read_abi_ebands
    export read_abi_eig
    export read_qe_bands_file
    export read_ks_from_qe_output
    export read_fermi_from_qe_output
    export read_qe_kpdos
    export read_qe_pdos
    export read_qe_polarization
    export read_qe_input
    export read_qe_output
    export read_wannier_input
    export write_input
    export write_job_files

    include("plotting.jl")
    export plot_qe_bands
    export plot_qe_kpdos

    include("server_comm.jl")
    export read_errors
    export pull_outputs
    export pull_file
    export pull_files
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

    #no extra functionality, for faster scripting
    include("shortnames.jl")

    if Pkg.installed("Atom") != nothing
        include("display/printing_juno.jl")
    end
    dfprintln(s::String) = println(s)
    function __init__()
        init_defaults(default_file)
    end

    include("display/printing.jl")
    export print_run_command
    export print_flow
    export print_block
    export print_blocks
    export print_data
    export print_filename
    export print_info
    export print_flags
    export print_flag
    export dfprintln
    const UNDO_JOBS = DFJob[]

end
