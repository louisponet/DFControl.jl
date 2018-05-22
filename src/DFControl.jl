__precompile__()
module DFControl
# using Reexport
    import Base.Iterators.flatten
    using RecipesBase
    using StaticArrays
    using GeometryTypes

    include("atom.jl")
    export Atom
    export orbital2atom

    include("structure.jl")
    export AbstractStructure
    export Structure
    export cif2structure

    include("types.jl")
    export Band
    export DFBand
    export Exec

    include("utils.jl")
    export getfirst

    include("input.jl")
    export QEDataBlock
    export WannierDataBlock
    export AbinitDataBlock
    export QEInput
    export AbinitInput
    export WannierInput

    include("job.jl")
    export DFJob
    export save
    export submit

    export setflags!
    export flag
    export setdata!
    export data
    export block
    export setblock!
    export setdata!
    export setflags!
    export rmflags!
    export setflow!
    export setflow!
    export add!
    export runcommand
    export setruncommand!
    export setrunflags!
    export runflags
    export setexecflags!
    export execflags
    export inputs
    export input
    export setfilename!
    export setkpoints!
    export setoption!
    export setheaderword!
    export setpseudos!
    export addbandscalculation!
    export path
    export undo!
    export undo
    export setserverdir!
    export atoms
    export setatoms!
    export cell
    export setcell!
    export path
    export addwancalc!
    export setlocaldir!
    export setserverdir!
    export setprojections!

    include("constants.jl")
    export qe_input_flags

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

    include("plotting.jl")
    export plot_qe_bands
    export plot_qe_kpdos

    include("server_comm.jl")
    export read_errors
    export outputs
    export pull_file
    export pull_files
    export qstat
    export watch_qstat

    include("defaults.jl")
    export setdefault_pseudodir
    export setdefault_server
    export configure_defaultpseudos
    export removedefault_pseudodir
    export default_pseudo
    export setdefault_jobheader
    export @add_default
    export load_defaults
    export setdefault_input
    export removedefault_input

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

    const cif2cellpath = joinpath(@__DIR__, "../deps/bin/cif2cell")
end
