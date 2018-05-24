__precompile__()
module DFControl
# using Reexport
    import Base.Iterators.flatten
    using RecipesBase
    using StaticArrays
    using GeometryTypes
    using Parameters

    abstract type Package end
    struct Wannier90 <: Package end
    struct QE <: Package end
    struct Abinit <: Package end

    include("atom.jl")
    include("structure.jl")
    include("types.jl")
    include("input.jl")
    include("utils.jl")
    include("job.jl")

    export DFJob
    export save
    export submit

    export flag
    export setflags!
    export rmflags!

    export data
    export setdata!

    export setflow!

    export runcommand
    export setruncommand!

    export setrunflags!
    export runflags
    export setexecflags!
    export inputs
    export input
    export setfilename!
    export setkpoints!
    export setdataoption!
    export setheaderword!
    export setpseudos!
    export addbandscalculation!
    export path
    export undo
    export undo!
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
    export print_info

    include("constants.jl")
    export qe_input_flags

    include("fileio.jl")
    # export read_abi_input
    # export read_abi_output
    # export read_abi_fatbands
    # export read_abi_ebands
    # export read_abi_eig
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
    export removedefault_pseudos
    export setdefault_jobheader
    export @add_default
    export setdefault_input
    export removedefault_input
    export getdefault_pseudo
    export getdefault_server
    export getdefault_jobheader
    export getdefault_pseudodir
    export getdefault_pseudodirs
    export removedefault

    if Pkg.installed("Atom") != nothing
        include("display/printing_juno.jl")
    end
    dfprintln(s::String) = println(s)
    function __init__()
        init_defaults(default_file)
    end


    const UNDO_JOBS = DFJob[]

    const cif2cellpath = joinpath(@__DIR__, "../deps/bin/cif2cell")
end
