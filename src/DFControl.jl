__precompile__()
module DFControl
# using Reexport
    import Base.Iterators.flatten
    using RecipesBase, GeometryTypes, Parameters, Nullables

    abstract type Package end
    struct Wannier90 <: Package end
    struct QE <: Package end
    struct Abinit <: Package end

    include("atom.jl")
    include("structure.jl")
    include("types.jl")
    export kpoints

    include("input.jl")
    include("utils.jl")
    export kgrid
    include("job.jl")

    export DFJob, Exec, DFInput
    export flag, setflags!, rmflags!,
           data, setdata!,
           cell, setcell!,
           setflow!,
           execs, setexecflags!, setexecdir!, rmexecflags!,
           input, inputs,
           path, outpath, setfilename!, setkpoints!, setdataoption!, setpseudos!,
           atoms, atom, setatoms!, setprojections!, projections,
           addwancalc!, addcalc!,
           setwanenergies!,
           outputdata,

           setheaderword!, setserverdir!, setlocaldir!,
           save, submit,
           undo, undo!,
           print_info

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

    include("server_comm.jl")
    export outputs
    export pullfile
    export pullfiles
    export qstat
    export watch_qstat

    include("defaults.jl")
    export setdefault_pseudodir
    export setdefault_server
    export configuredefault_pseudos
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
    include("display/overloads.jl")
    if Pkg.installed("Atom") != nothing
        include("display/printing_juno.jl")
    end
    dfprintln(io::IO, s::String) = println(io, s)
    dfprintln(s::String) = println(s)
    function __init__()
        init_defaults(default_file)
    end


    const UNDO_JOBS = DFJob[]

    const cif2cellpath = joinpath(@__DIR__, "../deps/bin/cif2cell")
end
