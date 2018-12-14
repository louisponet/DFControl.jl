module DFControl
    using LinearAlgebra
    using Statistics
    using Media
    import Base.Iterators.flatten
    using RecipesBase
    include("FixedSizeArrays.jl")
    using .FixedSizeArrays
    using Parameters

    abstract type Package end
    struct Wannier90 <: Package end
    struct QE <: Package end
    struct Abinit <: Package end

    const depsdir = joinpath(dirname(@__DIR__), "deps")

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
           inpath, outpath, setname!, setkpoints!, setdataoption!, setpseudos!, setcutoffs!,
           atoms, atom, setatoms!, setprojections!, projections,
           addwancalc!, addcalc!,
           setwanenergies!,
           outputdata,

           setheaderword!, setserverdir!, setlocaldir!,
           save, submit, abort, isrunning, progressreport,
           undo, undo!

    include("constants.jl")
    export qe_input_flags

    include("fileio.jl")
    # export read_abi_input
    export read_abi_output
    # export read_abi_fatbands
    export read_abi_ebands
    # export read_abi_eig
    export read_qe_bands_file
    export read_ks_from_qe_output
    export read_fermi_from_qe_output
    export read_qe_kpdos
    export read_qe_pdos
    export read_qe_polarization
    export read_qe_input
    export read_qe_output
    export read_qe_projwfc
    export read_wannier_input
    export read_wannier_output
    export readbands

    include("plotting.jl")

    include("server.jl")
    export pulloutputs
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

    using Crayons
    const dfprintln = println
    const dfprint = print
    include("display/overloads.jl")
    using Requires
    function __init__()
        @require Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d" include("display/printing_juno.jl")
        init_defaults(default_file)
    end

    const pythonpath = Sys.iswindows() ? joinpath(depsdir, "python2", "python") : joinpath(dirname(@__DIR__), "deps", "python2", "bin", "python")
    const cif2cellpath = Sys.iswindows() ? joinpath(depsdir, "python2", "Scripts", "cif2cell") : joinpath(dirname(@__DIR__), "deps", "python2", "bin", "cif2cell")

end
