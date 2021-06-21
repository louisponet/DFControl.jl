"""
    Orbital(name::String, size::Int, l::Int, mr::Int)

Wannier90 orbital as defined in the [Wannier90 User Guide](https://github.com/wannier-developers/wannier90/raw/v3.1.0/doc/compiled_docs/user_guide.pdf). The `name` will be written in the `projections` block of the Wannier90 input.
"""
struct Orbital
    name::String
    size::Int
    l   ::Int
    mr  ::Int
end
Orbital() = Orbital("none", 0, 0, 0)


"""
    Projection(orb::Orbital, start::Int, last::Int)

A Wannier90 `Projection`, representing an `Orbital` with indices from `start` to `last`.
"""
struct Projection
    orb   ::Orbital
    start ::Int
    last  ::Int
end
Projection() = Projection(Orbital(), 0, 0)


"""
    DFTU(;l ::Int = -1,
          U ::T   = zero(T),
          J0::T  = zero(T),
          α ::T   = zero(T),
          β ::T   = zero(T),
          J ::Vector{T} = T[zero(T)])

DFT+U parameters for a given [`Atom`](@ref).
"""
@with_kw mutable struct DFTU{T}
    l::Int = -1
    U::T   = zero(T)
    J0::T  = zero(T)
    #QE params
    α::T   = zero(T)
    β::T   = zero(T)
    J::Vector{T} = T[zero(T)]
end


"""
    Pseudo(name::String, dir::String)

A pseudo potential file.
"""
mutable struct Pseudo
    name::String
    dir ::String
    Pseudo() =
        new("", "")
    Pseudo(name::AbstractString, dir::AbstractString) =
        new(name, abspath(dir))
end


"""
    Element(symbol::Symbol, Z::Int, name::String, atomic_weight::Float64, color::NTuple{3, Float64})
    
Represents an element. Most conveniently used trough the function [`element`](@ref).
"""
struct Element
    symbol        ::Symbol
    Z             ::Int64
    name          ::String
    atomic_weight ::Float64
    color         ::NTuple{3, Float64}
end


abstract type AbstractAtom{T, LT<:Length{T}} end
# TODO Multiple l per atom in Elk??
#We use angstrom everywhere
"""
    Atom(name::Symbol, element::Element, position_cart::Point3{Length}, position_cryst::Point3;
         pseudo::Pseudo = Pseudo(),
         projections::Vector{Projection} = Projection[],
         magnetization::Vec3 = Vec3(0.0, 0.0, 0.0),
         dftu::DFTU = DFTU())
         
Representation of an `atom`.
    
The `name` of the `atom` is used as an identifier for the `atom` type, in the sense that atoms with the same `pseudo`, `projections`, `magnetization` and [`dftu`](@ref DFTU) attributes should belong to the same type. This also means that during sanity checks atoms that are not of the same type will be given different names. This is done in this way because it often makes sense to change these parameters on all atoms of the same kind at the same time, but still allow the flexibility to change them for individual atoms as well.

`position_cart` should have a valid `Unitful.Length` type such as `Ang`.

See documentation for [`Element`](@ref) and [`Pseudo`](@ref) for further information on these attributes.
"""
@with_kw_noshow mutable struct Atom{T<:AbstractFloat, LT<:Length{T}} <: AbstractAtom{T, LT}
    name          ::Symbol
    element       ::Element
    position_cart ::Point3{LT}
    position_cryst::Point3{T}
    pseudo        ::Pseudo=Pseudo()
    projections   ::Vector{Projection} = Projection[]
    magnetization ::Vec3{T} = zero(Vec3{T})
    dftu          ::DFTU{T} = DFTU{T}()
end
Atom(name::Symbol, el::Element, pos_cart::Point3{LT}, pos_cryst::Point3{T}; kwargs...) where {T, LT<:Length{T}} =
    Atom{T, LT}(name=name, element=el, position_cart=pos_cart, position_cryst=pos_cryst; kwargs...)
Atom(name::Symbol, el::Symbol, args...; kwargs...) =
    Atom(name, element(el), args...; kwargs...)

#TODO this is a little iffy
Atom(orig_at::Atom, new_pos_cart::Point3, new_pos_cryst::Point3) = Atom(name(orig_at), element(orig_at), new_pos_cart, new_pos_cryst, pseudo(orig_at), projections(orig_at), magnetization(orig_at), dftu(orig_at))


abstract type AbstractStructure{T,LT} end
"""
    Structure([name::String,] cell::Mat3, atoms::Vector{<:AbstractAtom}[, data::Dict{Symbol,Any}])

The structure on which the [`DFCalculations`](@ref DFCalculation) will be performed.

    Structure(cif_file::String; name = "NoName")

Creates a [`Structure`](@ref) from the supplied cif file.
"""
mutable struct Structure{T <: AbstractFloat,AA <: AbstractAtom{T},LT <: Length{T}} <: AbstractStructure{T,LT}
    name  ::String
    cell  ::Mat3{LT}
    atoms ::Vector{AA}
    data  ::Dict{Symbol,Any}
end

Structure(name, cell::Mat3{LT}, atoms::Vector{<:Atom{T,LT}}) where {T <: AbstractFloat,LT <: Length{T}} =
    Structure{T,Atom{T,LT},LT}(name, cell, atoms, Dict{Symbol,Any}())

Structure(str::AbstractStructure{T}, atoms::Vector{AT}) where {T <: AbstractFloat,AT <: AbstractAtom{T}} =
    Structure{T,AT}(name(str), cell(str), atoms, data(str))

Structure(cell::Matrix{LT}, atoms::Vector{Atom{T,LT}}) where {T <: AbstractFloat,LT <: Length{T}} =
    Structure{T,Atom{T,LT},LT}("NoName", cell, atoms, Dict{Symbol,Any}())

function Structure(cif_file::String; name = "NoName")
    str = cif2structure(cif_file, structure_name = name)
    @info "Structure extracted from $cif_file\n\tcell parameters: \n\t a = $((str.cell[:,1]...,))\n\t b = $((str.cell[:,2]...,))\n\t c = $((str.cell[:,3]...,))\n\tnat = $(length(str.atoms))\n\telements = $(unique(getfield.(str.atoms, :name)))"
    return str
end


mutable struct ExecFlag
    symbol     ::Symbol
    name       ::String
    typ        ::Type
    description::String
    value
    minus_count::Int
end

ExecFlag(e::ExecFlag, value) = ExecFlag(e.symbol, e.name, e.typ, e.description, value, 1)
ExecFlag(p::Pair{Symbol, T}) where T = ExecFlag(first(p), String(first(p)), T, "", last(p), 1)
ExecFlag(p::Pair{Symbol, T}, count::Int) where T = ExecFlag(first(p), String(first(p)), T, "", last(p), count)

"""
    Exec(;exec::String = "", dir::String = "", flags::Vector{ExecFlag} = ExecFlag[])

Representation of an `executable` that will run the [`DFCalculation`](@ref DFCalculation).
Basically `dir/exec --<flags>` inside a job script.

    Exec(exec::String, dir::String, flags::Pair{Symbol}...)

Will first transform `flags` into a `Vector{ExecFlag}`, and construct the [`Exec`](@ref). 
"""
Parameters.@with_kw mutable struct Exec
    exec ::String = ""
    dir  ::String = ""
    flags::Vector{ExecFlag} = ExecFlag[]
end

Exec(exec::String, dir::String, flags::Pair{Symbol}...) =
    Exec(exec, dir, SymAnyDict(flags))

function Exec(exec::String, dir::String, flags::SymAnyDict)
    _flags = ExecFlag[]
    for (f, v) in flags
        if occursin("mpi", exec)
            mflag = mpi_flag(f)
            @assert mflag !== nothing "$f is not a recognized mpirun flag."
        end            
        push!(_flags, ExecFlag(f => v))
    end
    return Exec(exec, dir, _flags)
end


"""
    InputData(name::Symbol, option::Symbol, data::Any)

Represents a more structured block of input data.
e.g. `InputData(:k_points, :automatic, (6,6,6,1,1,1))`
would be translated for a QE calculation into
```
K_POINTS(automatic)
6 6 6 1 1 1
```
"""
mutable struct InputData
    name   ::Symbol
    option ::Symbol
    data   ::Any
end


"""
    DFCalculation{P<:Package}(name    ::String;
                              dir     ::String = "",
                              flags   ::AbstractDict = Dict{Symbol, Any}(),
                              data    ::Vector{InputData} = InputData[],
                              execs   ::Vector{Exec},
                              run     ::Bool = true,
                              outdata ::AbstractDict = Dict{Symbol, Any}(),
                              infile  ::String = P == Wannier90 ? name * ".win" : name * ".in",
                              outfile ::String = name * ".out")

The representation of a *DFT* calculation of package `P`,
holding the `flags` that will be written to the `infile`,
the executables in `execs` and the output written by the calculation to the `outfile`.
It essentially represents a line in a job script similar to `exec1 exec2 < infile.in > outfile.out`. 
`outdata` stores the parsed calculation output after it was read at least once.
The `run` field indicates whether the calculation should be actually performed,
e.g. if `run=false` the corresponding line will be commented out in the job script.

    DFCalculation{P<:Package}(name::AbstractString, flags::Pair{Symbol, Any}...; kwargs...)

Create a [`DFCalculation`](@ref) from `name` and `flags`, other `kwargs...` will be passed to the constructor.

    DFCalculation(template::DFCalculation, name::AbstractString, flags::Pair{Symbol, Any}...; excs=deepcopy(execs(template)), run=true, data=nothing, dir=copy(template.dir))

Creates a new [`DFCalculation`](@ref) from the `template`, setting the `flags` of the newly created one to the specified ones.
"""
@with_kw_noshow mutable struct DFCalculation{P <: Package}
    name     ::String
    dir      ::String = ""
    flags    ::AbstractDict = SymAnyDict()
    data     ::Vector{InputData} = InputData[]
    execs    ::Vector{Exec}
    run      ::Bool = true
    outdata  ::SymAnyDict=SymAnyDict()
    infile   ::String = P == Wannier90 ? name * ".win" : name * ".in"
    outfile  ::String = name * ".out"
    function DFCalculation{P}(name, dir, flags, data, execs, run, outdata, infile, outfile) where {P <: Package}
        out = new{P}(name, dir, SymAnyDict(), data, execs, run, outdata, infile, outfile)
        set_flags!(out, flags...; print=false)
        for (f, v) in flags
            if !hasflag(out, f)
                @warn "Flag $f was not valid for calculation $name."
            end
        end
        return out
    end
end
DFCalculation{P}(name, dir, flags, data, execs, run) where P<:Package = DFCalculation{P}(name, abspath(dir), flags, data, execs, run, SymAnyDict(), P == Wannier90 ? name * ".win" : name * ".in", P == Wannier90 ? name * ".wout" : name * ".out")

DFCalculation{P}(name, flags...; kwargs...) where P<:Package =
    DFCalculation{P}(name=name, flags=flags; kwargs...)

function DFCalculation(template::DFCalculation, name, newflags...;
                 excs = deepcopy(execs(template)),
                 run  = true,
                 data = nothing,
                 dir  = deepcopy(template.dir))
    newflags = Dict(newflags...)

    calculation          = deepcopy(template)
    calculation.name     = name
    calculation.execs    = excs
    calculation.run      = run
    calculation.dir      = dir
    set_flags!(calculation, newflags..., print=false)

    if data != nothing
        for (name, (option, data)) in data
            set_data!(calculation, name, data, option=option, print=false)
        end
    end
    return calculation
end


#TODO should we also create a config file for each job with stuff like server etc? and other config things,
#      which if not supplied could contain the default stuff?
"""
    DFJob(name::String, structure::AbstractStructure;
          calculations      ::Vector{DFCalculation} = DFCalculation[],
          local_dir         ::String = pwd(),
          header            ::Vector{String} = getdefault_jobheader(),
          metadata          ::Dict = Dict(),
          version           ::Int = last_job_version(local_dir),
          copy_temp_folders ::Bool = false, 
          server            ::String = getdefault_server(),
          server_dir        ::String = "")

A [`DFJob`](@ref) embodies a set of [`DFCalculations`](@ref DFCalculation) to be ran in directory `local_dir`, with the [`Structure`](@ref) as the subject.
## Keywords/further attributes
- `calculations`: calculations to calculations that will be run sequentially.
- `local_dir`: the directory where the calculations will be run.
- `header`: lines that will be pasted at the head of the job script, e.g. exports `export OMP_NUM_THREADS=1`, slurm settings`#SBATCH`, etc.
- `metadata`: various additional information, will be saved in `.metadata.jld2` in the `local_dir`.
- `version`: the current version of the job.
- `copy_temp_folders`: whether or not the temporary directory associated with intermediate calculation results should be copied when storing a job version. *CAUTION* These can be quite large.
- The `server` and `server_dir` keywords should be avoided for the time being, as this functionality is not well tested.
 
    DFJob(job_name::String, structure::AbstractStructure, calculations::Vector{<:DFCalculation}, common_flags::Pair{Symbol, <:Any}...; kwargs...)

Creates a new job. The common flags will be attempted to be set in each of the `calculations`. The `kwargs...` are passed to the [`DFJob`](@ref) constructor. 

    DFJob(job_dir::String, job_script="job.tt"; version=nothing, kwargs...)

Loads the job in the `local_dir`.
If `job_dir` is not a valid job path, the previously saved jobs will be scanned for a job with a `local_dir` that
partly includes `job_dir`. If `version` is specified the corresponding job version will be returned if it exists. 
The `kwargs...` will be passed to the [`DFJob`](@ref) constructor.
"""
@with_kw_noshow mutable struct DFJob
    name         ::String
    structure    ::AbstractStructure
    calculations ::Vector{DFCalculation} = DFCalculation[]
    local_dir    ::String = pwd()
    header       ::Vector{String} = getdefault_jobheader()
    metadata     ::Dict = Dict()
    version      ::Int = last_job_version(local_dir)
    copy_temp_folders::Bool = false
    server       ::String = getdefault_server()
    server_dir   ::String = ""
    function DFJob(name, structure, calculations, local_dir, header, metadata, version, copy_temp_folders, server, server_dir)
        if local_dir[end] == '/'
            local_dir = local_dir[1:end-1]
        end
        if !isabspath(local_dir)
            local_dir = abspath(local_dir)
        end
        if isempty(structure.name)
            structure.name = split(name, "_")[1]
        end
        last_version = last_job_version(local_dir)
        if isempty(metadata)
            mpath = joinpath(local_dir, ".metadata.jld2")
            if ispath(mpath)
                stored_data = load(mpath)
                metadata = haskey(stored_data, "metadata") ? stored_data["metadata"] : metadata
                version = haskey(stored_data, "version") ? stored_data["version"] : version
            end
        end
        out = new(name, structure, calculations, local_dir, header, metadata, version, copy_temp_folders, server, server_dir)
        return out
    end
end

#TODO implement abinit
function DFJob(job_name::String, structure::AbstractStructure, calculations::Vector{<:DFCalculation}, common_flags::Pair{Symbol, <:Any}...;
                    kwargs...)

    shared_flags = typeof(common_flags) <: Dict ? common_flags : Dict(common_flags...)
    for i in calculations
        i.flags = merge(shared_flags, i.flags)
    end
    out = DFJob(name = job_name, structure = structure, calculations = calculations; kwargs...)
    return out
end

function DFJob(job_dir::String, job_script="job.tt"; version = nothing, kwargs...)
    if !isempty(job_dir) && ispath(abspath(job_dir)) && !isempty(searchdir(abspath(job_dir), job_script))
        real_path = abspath(job_dir)
    else
        real_path = request_job(job_dir)
        real_path === nothing && return
    end
    real_version = version === nothing ? last_job_version(real_path) : version
    return DFJob(;merge(merge((local_dir=real_path,version=real_version), read_job_calculations(joinpath(real_path,  job_script))), kwargs)...)

end        


abstract type Band end

"""
Energy band from DFT calculation.
"""
mutable struct DFBand{T<:AbstractFloat,K} <: Band
    k_points_cart  ::Vector{Vec3{K}}
    k_points_cryst ::Vector{Vec3{T}}
    eigvals        ::Vector{T}
    extra          ::Dict{Symbol, Any}
end
DFBand(k_points_cart, k_points_cryst, eigvals) = DFBand(k_points_cart, k_points_cryst, eigvals, Dict{Symbol,Any}())
DFBand(::Type{T}, vlength::Int) where T = DFBand(Vector{Vec3{T}}(undef, vlength), Vector{Vec3{T}}(undef, vlength), Vector{T}(undef, vlength), SymAnyDict())
DFBand(vlength::Int) = DFBand(Float64, vlength)

kpoints(band::DFBand, kind=:cryst) = kind == :cart ? band.k_points_cart : band.k_points_cryst
eigvals(band::DFBand) = band.eigvals

"""
    bandgap(bands::AbstractVector{DFBand}, fermi=0.0)

Calculates the bandgap (possibly indirect) around the fermi level.
"""
function bandgap(bands::Union{Iterators.Flatten, AbstractVector{<:Band}}, fermi=0.0)
    max_valence = -Inf
    min_conduction = Inf
    for b in bands
        max = maximum(eigvals(b).-fermi)
        min = minimum(eigvals(b).-fermi)
        if max_valence <= max <= 0.0
            max_valence = max
        end
        if 0.0 <= min <= min_conduction
            min_conduction = min
        end
    end
    return min_conduction - max_valence
end
bandgap(u_d_bands::Union{NamedTuple, Tuple}, args...) = bandgap(Iterators.flatten(u_d_bands), args...)

mutable struct TimingData
    name::String
    cpu::Dates.AbstractTime
    wall::Dates.AbstractTime
    calls::Int
    children::Vector{TimingData}
end
