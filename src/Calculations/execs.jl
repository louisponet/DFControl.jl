mutable struct ExecFlag{T}
    symbol      :: Symbol
    name        :: String
    description :: String
    value       :: T
    minus_count :: Int
end
ExecFlag(e::ExecFlag, value) = ExecFlag(e.symbol, e.name, e.description, value, e.minus_count)
function ExecFlag(p::Pair)
    return ExecFlag(first(p), String(first(p)), "", last(p), 1)
end
function ExecFlag(p::Pair, count::Int)
    return ExecFlag(first(p), String(first(p)), "", last(p), count)
end

StructTypes.StructType(::Type{<:ExecFlag}) = StructTypes.Struct()

Base.hash(e::ExecFlag, h::UInt)        = hash(e.symbol, hash(e.value, h))
Base.:(==)(e1::ExecFlag, e2::ExecFlag) = e1.symbol == e2.symbol && e1.value == e2.value
Base.eltype(::ExecFlag{T}) where {T}   = T

"""
    Exec(;exec::String = "", dir::String = "", flags::Vector{ExecFlag} = ExecFlag[])

Representation of an `executable` that will run the [`Calculation`](@ref Calculation).
Basically `dir/exec --<flags>` inside a job script.

    Exec(exec::String, dir::String, flags::Pair{Symbol}...)

Will first transform `flags` into a `Vector{ExecFlag}`, and construct the [`Exec`](@ref). 
"""
@with_kw mutable struct Exec
    exec::String = ""
    dir::String = ""
    flags::Vector{ExecFlag} = ExecFlag[]
end

function Exec(exec::String, dir::String, flags...)
    _flags = ExecFlag[]
    ismpi = occursin("mpi", exec)
    for (f, v) in flags
        if ismpi 
            mflag = mpi_flag(f)
            @assert mflag !== nothing "$f is not a recognized mpirun flag."
        end
        push!(_flags, ExecFlag(f => v))
    end
    return Exec(exec, dir, _flags)
end

StructTypes.StructType(::Type{Exec}) = StructTypes.Mutable()

const RUN_EXECS = ["mpirun", "mpiexec", "srun"]

hasflag(exec::Exec, s::Symbol) = findfirst(x -> x.symbol == s, exec.flags) != nothing

path(exec::Exec)                       = joinpath(exec.dir, exec.exec)

function Base.:(==)(e1::Exec, e2::Exec)
    if e1.exec != e2.exec || e1.dir != e2.dir
        return false
    else
        for f in e1.flags
            f2 = getfirst(isequal(f), e2.flags)
            f2 === nothing && return false
            f == f2 && return true
        end
        return true
    end
end

function Base.string(e::Exec)
    direxec = joinpath(e.dir, e.exec)
    str = "$direxec"
    for flag in e.flags
        str *= " $(join(fill('-', flag.minus_count)))$(flag.symbol)"
        if !isa(flag.value, AbstractString)
            for v in flag.value
                str *= " $v"
            end
        else
            str *= " $(flag.value)"
        end
    end
    return str
end

function set_flags!(exec::Exec, flags...)
    for (f, val) in flags
        flag = isa(f, String) ? getfirst(x -> x.name == f, exec.flags) :
               getfirst(x -> x.symbol == f, exec.flags)
        if flag != nothing
            flag.value = convert(eltype(flag), val)
        else
            found = false
            #TODO generalize this
            for (f1, f2) in zip((is_qe_exec, is_wannier_exec, is_mpi_exec),
                    (qe_execflag, wan_execflag, mpi_flag))

                if f1(exec)
                    def_flag = f2(f)
                    if def_flag !== nothing
                        push!(exec.flags, ExecFlag(def_flag, val))
                        found = true
                        break
                    end
                end
            end
            if !found
                error("Flag $f was not found in allowed executable flags.")
            end
        end
    end
    return exec.flags
end

function rm_flags!(exec::Exec, flags...)
    for f in flags
        if isa(f, String)
            filter!(x -> x.name != f, exec.flags)
        elseif isa(f, Symbol)
            filter!(x -> x.symbol != f, exec.flags)
        end
    end
    return exec.flags
end

#### MPI Exec Functionality ###
include(joinpath(depsdir, "mpirunflags.jl"))
const MPI_FLAGS = _MPI_FLAGS()

mpi_flag(flag::AbstractString) = getfirst(x -> x.name == flag, MPI_FLAGS)
mpi_flag(flag::Symbol) = getfirst(x -> x.symbol == flag, MPI_FLAGS)

function mpi_flag_val(::Type{String}, line, i)
    v = line[i+1]
    return v, i + 2
end

function mpi_flag_val(::Type{Vector{String}}, line, i)
    tval = String[]
    while i + 1 <= length(line) && !occursin('-', line[i+1])
        push!(tval, line[i+1])
        i += 1
    end
    return tval, i + 1
end

function mpi_flag_val(::Type{Int}, line, i)
    if line[i+1][1] == '$'
        tval = line[i+1]
    else
        tval = parse(Int, line[i+1])
    end
    return tval, i + 2
end

function mpi_flag_val(::Type{Vector{Int}}, line, i)
    tval = Union{Int,String}[]
    while i + 1 <= length(line) && !occursin('-', line[i+1])
        if line[i+1][1] == '$'
            push!(tval, line[i+1])
        else
            push!(tval, parse(Int, line[i+1]))
        end
        i += 1
    end
    return tval, i + 1
end

function mpi_flag_val(::Type{Pair{Int,String}}, line, i)
    tval = Union{Int,String}[]
    while i + 1 <= length(line) && !occursin('-', line[i+1])
        if line[i+1][1] == '$'
            push!(tval, line[i+1])
        else
            tparse = tryparse(Int, line[i+1])
            if tparse != nothing
                push!(tval, tparse)
            else
                push!(tval, string(line[i+1]))
            end
        end
        i += 1
    end
    return tval, i + 1
end
function parse_mpi_flags(line::Vector{<:SubString})
    eflags = ExecFlag[]
    i = 1
    while i <= length(line)
        s = line[i]
        if s[1] != '-'
            break
        end
        if s[2] == '-'
            mflag = mpi_flag(strip(s, '-'))
        else
            mflag = mpi_flag(Symbol(strip(s, '-')))
        end
        @assert mflag !== nothing "$(strip(s, '-')) is not a recognized mpi_flag"
        val, i = mpi_flag_val(eltype(mflag), line, i)
        push!(eflags, ExecFlag(mflag, val))
    end
    return eflags
end

is_mpi_exec(exec::Exec) = occursin("mpi", exec.exec)


#Wannier
const WAN_EXECS = ["wannier90.x"]

const WAN_EXECFLAGS = ExecFlag[ExecFlag(:pp, "preprocess",
                                        "Whether or not to preprocess the wannier calculation",
                                        false, 1),]

wan_execflag(flag::AbstractString) = getfirst(x -> x.name == flag, WAN_EXECFLAGS)
wan_execflag(flag::Symbol) = getfirst(x -> x.symbol == flag, WAN_EXECFLAGS)
function parse_wan_execflags(line::Vector{<:AbstractString})
    flags = ExecFlag[]
    i = 1
    while i <= length(line)
        s = strip(line[i], '-')
        push!(flags, ExecFlag(wan_execflag(Symbol(s)), nothing))
        i += 1
    end
    return flags
end

is_wannier_exec(exec::Exec) = exec.exec ∈ WAN_EXECS


#QE
const QE_EXECS = ["pw.x", "projwfc.x", "pp.x", "ld1.x", "ph.x", "pw2wannier90.x", "hp.x"]

const QE_EXECFLAGS = ExecFlag[ExecFlag(:nk, "kpoint-pools", 
                                       "groups k-point parallelization into nk processor pools",
                                       0, 1),
                              ExecFlag(:ntg, "task-groups",  "FFT task groups", 0, 1),
                              ExecFlag(:ndiag, "diag", 
                                       "Number of processes for linear algebra", 0, 1),
                              ExecFlag(:ni, "images", 
                                       "Number of processes used for the images", 0, 1)]

qe_execflag(flag::AbstractString) = getfirst(x -> x.name == flag, QE_EXECFLAGS)
qe_execflag(flag::Symbol) = getfirst(x -> x.symbol == flag, QE_EXECFLAGS)

function parse_qe_execflags(line::Vector{<:AbstractString})
    flags = ExecFlag[]
    i = 1
    while i <= length(line)
        s = strip(line[i], '-')
        push!(flags, ExecFlag(qe_execflag(Symbol(s)), parse(Int, line[i+1])))
        i += 2
    end
    return flags
end

is_qe_exec(exec::Exec) = exec.exec ∈ QE_EXECS

# Elk
const ELK_EXECS = ["elk", "elk-omp"]

is_elk_exec(exec::Exec) = exec.exec ∈ ELK_EXECS