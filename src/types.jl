const SymAnyDict = Dict{Symbol, Any}
const Vec{N, T} = SVector{N, T}
const Point = Vec
const Point3 = Point{3}
const Vec3   = Vec{3}
const Mat3{T} = SMatrix{3, 3, T}
const Mat4{T} = SMatrix{4, 4, T}

Point{N, T}(x::T) where {N, T} = Point{N, T}(x, x, x)

Base.convert(::Type{Point3{T}}, x::Vector{T}) where T<:AbstractFloat = Point3{T}(x[1], x[2], x[3])

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

mutable struct ExecFlag
    symbol     ::Symbol
    name       ::String
    typ        ::Type
    description::String
    value
end

ExecFlag(e::ExecFlag, value) = ExecFlag(e.symbol, e.name, e.typ, e.description, value)
ExecFlag(p::Pair{Symbol, T}) where T = ExecFlag(first(p), String(first(p)), T, "", last(p))

const QE_EXECFLAGS = ExecFlag[
    ExecFlag(:nk, "kpoint-pools", Int, "groups k-point parallelization into nk processor pools", 0),
    ExecFlag(:ntg, "task-groups", Int, "FFT task groups", 0),
    ExecFlag(:ndiag, "diag", Int, "Number of processes for linear algebra", 0)
]

qeexecflag(flag::AbstractString) = getfirst(x -> x.name==flag, QE_EXECFLAGS)
qeexecflag(flag::Symbol) = getfirst(x -> x.symbol==flag, QE_EXECFLAGS)

function parse_qeexecflags(line::Vector{<:AbstractString})
    flags = ExecFlag[]
    i=1
    while i<=length(line)
        s = strip(line[i], '-')
        push!(flags, ExecFlag(qeexecflag(Symbol(s)), parse(Int, line[i+1])))
        i += 2
    end
    flags
end

const WAN_EXECFLAGS = ExecFlag[
    ExecFlag(:pp, "preprocess", Nothing, "Whether or not to preprocess the wannier input", nothing),
]

wan_execflag(flag::AbstractString) = getfirst(x -> x.name==flag, WAN_EXECFLAGS)
wan_execflag(flag::Symbol) = getfirst(x -> x.symbol==flag, WAN_EXECFLAGS)
function parse_wan_execflags(line::Vector{<:AbstractString})
    flags = ExecFlag[]
    i=1
    while i<=length(line)
        s = strip(line[i], '-')
        push!(flags, ExecFlag(wan_execflag(Symbol(s)), nothing))
        i += 1
    end
    flags
end

include(joinpath(depsdir, "mpirunflags.jl"))
const MPI_FLAGS = _MPI_FLAGS()

mpi_flag(flag::AbstractString) = getfirst(x -> x.name==flag, MPI_FLAGS)
mpi_flag(flag::Symbol) = getfirst(x -> x.symbol==flag, MPI_FLAGS)

function mpi_flag_val(::Type{String}, line, i)
    v = line[i+1]
    return v, i + 2
end

function mpi_flag_val(::Type{Vector{String}}, line, i)
    tval = String[]
    while i+1 <= length(line) && !occursin('-', line[i+1])
        push!(tval, line[i+1])
        i += 1
    end
    return tval, i + 1
end

function mpi_flag_val(::Type{Int}, line, i)
    if line[i+1][1] == '\$'
        tval = line[i+1]
    else
        tval = parse(Int, line[i+1])
    end
    return tval, i + 2
end

function mpi_flag_val(::Type{Vector{Int}}, line, i)
    tval = Union{Int, String}[]
    while i+1 <= length(line) && !occursin('-', line[i+1])
        if line[i+1][1] == '\$'
            push!(tval, line[i+1])
        else
            push!(tval, parse(Int, line[i+1]))
        end
        i += 1
    end
    return tval, i + 1
end

function mpi_flag_val(::Type{Pair{Int, String}}, line, i)
    tval = Union{Int, String}[]
    while i+1<=length(line) && !occursin('-', line[i+1])
        if line[i+1][1] == '\$'
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

        @assert mflag != nothing "$(strip(s, '-')) is not a recognized mpi_flag"
        val, i = mpi_flag_val(mflag.typ, line, i)
        push!(eflags, ExecFlag(mflag, val))
    end
    eflags
end

const RUN_EXECS = ["mpirun"]
allexecs() = vcat(RUN_EXECS, QE_EXECS, WAN_EXECS, ELK_EXECS)
parseable_execs() = vcat(QE_EXECS, WAN_EXECS, ELK_EXECS)
has_parseable_exec(l::String) = occursin(">", l) && any(occursin.(parseable_execs(), (l,)))

mutable struct Exec
    exec ::String
    dir  ::String
    flags::Vector{ExecFlag}
end

Exec() = Exec("")
Exec(exec::String) = Exec(exec, "")
Exec(exec::String, dir::String) = Exec(exec, dir, ExecFlag[])
Exec(exec::String, dir::String, flags...) = Exec(exec, dir, SymAnyDict(flags))
function Exec(exec::String, dir::String, flags::SymAnyDict)
    _flags = ExecFlag[]
    for (f, v) in flags
        push!(_flags, ExecFlag(f => v))
    end
    return Exec(exec, dir, _flags)
end

isparseable(exec::Exec) = exec.exec ∈ parseable_execs()

is_qe_exec(exec::Exec)      = exec.exec ∈ QE_EXECS
is_wannier_exec(exec::Exec) = exec.exec ∈ WAN_EXECS
is_elk_exec(exec::Exec)     = exec.exec ∈ ELK_EXECS
is_mpi_exec(exec::Exec)     = exec.exec == "mpirun"

function inputparser(exec::Exec)
    if is_qe_exec(exec) 
        qe_read_input
    elseif is_wannier_exec(exec)
        wan_read_input
    elseif is_elk_exec(exec)
        elk_read_input
    end
end

function setflags!(exec::Exec, flags...)
    for (f, val) in flags
        flag = isa(f, String) ? getfirst(x -> x.name == f, exec.flags) : getfirst(x -> x.symbol == f, exec.flags)
        if flag != nothing
            flag.value = convert(flag.typ, val)
        else
	        for (f1, f2) in zip((is_qe_exec, is_wannier_exec, is_mpi_exec), (qeexecflag, wan_execflag, mpi_flag))
	        	if f1(exec)
		        	def_flag = f2(f)
			        if def_flag != nothing
				        push!(exec.flags, ExecFlag(def_flag, val))
				        break
			        end
		        end
	        end
        end
    end
    exec.flags
end

function rmflags!(exec::Exec, flags...)
    for f in flags
        if isa(f, String)
            filter!(x -> x.name != f, exec.flags)
        elseif isa(f, Symbol)
            filter!(x -> x.symbol != f, exec.flags)
        end
    end
    exec.flags
end
hasflag(exec::Exec, s::Symbol) = findfirst(x->x.symbol == s, exec.flags) != nothing
setexecdir!(exec::Exec, dir) = exec.dir = dir
