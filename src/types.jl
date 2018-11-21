const SymAnyDict = Dict{Symbol, Any}

Base.convert(::Type{Point3{T}}, x::Vector{T}) where T<:AbstractFloat = Point3{T}(x[1], x[2], x[3])


abstract type Band end

"""
Energy band from DFT calculation.
"""
mutable struct DFBand{T<:AbstractFloat} <: Band
    k_points_cart  ::Vector{Vec3{T}}
    k_points_cryst ::Vector{Vec3{T}}
    eigvals        ::Vector{T}
    extra          ::Dict{Symbol, Any}
end
DFBand(k_points_cart, k_points_cryst, eigvals) = DFBand(k_points_cart, k_points_cryst, eigvals, Dict{Symbol,Any}())

kpoints(band::DFBand, kind=:cryst) = kind == :cart ? band.k_points_cart : band.k_points_cryst


mutable struct ExecFlag
    symbol     ::Symbol
    name       ::String
    type       ::Type
    description::String
    value
end

ExecFlag(e::ExecFlag, value) = ExecFlag(e.symbol, e.name, e.type, e.description, value)
ExecFlag(p::Pair{Symbol, T}) where T = ExecFlag(first(p), String(first(p)), T, "", last(p))

const QEFLAGS = ExecFlag[
    ExecFlag(:nk, "kpoint-pools", Int, "groups k-point parallelization into nk processor pools", 0),
    ExecFlag(:ntg, "task-groups", Int, "FFT task groups", 0)
]

qeflag(flag::AbstractString) = getfirst(x -> x.name==flag, QEFLAGS)
qeflag(flag::Symbol) = getfirst(x -> x.symbol==flag, QEFLAGS)
function parse_qeflags(line::Vector{<:AbstractString})
    flags = ExecFlag[]
    i=1
    while i<=length(line)
        s = strip(line[i], '-')
        push!(flags, ExecFlag(qeflag(Symbol(s)), parse(Int, line[i+1])))
        i += 2
    end
    flags
end

const MPIFLAGS = ExecFlag[]

function parse_manfile(file)
    flags = ExecFlag[]
    open(file, "r") do f
        line = readline(f)
        while line != "OPTIONS"
            line = readline(f)
        end
        while line != "Environment Variables"
            line = strip(readline(f))
            if !isempty(line) && line[1] == '-'
                name        = ""     #--
                symbols     = Symbol[] #-
                type        = Nothing
                description = ""
                sline = strip.(split(line), ',')
                for s in sline
                    if s[2] == '-' #--npernode
                        name = strip(s, '-')
                    elseif occursin('<', s) # <#persocket>
                        type = if occursin("#", s)
                                   occursin(',', s) ? Vector{Int} : Int
                               elseif occursin("ppr", s) #ppr:N:<object>
                                   Pair{Int, String}
                               else
                                   occursin(',', s) ? Vector{String} : String
                               end
                        break
                    else  #-np
                        push!(symbols, Symbol(strip(s, '-')))
                    end
                end
                line = strip(readline(f))
                while !isempty(line)
                    description *= " " * line
                    line = strip(readline(f))
                end
                if name != "" && isempty(symbols)
                    symbols = [Symbol(name)]
                end
                for symbol in symbols
                    push!(flags, ExecFlag(symbol, name, type, strip(description), nothing))
                end
            end
        end
    end
    return flags
end

init_mpiflags() = append!(MPIFLAGS, parse_manfile(joinpath(@__DIR__, "..", "assets","mpirun_man.txt")))
mpiflag(flag::AbstractString) = getfirst(x -> x.name==flag, MPIFLAGS)
mpiflag(flag::Symbol) = getfirst(x -> x.symbol==flag, MPIFLAGS)

function parse_mpiflags(line::Vector{<:AbstractString})
    eflags = ExecFlag[]
    i = 1
    while i <= length(line)
        s = line[i]
        if s[1] != '-'
            break
        end
        mflag = if s[2] == '-'
            mpiflag(strip(s, '-'))
        else
            mpiflag(Symbol(strip(s, '-')))
        end

        @assert mflag != nothing "$(strip(s, '-')) is not a recognized mpiflag"
        flagtype = mflag.type
        val = if flagtype ==  String
                  val = line[i+1]
                  i+=2
                  val
              elseif flagtype == Vector{String}
                  tval = []
                  while i+1 <= length(line) && !occursin('-', line[i+1])
                       push!(tval, line[i+1])
                       i += 1
                   end
                   i+=1
                   tval
               elseif flagtype == Int
                   if line[i+1][1] == '\$'
                       line[i+1]
                   else
                       parse(Int, line[i+1])
                   end
                   i += 2
               elseif flagtype == Vector{Int}
                   tval = []
                   while i+1 <= length(line) && !occursin('-', line[i+1])
                       if line[i+1][1] == '\$'
                           push!(tval, line[i+1])
                       else
                           push!(tval, parse(Int,line[i+1]))
                       end
                       i += 1
                   end
                   i+=1
                   tval
               elseif flagtype <: Pair
                   tval = []
                   while i+1<=length(line) && !occursin('-', line[i+1])
                       if line[i+1][1] == '\$'
                           push!(tval, line[i+1])
                       else
                           for p in flagtype.parameters
                               if p <: AbstractString
                                   push!(tval, line[i+1])
                                   break
                               end
                               tparse = tryparse(p, line[i+1])
                               if tparse != nothing
                                   push!(tval, tparse)
                                   break
                               end
                           end
                       end
                       i += 1
                   end
                   i+=1
                   tval
               end
        push!(eflags, ExecFlag(mflag, val))
    end
    eflags
end

const RUNEXECS = ["mpirun"]
allexecs() = vcat(RUNEXECS, QEEXECS, WANEXECS)
parseable_execs() = vcat(QEEXECS, WANEXECS)

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

function inputparser(exec::Exec)
    if exec.exec ∈ QEEXECS
        read_qe_input
    elseif exec.exec ∈ WANEXECS
        read_wannier_input
    end
end

function setflags!(exec::Exec, flags...)
    for (f, val) in flags
        flag = isa(f, String) ? getfirst(x -> x.name == f, exec.flags) : getfirst(x -> x.symbol == f, exec.flags)
        if flag != nothing
            flag.value = convert(flag.type, val)
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

setexecdir!(exec::Exec, dir) = exec.dir = dir
