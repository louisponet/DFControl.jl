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

mutable struct Exec
    exec ::String
    dir  ::String
    flags::Dict{Symbol, Any}
    function Exec(exec::String, dir::String, flags::Dict)
        return new(exec, dir, convert(SymAnyDict, flags))
    end
end

Exec() = Exec("")
Exec(exec::String) = Exec(exec, "")
Exec(exec::String, dir::String) = Exec(exec, dir, SymAnyDict())
Exec(exec::String, dir::String, flags...) = Exec(exec, dir, SymAnyDict(flags))

function setflags!(exec::Exec, flags...)
    for (f, val) in flags
        exec.flags[f] = val
    end
    exec.flags
end
function rmflags!(exec::Exec, flags...)
    for f in flags
        pop!(exec.flags, f)
    end
    exec.flags
end

setexecdir!(exec::Exec, dir) = exec.dir = dir

struct ExecFlag
    symbol     ::Symbol
    name       ::String
    type       ::Type
    description::String
end

const mpiflags = ExecFlag[]

function parse_execflags(file)
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
                symbol      = :nosym #-
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
                        symbol = Symbol(strip(s, '-'))
                    end
                end
                line = strip(readline(f))
                while !isempty(line)
                    description *= " " * line
                    line = strip(readline(f))
                end
                push!(flags, ExecFlag(symbol, name, type, strip(description)))
            end
        end
    end
    return flags
end

init_mpiflags() = append!(mpiflags, parse_execflags(joinpath(@__DIR__, "..", "assets","mpirun_man.txt")))
