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
