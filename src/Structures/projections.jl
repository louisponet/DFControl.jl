"""
    Orbital(name::String, size::Int, l::Int, mr::Int)

Wannier90 orbital as defined in the [Wannier90 User Guide](https://github.com/wannier-developers/wannier90/raw/v3.1.0/doc/compiled_docs/user_guide.pdf). The `name` will be written in the `projections` block of the Wannier90 input.
"""
struct Orbital
    name :: String
    size :: Int
    l    :: Int
    mr   :: Int
end
Orbital() = Orbital("none", 0, 0, 0)
Orbital(dict::JSON3.Object) = Orbital(dict[:name], dict[:size], dict[:l], dict[:mr])
StructTypes.StructType(::Type{Orbital}) = StructTypes.Struct()

const ORBITALS = [Orbital("s", 1, 0, 1), Orbital("p", 3, 1, 0), Orbital("px", 1, 1, 2),
                  Orbital("py", 1, 1, 3), Orbital("pz", 1, 1, 1), Orbital("d", 5, 2, 0),
                  Orbital("dz2", 1, 2, 1), Orbital("dxz", 1, 2, 2), Orbital("dyz", 1, 2, 3),
                  Orbital("dx2-y2", 1, 2, 4), Orbital("dxy", 1, 2, 5),
                  Orbital("f", 7, 3, 0), Orbital("fz3", 1, 3, 1), Orbital("fxz2", 1, 3, 2),
                  Orbital("fyz2", 1, 3, 3), Orbital("fz(x2-y2)", 1, 3, 4),
                  Orbital("fxyz", 1, 3, 5), Orbital("fx(x2-3y2)", 1, 3, 6),
                  Orbital("fy(3x2-y2)", 1, 3, 7), Orbital("sp", 2, -1, 0),
                  Orbital("sp2", 3, -2, 0), Orbital("sp3", 4, -3, 0),
                  Orbital("sp3d2", 6, -5, 0)]

"""
    orbital(s::String)

Returns the [Orbital](@ref Orbitals) identified by `s`. 
"""
orbital(s::String) = getfirst(x -> x.name == s, ORBITALS)
orbital(l::Number) = getfirst(x -> x.l == l, ORBITALS)

Base.convert(::Type{String}, x::Orbital) = x.name
Base.length(o::Orbital) = o.size

Base.hash(orbital::Orbital, h::UInt) = hash(orbital.mr, hash(orbital.l, h))
Base.:(==)(o1::Orbital, o2::Orbital) = o1.l == o2.l && o1.mr == o2.mr

"""
    Projection(orbital::Orbital, start::Int, last::Int)

A Wannier90 `Projection`, representing an `Orbital` with indices from `start` to `last`.
"""
mutable struct Projection
    orbital::Orbital
    start::Int
    last::Int
end
Projection() = Projection(Orbital(), 0, 0)
Projection(o::Orbital) = Projection(o, 0, 0)
Projection(o::String) = Projection(orbital(o), 0, 0) 
Projection(o::Orbital, start::Int) = Projection(o, start, start + o.size - 1)
Projection(ostr::String, start::Int) = (o = orbital(ostr); Projection(o, start))
StructTypes.StructType(::Type{Projection}) = StructTypes.Struct()
Projection(dict::JSON3.Object) = Projection(Orbital(dict[:orbital]), dict[:start], dict[:last])

function sanitize!(projs::Vector{Projection}, soc::Bool)
    id = 1
    for proj in projs
        size = soc ? 2 * proj.orbital.size : proj.orbital.size
        proj.start = id
        proj.last = id + size - 1
        id = proj.last + 1
    end
    return projs
end

Base.range(proj::Projection) = proj.start:proj.last
Base.length(proj::Projection) = length(proj.orbital)

function Base.hash(data::Projection, h::UInt)
    for f in fieldnames(Projection)
        h = hash(getfield(data, f), h)
    end
    return h
end
