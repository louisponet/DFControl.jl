
@enum Orbital s p d f
function Orbital(s::Symbol)
    t = 0
    while Symbol(Orbital(t)) != s
        t += 1
        if t > Int(f)
            error("Orbital $s not defined.")
        end
    end
    return t
end
Base.convert(::Type{Symbol}, x::Orbital) = Symbol(x)
orbsize(orbital::Orbital) = Int(orbital) * 2 + 1
orbsize(orbital::Symbol)  = Orbital(orbital) * 2 + 1

@with_kw struct Projection
    orb   ::Orbital = s
    start ::Int = 0
    size  ::Int = 0
    last  ::Int = 0
end

"""
    add_projections(projections, atoms)

Takes an array of `Pair{Symbol, Orbital}` where Symbol signifies the atom symbol for the projections, and an Array of `Atom` and then assigns the correct `Projection` arrays to each atom.
"""
function add_projections(projections, atoms)
    t_start = 1
    for (proj_at, projs) in projections
        for proj in projs
            for at in atoms
                if id(at) == proj_at
                    size = orbsize(proj)
                    t_proj = Projection(Orbital(proj), t_start, size, t_start + size - 1)
                    if !isdefined(at, :projections)
                        setprojections!(at, [t_proj])
                    else
                        push!(projections(at), t_proj)
                    end
                    t_start += size
                end
            end
        end
    end
end

Base.zero(::Type{Projection}) = Projection(s, 0, 0, 0)
Base.range(proj::Projection)  = proj.start:proj.last
