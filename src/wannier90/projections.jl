struct Orbital
    id::Symbol
    size::Int
end
const orbitals = [
    Orbital(:s, 1),
    Orbital(:p, 3),
    Orbital(:d, 5),
    Orbital(:f, 7),
    Orbital(:dz2, 1),
    Orbital(:dxz, 1),
    Orbital(:dyz, 1),
    Orbital(Symbol("dx2-y2"), 1),
    Orbital(:dxy,1),
]
orbital(s::Symbol) = getfirst(x->x.id==s, orbitals)

Base.convert(::Type{Symbol}, x::Orbital) = x.id
orbsize(orb::Orbital) = orb.size
orbsize(orb::Symbol)  = orbsize(orbital(orb))

@with_kw struct Projection
    orb   ::Orbital = s
    start ::Int = 0
    size  ::Int = 0
    last  ::Int = 0
end

orbsize(proj::Projection) = proj.size

"""
Adds projections to atoms.
"""
function addprojections!(projections_, atoms)
    t_start = 1
    for (proj_at, projs) in projections_
        for at in atoms
            if length(projs) <= length(projections(at))
                continue
            end
            for proj in projs
                if id(at) == proj_at
                    size = orbsize(proj)
                    t_proj = Projection(orbital(proj), t_start, size, t_start + size - 1)
                    push!(projections(at), t_proj)
                    t_start += size
                end
            end
        end
    end
end

Base.zero(::Type{Projection}) = Projection(s, 0, 0, 0)
Base.range(proj::Projection)  = proj.start:proj.last
