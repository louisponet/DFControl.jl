struct Orbital
    name  ::Symbol
    size::Int
    l   ::Int
end
const orbitals = [
    Orbital(:s, 1, 0),
    Orbital(:p, 3, 1),
    Orbital(:d, 5, 2),
    Orbital(:f, 7, 3),
    Orbital(:dz2, 1, 2),
    Orbital(:dxz, 1, 2),
    Orbital(:dyz, 1, 2),
    Orbital(Symbol("dx2-y2"), 1, 2),
    Orbital(:dxy, 1, 2),
]
orbital(s::Symbol) = getfirst(x -> x.name == s, orbitals)

Base.convert(::Type{Symbol}, x::Orbital) = x.name
orbsize(orb::Orbital) = orb.size
orbsize(orb::Symbol)  = orbsize(orbital(orb))

@with_kw struct Projection
    orb   ::Orbital = orbitals[1]
    start ::Int = 0
    last  ::Int = 0
end

orbital(proj::Projection) = proj.orb
orbsize(proj::Projection) = proj.last - proj.start + 1

"""
Adds projections to atoms.
"""
function addprojections!(atoms, projections_)
    t_start = 1
    for (proj_at, projs) in projections_
        for at in atoms
            if length(projs) <= length(projections(at))
                continue
            end
            for proj in projs
                if name(at) == proj_at
                    size   = orbsize(proj)
                    t_proj = Projection(orbital(proj), t_start, t_start + size - 1)
                    push!(projections(at), t_proj)
                    t_start += size
                end
            end
        end
    end
end

Base.zero(::Type{Projection}) = Projection(Orbital(:zero, 0, 0), 0, 0)
Base.zero(::Projection) = Projection(Orbital(:zero, 0, 0), 0, 0)

Base.range(proj::Projection)  = proj.start:proj.last
