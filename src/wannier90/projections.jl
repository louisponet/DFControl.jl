
@enum Orbital s p d f
function orbital(s::Symbol)
    t = 0
    while Symbol(Orbital(t)) != s
        t += 1
        if t > Int(f)
            error("Orbital $s not defined.")
        end
    end
    return Orbital(t)
end
Base.convert(::Type{Symbol}, x::Orbital) = Symbol(x)
orbsize(orb::Orbital) = Int(orb) * 2 + 1
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
