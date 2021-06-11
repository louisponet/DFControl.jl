struct Orbital
    name::String
    size::Int
    l   ::Int
    mr  ::Int
end
const ORBITALS = [
    Orbital("s", 1, 0, 1),
    Orbital("p", 3, 1, 0),
    Orbital("px", 1, 1, 2),
    Orbital("py", 1, 1, 3),
    Orbital("pz", 1, 1, 1),
    Orbital("d", 5, 2, 0),
    Orbital("dz2", 1, 2, 1), 
    Orbital("dxz", 1, 2, 2), 
    Orbital("dyz", 1, 2, 3), 
    Orbital("dx2-y2", 1, 2, 4),
    Orbital("dxy", 1, 2, 5), 
    Orbital("f", 7, 3, 0),
    Orbital("fz3", 1, 3, 1),
    Orbital("fxz2", 1, 3, 2),
    Orbital("fyz2", 1, 3, 3),
    Orbital("fz(x2-y2)", 1, 3, 4),
    Orbital("fxyz", 1, 3, 5),
    Orbital("fx(x2-3y2)", 1, 3, 6),
    Orbital("fy(3x2-y2)", 1, 3, 7),
    Orbital("sp", 2, -1, 0),
    Orbital("sp2", 3, -2, 0),
    Orbital("sp3", 4, -3, 0),
    Orbital("sp3d2", 6, -5, 0),
]
orbital(s::String) = getfirst(x -> x.name == s, ORBITALS)
orbital(l::Number) = getfirst(x -> x.l == l, ORBITALS)

Base.convert(::Type{String}, x::Orbital) = x.name
orbsize(orb::Orbital) = orb.size
orbsize(orb::String)  = orbsize(orbital(orb))

@with_kw_noshow struct Projection
    orb   ::Orbital = ORBITALS[1]
    start ::Int = 0
    last  ::Int = 0
end

orbital(proj::Projection) = proj.orb
orbsize(proj::Projection) = proj.last - proj.start + 1

function Base.show(io::IO, proj::Projection)
    dfprintln(io, crayon"cyan", "Orbital: ", crayon"reset", "$(proj.orb.name)") 
    dfprintln(io, crayon"red", "start index: ", crayon"reset", "$(proj.start)")
    dfprintln(io, crayon"red", "last index: ", crayon"reset", "$(proj.last)")
end

"""
Adds projections to atoms.
"""
function addprojections!(atoms, projections_, soc; print=true)
    t_start = 1
    for (proj_at, projs) in projections_
        for proj in projs
            for at in atoms
                if length(projs) <= length(projections(at))
                    continue
                end
                if name(at) == proj_at
                    size   = soc ? 2*orbsize(proj) : orbsize(proj)
                    t_proj = Projection(orbital(proj), t_start, t_start + size - 1)
                    push!(projections(at), t_proj)
                    t_start += size
                end
            end
        end
    end
    print && @info "Projections are now:"
    for at in atoms
        if !isempty(projections(at))
            print && @info "$(name(at));$(position_cryst(at)) :"
            for p in projections(at)
                print && Base.show(p)
            end
        end
    end
end

Base.zero(::Type{Projection}) = Projection(Orbital("zero", 0, 0, 0), 0, 0)
Base.zero(::Projection) = Projection(Orbital("zero", 0, 0, 0), 0, 0)


Base.range(proj::Projection)  = proj.start:proj.last
