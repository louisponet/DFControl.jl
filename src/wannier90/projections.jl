
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
orbsize(orbital::Orbital) = Int(orbital) * 2 + 1
orbsize(orbital::Symbol)  = Orbital(orbital) * 2 + 1 

struct WannProjection
    orb::Symbol
    start::Int
    size::Int
    last::Int
end

function add_projections(projections, atoms)
    t_start = 1
    for (proj_at, projs) in projections
        for at in atoms
            if at.id != proj_at
                continue
            end
            t_projs = WannProjection[]
            for proj in projs 
                size = orbsize(proj)
                push!(t_projs, WannProjection(proj, t_start, size, t_start + size - 1))
                t_start += size
            end
            at.data[:projections] = t_projs
        end
    end
end