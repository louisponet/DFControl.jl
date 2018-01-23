
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