
struct WannProjection{T <: AbstractFloat}
    orb::Symbol
    atom::Symbol
    start::Int
    size::Int
    last::Int
    position::Point3D{T}
end


function get_wan_projections(filename::String, T=Float64)
    projections = OrderedDict()
    atoms       = Tuple[]
    n_proj      = 0
    open(filename, "r") do f
        while !eof(f)
            line = lowercase(readline(f))
            if contains(line, "begin projections")
                line = readline(f)
                while !contains(lowercase(line), "end")
                    if contains(line, "!") || isempty(line)
                        line = readline(f)
                        continue
                    end
                    if contains(line, "random")
                        error("Can't read the atomic info when projections are random!")
                    end
                    split_line   = DFControl.strip_split(line, ':')
                    atom         = Symbol(split_line[1])
                    _projections = [Symbol(proj) for proj in DFControl.strip_split(split_line[2], ';')]
                    projections[atom] = _projections
                    line = readline(f)
                end
            elseif contains(line, "begin") && contains(line, "atoms")
                line = readline(f)
                while !contains(lowercase(line), "end")
                    if contains(line, "!")
                    line = readline(f)
                        continue
                    end
                    split_line = DFControl.strip_split(line)
                    atom       = Symbol(split_line[1])
                    position   = Point3D(DFControl.parse_string_array(T, split_line[2:4]))
                    push!(atoms,(atom, position))
                    n_proj += length(projections[atom])
                    line = readline(f)
                end
            end
        end
    end

    t_out = Array{WannProjection, 1}(n_proj)
    t_start = 1
    i = 1 
    for (proj_at, projs) in projections
        for proj in projs
            for (pos_at, pos) in atoms
                if proj_at != pos_at
                    continue
                end
                size = orbsize(proj)
                t_out[i] = WannProjection(proj, pos_at, t_start, size, t_start + size - 1, pos)
                i += 1
                t_start += size
            end
        end
    end
    out = WannProjection[] 
    i = 1
    for (pos_at, pos) in atoms
        for proj in t_out
            if proj.position == pos
                push!(out, proj)
            end
        end
    end
    return out 
end

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

