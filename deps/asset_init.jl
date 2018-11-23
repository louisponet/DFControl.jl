function writefbodyline(f, indent, s)
    for i=1:indent
        write(f, "\t")
    end
    write(f, "$s\n")
end

function fort2julia(f_type)
    f_type = lowercase(f_type)
    if f_type == "real"
        return Float32
    elseif f_type == "real(kind=dp)"
        return Float64
    elseif f_type == "complex(kind=dp)"
        return Complex{Float64}
    elseif occursin("character", f_type)
        return String
    elseif f_type == "string"
        return String
    elseif f_type == "integer"
        return Int
    elseif f_type == "logical"
        return Bool
    elseif occursin(".D", f_type)
        return replace(f_type, "D" => "e")
    else
        return Nothing
    end
end

open(joinpath(@__DIR__, "mpirunflags.jl"), "w") do wf
    write(wf, "_MPIFLAGS() = ExecFlag[\n")
    # writefbodyline(1, "flags = ")
    open(joinpath(@__DIR__, "..", "assets","mpirun_man.txt"), "r") do f
        line = readline(f)
        while line != "OPTIONS"
            line = readline(f)
        end
        while line != "Environment Variables"
            line = strip(readline(f))
            if !isempty(line) && line[1] == '-'
                name        = ""     #--
                symbols     = Symbol[] #-
                type        = Nothing
                description = ""
                sline = strip.(split(line), ',')
                for s in sline
                    if s[2] == '-' #--npernode
                        name = strip(s, '-')
                    elseif occursin('<', s) # <#persocket>
                        type = if occursin("#", s)
                                   occursin(',', s) ? Vector{Int} : Int
                               elseif occursin("ppr", s) #ppr:N:<object>
                                   Pair{Int, String}
                               else
                                   occursin(',', s) ? Vector{String} : String
                               end
                        break
                    else  #-np
                        push!(symbols, Symbol(strip(s, '-')))
                    end
                end
                line = strip(readline(f))
                while !isempty(line)
                    description *= " " * line
                    line = strip(readline(f))
                end
                description = replace(strip(description), "\"" => "'")
                if name != "" && isempty(symbols)
                    symbols = [Symbol(name)]
                end
                for symbol in symbols
                    writefbodyline(wf, 1,"""ExecFlag(Symbol("$symbol"), "$name", $type, "$description", nothing),""")
                end
            end
        end
    end
    writefbodyline(wf, 0,"]")
end

open(joinpath(@__DIR__, "wannier90flags.jl"), "w") do wf
    write(wf, "_WANFLAGS() = Dict{Symbol, Type}(\n")
    open(joinpath(@__DIR__, "..", "assets", "inputs", "wannier", "input_flags.txt"), "r") do f
        while !eof(f)
            line = readline(f)
            if isempty(line) || line[1] == '!'
                continue
            else
                s_line    = split(line, "::")
                flag      = Symbol(strip(strip(split(split(split(s_line[end], "=")[1],"(")[1],"!")[1], ':')))
                fl_type   = fort2julia(strip(split(s_line[1])[1],','))
                writefbodyline(wf, 1, """Symbol("$flag") => $fl_type,""")
            end
        end
    end
    write(wf, ")")
end
