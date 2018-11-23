
open(joinpath(@__DIR__, "mpirunflags.jl"), "w") do wf
    function writefbodyline(indent, s)
        for i=1:indent
            write(wf, "\t")
        end
        write(wf, "$s\n")
    end
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
                    writefbodyline(1,"""ExecFlag(Symbol("$symbol"), "$name", $type, "$description", nothing),""")
                end
            end
        end
    end
    writefbodyline(1,"""ExecFlag(:dummy, "dummy", Nothing, "dummy", nothing)""")
    writefbodyline(0,"]")
end
