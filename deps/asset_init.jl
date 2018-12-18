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
indentabs(indent) = prod(["\t" for i=1:indent])

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
                    description *= join(split(line), " ") * "\n" * indentabs(2)
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
    writefbodyline(wf, 1, ":preprocess => Bool") #temporary hack
    write(wf, ")")
end

strip_split(line, args...) = strip.(split(line, args...))
function read_block(f, startstr::String, endstr::String)
    block = [startstr]
    line = readline(f)
    while !eof(f)
        while !occursin(endstr, line)
            line = readline(f)
            push!(block, line)
        end
        return block
    end
    error("Block not found: start = $startstr, end = $endstr.")
end

function vardim(line)
    if all(occursin.(['(', ')'], (line,)))
        dim = length(split(split(split(line, '(')[2], ')')[1], ','))
    else
        dim = 0
    end
    return dim
end

function qe_write_variable(wf, indent, lines, i)
    name = gensym()
    var_i = i
    i += 2
    typ = fort2julia(strip_split(lines[i])[2])
    description = ""
    # default = nothing
    i += 1
    line = lines[i]
    while !occursin("+------", line)
        # if occursin("Default", line)
        #     _t = strip_split(line)[2]
        #     _t = strip(strip(_t,'('),')')
        #     if occursin("D", _t)
        #         default = Meta.parse(typ, replace(_t,"D" => "e"))
        #     else
        #         _t = occursin("=",_t) ?split(_t,"=")[end] : _t
        #         default = typ ==String ? _t : Meta.parse(_t)
        #         println(_t)
        #         if typeof(default) != Symbol
        #             default = convert(typ, default)
        #         end
        #     end
        if occursin("Description", line)
            description *= strip_split(line,":")[2] * "\n" * indentabs(indent+1)
            i += 1
            line = lines[i]
            while !occursin("+------", line)
                description *= strip(lines[i]) * "\n" * indentabs(indent+1)
                i += 1
                line = lines[i]
            end
            @goto break_label
        end
        i += 1
        line = lines[i]
    end
    @label break_label
    line = lines[var_i]
    dim = vardim(line)
    typ = if dim == 2
            Matrix{typ}
        elseif dim == 1
            Vector{typ}
        else
            typ
        end
    description = replace(replace(replace(replace(replace(description, "\"" => "'"), "\\" => "\\\\"), "\$" => "\\\$"), "\\\\t" => "\\t"), "\\\\n" => "\\n")
    if occursin("Variables", line)
        spl = [split(x,"(")[1] for x in strip.(filter(x -> !occursin("=", x), split(line)[2:end]), ',')]
        names = Symbol.(spl)
        for name in names
            writefbodyline(wf, indent,  """QEFlagInfo{$typ}(Symbol("$name"), "$description"),""")
        end
        return i
    else
        if occursin("(", line) && occursin(")", line)
            name = Symbol(split(strip_split(line, ":")[end],"(")[1])
        else
            name = Symbol(strip_split(line,":")[end])
        end
        writefbodyline(wf, indent, """QEFlagInfo{$typ}(Symbol("$name"), "$description"),""")
        return i
    end
end

function write_QEControlBlockInfo(wf, indent, lines)
    name  = Symbol(lowercase(strip_split(lines[1], "&")[2]))
    writefbodyline(wf, indent, """QEControlBlockInfo(Symbol("$name"), [""")
    for i=1:length(lines)
        line = lines[i]
        if occursin("Variable", line)
            i += qe_write_variable(wf, indent+1, lines, i)
        end
    end
    writefbodyline(wf, indent + 1, "]),")
    # writefbodyline(wf, indent, "),")
end

function write_QEDataBlockInfo(wf, indent, lines)
    spl                 = split(strip(lines[1]))
    name                = length(spl) > 1 ? lowercase(spl[2]) : "noname"
    options             = Symbol.(spl[4:2:end])
    description         = ""
    options_description = ""
    istop = findfirst(x -> occursin("DESCRIPTION", x), lines)
    for i=1:istop-1
        line = strip(replace(replace(lines[i], "_" => ""), "-" => ""))
        if isempty(line)
            continue
        end
        description *= strip(line) * "\n" * indentabs(indent+2)
    end
    description = replace(description, "\"" => "'")


    istart = findnext(x -> occursin("Card's flags", x), lines, istop)
    if istart != nothing
        istop = findnext(x -> occursin("+------", x), lines, istart)
        for i=istart+1:istop - 1
            line = lines[i]
            if occursin("Default", line) || occursin("Description", line)
                continue
            end
            options_description *= strip(line) * "\n" * indentabs(indent+2)
        end
    end

    options_description = replace(options_description, "\"" => "'")
    writefbodyline(wf, indent+1, """QEDataBlockInfo(Symbol("$name"),"$description", $options, "$options_description", [""")
    # writefbodyline(wf, indent, """QEDataBlockInfo(Symbol("$name"),""")
    i = istop
    while i <= length(lines) - 1
        line = strip(lines[i])
        if occursin("Variable", line)
            i = qe_write_variable(wf, indent+2, lines, i)
        end
        i += 1
    end
    writefbodyline(wf, indent+1, "]),")
end
function write_QEInputInfo(wf, filename, indent, exec)
    writefbodyline(wf, indent, """QEInputInfo("$exec", [""")
    allcontrolwritten = false
    anydatawritten = false
    open(filename, "r") do f
        while !eof(f)
            line = readline(f)
            if occursin("NAMELIST", line)
                write_QEControlBlockInfo(wf, indent+1, read_block(f, line, "END OF NAMELIST"))
            elseif occursin("CARD:", line)
                anydatawritten = true
                if !allcontrolwritten
                    writefbodyline(wf, indent+1, "],")
                    writefbodyline(wf, indent+1, "[")
                    allcontrolwritten = true
                end
                write_QEDataBlockInfo(wf, indent+1, read_block(f, line, "END OF CARD"))
            end
        end
    end
    if !anydatawritten
        if !allcontrolwritten
            writefbodyline(wf, indent+1, "],")
        end
        writefbodyline(wf, indent+1, "QEDataBlockInfo[]")
    else
        writefbodyline(wf, indent+1, "]")
    end
    writefbodyline(wf, indent, "),")
end

searchdir(path::String, key) = filter(x -> occursin(key, x), readdir(path))
open(joinpath(@__DIR__, "qeflags.jl"), "w") do wf
    write(wf, "_QEINPUTINFOS() = QEInputInfo[\n")
    input_files = searchdir(joinpath(@__DIR__,  "..", "assets", "inputs", "qe"), "INPUT")
    filepaths  = joinpath.(Ref(joinpath(@__DIR__, "..", "assets", "inputs", "qe")), input_files)
    for _f in filepaths
        exec_name = join([lowercase(splitext(split(_f, "_")[end])[1]),".x"],"")
        write_QEInputInfo(wf, _f, 1, exec_name)
    end
    write(wf, "]")
end
