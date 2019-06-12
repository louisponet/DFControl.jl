using StaticArrays
const Vec = SVector

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
    write(wf, "_MPI_FLAGS() = ExecFlag[\n")
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
    writefbodyline(wf, 0, "]")
end

open(joinpath(@__DIR__, "wannier90flags.jl"), "w") do wf
    write(wf, "_WAN_FLAGS() = Dict{Symbol, Type}(\n")
    open(joinpath(@__DIR__, "..", "assets", "inputs", "wannier", "input_flags.txt"), "r") do f
        while !eof(f)
            line = readline(f)
            if isempty(line) || line[1] == '!'
                continue
            else
                s_line    = split(line, "::")
                flagpart = s_line[end]
                flag      = Symbol(strip(strip(split(split(split(s_line[end], "=")[1],"(")[1],"!")[1], ':')))
                fl_type   = fort2julia(strip(split(s_line[1])[1],','))
                flagstring = occursin("(", flagpart) && occursin(")", flagpart) ? """Symbol("$flag") => Vector{$fl_type},""" : """Symbol("$flag") => $fl_type,"""
                writefbodyline(wf, 1, flagstring)
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
            description *= strip_split(line,":")[2] * "\n"
            i += 1
            line = lines[i]
            while !occursin("+------", line)
                description *= strip(lines[i]) * "\n"
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

open(joinpath(@__DIR__, "abinitflags.jl"), "w") do wf
    flaglines = split.(readlines(joinpath(@__DIR__,  "..", "assets", "inputs", "abinit", "input_variables.txt")))
    write(wf, "_ABINITFLAGS() = Dict{Symbol, Type}(\n")
    for (fl, typ) in flaglines
        write(wf, "\t:$fl => $typ,\n")
    end
    write(wf, ")")
end


#all this so I can be lazy and just use repr
struct ElkFlagInfo{T}
    name::Symbol
    default::Union{T, Nothing}  #check again v0.7 Some
    description::String
end
ElkFlagInfo() = ElkFlagInfo{Nothing}(:error, "")
Base.eltype(x::ElkFlagInfo{T}) where T = T


struct ElkControlBlockInfo
    name::Symbol
    flags::Vector{<:ElkFlagInfo}
    description::String
end

function closing(c::Char)
	if c == '{'
		return '}'
	elseif c == '('
		return ')'
	elseif c == '['
		return ']'
	end
end

function readuntil_closing(io::IO, opening_char::Char)
	out = ""
	counter = 1
	closing_char = closing(opening_char)
	while counter > 0
		c = read(io, Char)
		out *= c
		if c == opening_char
			counter += 1
		elseif c == closing_char
			counter -= 1
		end
	end
	return out
end

function elk2julia_type(s::AbstractString)
	if s == "integer"
		return Int
	elseif s == "real"
		return Float64
	elseif s == "logical"
		return Bool
	elseif s == "string"
		return String
	elseif s == "complex"
		return ComplexF64
	end
end

function replace_multiple(str, replacements::Pair{String, String}...)
    tstr = deepcopy(str)
    for r in replacements
        tstr = replace(tstr, r)
    end
    return tstr
end

function parse_elk_default(::Type{T}, s) where T
	s_cleaned =strip(strip(strip(replace_multiple(s, "{" => "", "}" => "", "."=>"", " "=>""), '\$'), ')'), '(')
	if occursin("-", s_cleaned)
		return nothing
	elseif occursin(",", s_cleaned)
		s_t = split(s_cleaned, ',')
		return parse.(eltype(T), s_t)
	else
		return parse(T, s_cleaned)
	end
end

function blockflags(s::String)
	flag_regex = r"\{([^.}]*)\}"
	flags = ElkFlagInfo[]
	lines = split(s, "\\\\")
	for l in lines
		sline = strip_split(l, '&')
		flag = match(flag_regex, sline[1])
        descr = sline[2]
        flag_i = split(replace_multiple(sline[3], "(" => " ", ")" => " "))
        typ = length(flag_i) > 1 && flag_i[1] != "string" ? Vec{parse(Int, flag_i[2]), elk2julia_type(flag_i[1])} : elk2julia_type(flag_i[1])
		#TODO default
        default = typ == String ? sline[4] : parse_elk_default(typ, sline[4])
        flagname = Symbol(strip_split(flag.captures[1], '(')[1])
        if flagname == :wann_bands
	        push!(flags, ElkFlagInfo{UnitRange{Int}}(:wann_bands, nothing, string(descr)))
		elseif flagname == :wann_projections
	        push!(flags, ElkFlagInfo{Vector{String}}(:wann_projections, nothing, string(descr)))
		elseif flagname == :wann_seedname
	        push!(flags, ElkFlagInfo{Symbol}(:wann_seedname, nothing, string(descr)))
		else
	        push!(flags, ElkFlagInfo{typ}(flagname, default, string(descr)))
        end
	end
	return flags
end

replace_latex_symbols(s::String) = replace_multiple(s, "\\_" => "_", "\\hline" => "", "\\tt "=> "", "^"=>"e")

function read_elk_doc(fn::String)
	blocks = ElkControlBlockInfo[]
	blockname_regex = r"\{(.*)\}"
	open(fn, "r") do f
		while !eof(f)
			line = readline(f)
			if occursin("\\block", line)
				m = match(blockname_regex, line)
				blockname = Symbol(m.captures[1])
				block_text = replace_latex_symbols(readuntil_closing(f, '{'))
				flags = blockflags(block_text)
				push!(blocks, ElkControlBlockInfo(blockname, flags, readuntil(f, "\\block")))
				seek(f, Base.position(f) - 6) #\\block offset
			end
		end
	end
	return blocks
end

ELK_CONTROLBLOCKS = read_elk_doc(joinpath(@__DIR__, "..","assets", "inputs", "elk", "elk.txt"))

open(joinpath(@__DIR__,"elkflags.jl"), "w") do f
	write(f, "_ELK_CONTROLBLOCKS() = $(repr(ELK_CONTROLBLOCKS))")
end
