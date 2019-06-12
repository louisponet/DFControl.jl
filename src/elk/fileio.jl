function read_atoms(fn::IO)
	sane_readline() = strip(strip_split(readline(fn), ":")[1], '\'')
	nspecies = parse(Int, sane_readline())
	species_counter = 0

	reading_atoms = false
	atname = nothing
	natoms = 0
	atcounter = 0
	positions = Point3{Float64}[] 
	bfcmts    = Vec3{Float64}[]
	atoms     = Atom{Float64}[]
	atname = Symbol(split(sane_readline(), ".")[1])

	while species_counter < nspecies
		l = sane_readline()
		isempty(l) && continue
		if !reading_atoms
			natoms = parse(Int, l)
			reading_atoms = true
		elseif atcounter < natoms
			sline = split(l)
			push!(positions, Point3{Float64}(parse.(Float64, sline[1:3])))
			push!(bfcmts, Vec3{Float64}(parse.(Float64, sline[4:end])))
			atcounter += 1
		else
			for (p, mag) in zip(positions, bfcmts)
				push!(atoms, Atom(name=atname, element=element(atname), position=p, magnetization=mag))
			end
			atname = Symbol(split(l, ".")[1])
			reading_atoms   = false
			natoms          = 0
			atcounter       = 0
			positions       = Point3{Float64}[] 
			bfcmts          = Vec3{Float64}[]
			species_counter += 1
		end
	end
	return atoms
end


"""
    elk_read_input(filename; execs=[Exec("elk")], run=true, structure_name="noname")

Reads an Elk input file. The `ELK_EXEC` inside execs gets used to find which flags are allowed in this input file, and convert the read values to the correct Types.
Returns a `DFInput{Elk}` and the `Structure` that is found in the input.
"""
function elk_read_input(fn::String; execs=[Exec("elk")], run=true, structure_name="noname")
	blocknames_flaglines = Dict{Symbol, Any}()
	atoms = Atom{Float64}[]
	inputs = DFInput[]
	dir, file = splitdir(fn)
	open(fn, "r") do f
		sane_readline() = strip_split(readline(f), ':')[1]
		while !eof(f)
			line = sane_readline()
			isempty(line) || line[1] == '!' && continue
			if !isempty(line)
				blockname = Symbol(line)
				info = elk_block_info(blockname)
				if blockname == :atoms
					atoms = read_atoms(f)
				elseif blockname == :wannierExtra
					@warn "Please supply the .win file instead of wannierExtra block for extracting the Wannier DFInput."

					wflags, wdata, ab, cb, proj_block = wan_read_input(f)
					continue
					# push!(inputs, DFInput{Wannier90}("wannier", dir, wflags, wdata, execs, run))
				elseif info != nothing
					blocklines = String[]
					line = sane_readline()
					while !isempty(line)
						push!(blocklines, strip(replace(strip(split(line, ':')[1]), "'" => ""), '.'))
						line = sane_readline()
					end
					if blockname != :tasks && blockname != :avec
						blocknames_flaglines[blockname] = parse_block_flags(blockname, collect(Iterators.flatten(split.(filter(x->!isempty(x), blocklines)))))
					else
						blocknames_flaglines[blockname] = filter(x->!isempty(x), blocklines)
					end
				end
			end
		end
	end
	# for i in 
	#cell
	scale = haskey(blocknames_flaglines, :scale) ? last(pop!(blocknames_flaglines[:scale])) : 1.0
	cell  = Mat3(reshape(scale .* parse.(Float64, collect(Iterators.flatten(split.(pop!(blocknames_flaglines, :avec))))), 3, 3))
	#structure #TODO make positions be in lattices coordinates
	newats = [Atom(at, cell' * position(at)) for at in atoms]
	structure = Structure(structure_name, cell, newats)

	#different tasks
	tasks_ = pop!(blocknames_flaglines, :tasks)
	should_run = map(x -> !occursin("!", x), tasks_)
	tasks      = map(x -> strip(strip(x, '!')), tasks_)

	#The wan tasks can prbably be skipped completely and added depending on the wannier input
	wan_tasks = String[]
	for (x, run) in zip(tasks, should_run)
		if x ∈ ["601", "602", "603", "604", "605"]
			push!(wan_tasks, x)
		else
			push!(inputs, input_from_task(x, blocknames_flaglines, dir, execs, run))
		end
	end
	if haskey(blocknames_flaglines, :wannier) && !isempty(wan_tasks)
		flags = pop!(blocknames_flaglines, :wannier)
		# flags[:elk2wan_tasks] = wan_tasks
		push!(inputs, DFInput{Elk}(name="elk2wannier", dir=dir, flags=flags, execs=execs, run=true))
	end
	for f in (:ngrid, :vkloff, :plot1d, :plot2d, :plot3d)
		haskey(blocknames_flaglines, f) && pop!(blocknames_flaglines, f)
	end

	#elk2wannier
	
	flags = SymAnyDict()
	for (k, v) in blocknames_flaglines
		merge!(flags, v)
	end
	for i in inputs
		if isempty(i.flags)
			i.flags = deepcopy(flags)
		end
	end
	return inputs, structure
end

function find_data(flags, blocknames_flaglines)
	data    = InputData[]
	for f in flags
		if haskey(blocknames_flaglines, f)
			push!(data, InputData(f, :nothing, blocknames_flaglines[f]))
		end
	end
	return data
end

function input_from_task(task, blocknames_flaglines, dir, execs, run)
	if task ∈ ["0", "1", "2", "3", "5"]
		data = find_data((:ngridk, :vkloff), blocknames_flaglines)
	elseif task ∈ ["20", "21"]
		data = find_data((:plot1d, :plot2d, :plot3d), blocknames_flaglines)
	end
	return DFInput{Elk}(task, dir, SymAnyDict(), data, execs, run)
end

function parse(::Type{UnitRange{Int}}, l::AbstractString)
	for delim in [":", ",", "-"]
		if occursin(delim, l)
			return UnitRange(parse.(Int, split(l, delim))...)
		end
	end
end

function parse_block_flags(blockname::Symbol, lines::Vector{<:AbstractString})
	flags = SymAnyDict()
	i = 1
	totlen = length(lines)
	for flaginfo in elk_block_info(blockname).flags
		flagname = flaginfo.name 
		i > length(lines) && break
		typ = eltype(flaginfo)
		if typ == String || typ == Vector{String}
			flags[flagname] = lines[i]
			i += 1
		elseif typ == Symbol
			flags[flagname] = Symbol(lines[i])
			i += 1
		else
			eltyp = eltype(typ)
			# The assumption is that the infinite repeating values are at the end of a block
			if flaginfo == last(elk_block_info(blockname).flags)
				typelen = length(typ)
				parsed = parse.(eltyp, lines[i:end])
				flags[flagname] = typ[]
				@assert length(parsed) % typelen == 0
				for j = 1:typelen:length(parsed)
					push!(flags[flagname], typ(parsed[j:j+typelen-1]...))
				end
			elseif typ <: UnitRange
				flags[flagname] = parse(typ, lines[i])
				i += 1
			else
				typelen = length(typ)
				parsed = parse.(eltyp, lines[i:i + typelen - 1])
				flags[flagname] = length(parsed) == 1 ? parsed[1] : typ(parsed)
				i += typelen 
			end
		end
	end
	return flags
end

function elk_write_structure(f, structure)
	write(f, "atoms\n\t")
	species = unique(map(x->x.name, atoms(structure)))
	nspecies = length(species)
	write(f, "$nspecies\n")
	for s in species
		write(f, "\t'$s.in'\n")
		ats = filter(x->x.name == s, atoms(structure))
		write(f, "\t$(length(ats))\n")
		for a in ats
			write(f, "\t")
			for i in inv(cell(structure)')*position(a)
				write(f, "$(round(i, digits=8)) ")
			end
			for i in magnetization(a)
				write(f, "$i ")
			end
			write(f, "\n")
		end
		write(f, "\n")
	end
	write(f, "\n")
	write(f, "avec\n\t")

	for i in eachindex(structure.cell)
		c = structure.cell[i]
		if i%3 == 0
			write(f, "$c\n\t")
		else
			write(f, "$c ")
		end
	end
	write(f, "\n")
end

function construct_flag_blocks(inputs::Vector{DFInput{Elk}})
	all_flags = SymAnyDict()
	for i in inputs
		merge!(all_flags, i.flags)
		for d in data(i)
			merge!(all_flags, d.data)
		end
	end
	block_flags = SymAnyDict()

	for b in ELK_CONTROLBLOCKS
		bvals = Any[]
		for f in b.flags
			if haskey(all_flags, f.name)
				push!(bvals, pop!(all_flags, f.name))
			end
		end
		if !isempty(bvals)
			block_flags[b.name] = bvals
		end
	end
	return block_flags
end

function elk_write_flag_val(f, val)
	write(f, "\t")
	if isa(val, UnitRange)
		write(f, "$(first(val))-$(last(val))")
	elseif isa(val, AbstractVector)
		for i in val
			write(f, "$i ")
		end
	else
		if isa(val, AbstractString)
			write(f, "'$val'")
		else
			write(f, "$val")
		end
	end
	write(f, "\n")
end

function save(inputs::Vector{DFInput{Elk}}, structure::Structure)
	tasks = map(x -> x.run ? "$(x.name)" : "!$(x.name)", filter(x->x.name != "elk2wannier", inputs))
	elk2wan_input = getfirst(x->x.name == "elk2wannier", inputs)
	if elk2wan_input != nothing && haskey(elk2wan_input.flags, :elk2wan_tasks)
		append!(tasks, pop!(elk2wan_input.flags, :elk2wan_tasks))
	end
	flags = construct_flag_blocks(inputs)
	open(joinpath(inputs[1].dir, "elk.in"), "w") do f
		elk_write_structure(f, structure)
		write(f, "tasks\n")
		for t in tasks
			write(f, "\t$t\n")
		end
		write(f, "\n")
		for (k, v) in flags
			write(f, "$k\n")
			for v_ in values(v)
				if isa(v_, Vector{<:AbstractVector})
					for v__ in v_
						elk_write_flag_val(f, v__)
					end
				else
					elk_write_flag_val(f, v_)
				end
			end
			write(f, "\n")
		end
	end
end
