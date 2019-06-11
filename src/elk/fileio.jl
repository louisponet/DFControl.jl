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
	atoms     = Atom[]
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


function elk_read_input(fn::String)
	blocknames_flaglines = Dict{Symbol, Vector{String}}()
	atoms = Atom[]
	open(fn, "r") do f
		while !eof(f)
			line = strip(readline(f))
			isempty(line) || line[1] == '!' && continue
			if !isempty(line)
				blockname = Symbol(line)
				info = elk_block_info(blockname)
				if blockname == :atoms
					atoms = read_atoms(f)
				elseif info != nothing
					blocklines = String[]
					line = readline(f)
					while !isempty(line)
						push!(blocklines, replace_multiple(strip(split(line, ':')[1]), "'" => "", "."=>""))
						line = readline(f)
					end
					blocknames_flaglines[blockname] = collect(Iterators.flatten(split.(filter(x->!isempty(x), blocklines))))
				end
			end
		end
	end
	flags = SymAnyDict()
	for (k, v) in blocknames_flaglines
		i = 1
		for flaginfo in elk_block_info(k).flags
			i > length(v) && break
			typ = eltype(flaginfo)
			if typ == String || typ == Vector{String}
				flags[flaginfo.name] = v[i]
				i += 1
			else
				parsed = parse.(eltype(typ), v[i:i + length(typ) - 1])
				flags[flaginfo.name] = length(parsed) == 1 ? parsed[1] : typ(parsed)
				i += length(typ)
			end
		end
	end
	@show flags

	@show blocknames_flaglines[:tasks]
end
