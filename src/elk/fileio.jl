function read_atoms(fn::IO)
    sane_readline() = strip(strip_split(readline(fn), ":")[1], '\'')
    nspecies = parse(Int, sane_readline())
    species_counter = 0

    reading_atoms = false
    atname = :nothing
    natoms = 0
    atcounter = 0
    positions = Point3{Float64}[]
    bfcmts = Vec3{Float64}[]
    atoms = Tuple{Symbol,Point3{Float64},Vec3{Float64}}[]
    while species_counter < nspecies
        l = sane_readline()
        if !reading_atoms
            atname = Symbol(split(l, ".")[1])
            natoms = parse(Int, sane_readline())
            reading_atoms = true
        elseif atcounter < natoms
            sline = split(l)
            push!(positions, Point3{Float64}(parse.(Float64, sline[1:3])))
            push!(bfcmts, Vec3{Float64}(parse.(Float64, sline[4:end])))
            atcounter += 1
        elseif isempty(l)
            for (p, mag) in zip(positions, bfcmts)
                push!(atoms, (atname, p, mag))
            end
            atname          = :nothing
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

function elk_parse_DFTU(atoms::Vector{<:Atom}, blocknames_flaglines::Dict{Symbol,Any})
    lines = pop!(blocknames_flaglines, Symbol("dft+u"))
    t_l = split(lines[1])
    species_start = 3
    if length(t_l) == 1
        dftu_type = parse(Int, t_l[1])
        inpdftu   = parse(Int, lines[2])
    else
        dftu_type, inpdftu = parse.(Int, t_l)
        species_start = 2
    end
    species_names = name.(unique(atoms))
    for line in lines[species_start:end]
        sline         = split(line)
        species_id, l = parse.(Int, sline[1:2])
        U, J0         = parse.(Float64, sline[3:4])
        for at in filter(x -> name(x) == species_names[species_id], atoms)
            at.dftu = DFTU{Float64}(; U = U, J0 = J0, l = l)
        end
    end
    return dftu_type, inpdftu
end

"""
    elk_read_calculation(filename; execs=[Exec(exec="elk")], run=true, structure_name="noname")

Reads an Elk calculation file. The `ELK_EXEC` inside execs gets used to find which flags are allowed in this calculation file, and convert the read values to the correct Types.
Returns a `DFCalculation{Elk}` and the `Structure` that is found in the calculation.
"""
function elk_read_calculation(fn::String; execs = [Exec(; exec = "elk")], run = true,
                              structure_name = "noname")
    blocknames_flaglines = Dict{Symbol,Any}()
    atoms = Atom{Float64}[]
    calculations = DFCalculation[]
    dir, file = splitdir(fn)

    #TODO it's a bit messy that the atoms block is not handled like this, not sure if we want to clean this up
    #blocks that require specific parsing
    special_blocks = (:avec, :tasks, Symbol("dft+u"))

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
                    @warn "Please supply the .win file instead of wannierExtra block for extracting the Wannier DFCalculation."

                    wflags, wdata, ab, cb, proj_block = wan_read_calculation(f)
                    continue
                    # push!(calculations, DFCalculation{Wannier90}("wannier", dir, wflags, wdata, execs, run))
                elseif info != nothing
                    blocklines = String[]
                    line = sane_readline()
                    while !isempty(line)
                        push!(blocklines,
                              strip(replace(strip(split(line, ':')[1]), "'" => ""), '.'))
                        line = sane_readline()
                    end
                    if !in(blockname, special_blocks)
                        blocknames_flaglines[blockname] = parse_block_flags(blockname,
                                                                            collect(Iterators.flatten(split.(filter(x -> !isempty(x),
                                                                                                                    blocklines)))))
                    else
                        blocknames_flaglines[blockname] = filter(x -> !isempty(x),
                                                                 blocklines)
                    end
                end
            end
        end
    end
    # for i in 
    #cell
    scale = haskey(blocknames_flaglines, :scale) ? last(pop!(blocknames_flaglines[:scale])) : 1.0
    cell  = uconvert.(Ang, Mat3(reshape(scale .* parse.(Float64, collect(Iterators.flatten(split.(pop!(blocknames_flaglines, :avec))))), 3, 3) .* 1a₀))'

    newats = [Atom{Float64,eltype(cell)}(; name = x[1], element = element(x[1]),
                                         position_cryst = x[2], position_cart = cell * x[2],
                                         magnetization = x[3]) for x in atoms]
    #structure #TODO make positions be in lattices coordinates
    structure = Structure(structure_name, cell, newats)

    #different tasks
    tasks_     = pop!(blocknames_flaglines, :tasks)
    should_run = map(x -> !occursin("!", x), tasks_)
    tasks      = map(x -> strip(strip(x, '!')), tasks_)

    #The wan tasks can prbably be skipped completely and added depending on the wannier calculation
    wan_tasks = String[]
    for (x, run) in zip(tasks, should_run)
        if x ∈ ["601", "602", "603", "604", "605"]
            push!(wan_tasks, x)
        else
            push!(calculations,
                  calculation_from_task(x, blocknames_flaglines, dir, execs, run))
        end
    end
    if haskey(blocknames_flaglines, :wannier) && !isempty(wan_tasks)
        flags = pop!(blocknames_flaglines, :wannier)
        # flags[:elk2wan_tasks] = wan_tasks
        push!(calculations,
              DFCalculation{Elk}(; name = "elk2wannier", dir = dir, flags = flags,
                                 execs = execs, run = true))
    end
    for f in (:ngrid, :vkloff, :plot1d, :plot2d, :plot3d)
        haskey(blocknames_flaglines, f) && pop!(blocknames_flaglines, f)
    end

    #elk2wannier

    flags = SymAnyDict()

    if haskey(blocknames_flaglines, Symbol("dft+u"))
        flags[:dftu], flags[:inpdftu] = elk_parse_DFTU(structure.atoms,
                                                       blocknames_flaglines)
    end

    for (k, v) in blocknames_flaglines
        merge!(flags, v)
    end

    for i in calculations
        if isempty(i.flags)
            i.flags = deepcopy(flags)
        end
    end
    return calculations, structure
end

function find_data(flags, blocknames_flaglines)
    data = InputData[]
    for f in flags
        if haskey(blocknames_flaglines, f)
            push!(data, InputData(f, :nothing, blocknames_flaglines[f]))
        end
    end
    return data
end

function calculation_from_task(task, blocknames_flaglines, dir, execs, run)
    if task ∈ ["0", "1", "2", "3", "5"]
        data = find_data((:ngridk, :vkloff), blocknames_flaglines)
    elseif task ∈ ["20", "21"]
        data = find_data((:plot1d, :plot2d, :plot3d), blocknames_flaglines)
    end
    return DFCalculation{Elk}(name = task, dir = dir, data = data, execs = execs, run = run)
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
                for j in 1:typelen:length(parsed)
                    push!(flags[flagname], typ(parsed[j:j+typelen-1]...))
                end
            elseif typ <: UnitRange
                flags[flagname] = parse(typ, lines[i])
                i += 1
            else
                typelen = length(typ)
                parsed = parse.(eltyp, lines[i:i+typelen-1])
                flags[flagname] = length(parsed) == 1 ? parsed[1] : typ(parsed)
                i += typelen
            end
        end
    end
    return flags
end

function elk_write_structure(f, structure)
    write(f, "atoms\n\t")
    species = unique(map(x -> x.name, atoms(structure)))
    nspecies = length(species)
    write(f, "$nspecies\n")
    for s in species
        write(f, "\t'$s.in'\n")
        ats = filter(x -> x.name == s, atoms(structure))
        write(f, "\t$(length(ats))\n")
        for a in ats
            write(f, "\t")
            for i in position_cryst(a)
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

    for i in eachindex(structure.cell')
        c = ustrip(uconvert(a₀, structure.cell[i]))
        if i % 3 == 0
            write(f, "$c\n\t")
        else
            write(f, "$c ")
        end
    end
    return write(f, "\n")
end

function construct_flag_blocks(calculations::Vector{DFCalculation{Elk}})
    all_flags = SymAnyDict()
    for i in calculations
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
    return write(f, "\n")
end

function elk_write_DFTU(f::IO, structure::AbstractStructure, dftu_vals::Vector)
    write(f, "dft+u\n")
    for v in dftu_vals
        elk_write_flag_val(f, v)
    end
    species = unique(atoms(structure))
    for (i, s) in enumerate(species)
        write(f, "\t$i $(string(Elk, dftu(s)))\n")
    end
    return write(f, "\n")
end

function save(calculations::Vector{DFCalculation{Elk}}, structure::Structure)
    tasks = map(x -> x.run ? "$(x.name)" : "!$(x.name)",
                filter(x -> x.name != "elk2wannier", calculations))
    elk2wan_calculation = getfirst(x -> x.name == "elk2wannier", calculations)
    if elk2wan_calculation != nothing && haskey(elk2wan_calculation.flags, :elk2wan_tasks)
        append!(tasks, pop!(elk2wan_calculation.flags, :elk2wan_tasks))
    end
    flags = construct_flag_blocks(calculations)
    open(joinpath(calculations[1].dir, "elk.in"), "w") do f
        write(f, "tasks\n")
        for t in tasks
            write(f, "\t$t\n")
        end
        write(f, "\n")

        elk_write_structure(f, structure)

        haskey(flags, Symbol("dft+u")) &&
            elk_write_DFTU(f, structure, pop!(flags, Symbol("dft+u")))

        #remaining generic blocks
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
