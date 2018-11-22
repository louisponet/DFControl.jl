import Base: parse

const QEEXECS = [
    "pw.x",
    "projwfc.x",
    "pp.x"
]

#this is all pretty hacky with regards to the new structure and atom api. can for sure be a lot better!
"Quantum espresso card option parser"
cardoption(line) = Symbol(match(r"((?:[a-z][a-z0-9_]*))", split(line)[2]).match)

"""
    read_qe_output(filename::String, T=Float64)

Reads a generic quantum espresso input, returns a dictionary with all found data in the file.
Possible keys:
 - `:fermi`
 - `:polarization`
 - `:pol_mod`
 - `:k_cryst`
 - `:k_cart`
 - `:alat`
 - `:cell_parameters`
 - `:pos_option`
 - `:atomic_positions`
 - `:total_force`
 - `:colin_mag_moments`
 - `:bands`
 - `:accuracy`
"""
function read_qe_output(filename::String, T=Float64)
    out = Dict{Symbol,Any}()
    open(filename, "r") do f
        prefac_k     = nothing
        k_eigvals    = Array{Array{T,1},1}()
        lowest_force = T(1000000)

        while !eof(f)
            line = readline(f)

            #polarization
            if occursin("C/m^2", line)
                s_line = split(line)
                P      = parse(T, s_line[3])
                mod    = parse(T, s_line[5][1:end-1])
                readline(f)
                s_line = parse.(T, split(readline(f))[6:2:10])
                out[:polarization] = Point3{T}(P * s_line[1], P * s_line[2], P * s_line[3])
                out[:pol_mod]      = mod

                #PseudoPot
            elseif occursin("PseudoPot", line)
                !haskey(out, :pseudos) && (out[:pseudos] = Dict{Symbol, String}())
                pseudopath = readline(f) |> strip |> splitdir
                out[:pseudos][Symbol(split(line)[5])] = pseudopath[2]
                !haskey(out, :pseudodir) && (out[:pseudodir] = pseudopath[1])
                #fermi energy
            elseif occursin("Fermi", line)
                out[:fermi]        = parse(T, split(line)[5])
            elseif occursin("lowest unoccupied", line) && occursin("highest occupied", line)
                out[:fermi]        = parse(T, split(line)[7])

            elseif occursin("lowest unoccupied", line) || occursin("highest occupied", line)
                out[:fermi]        = parse(T, split(line)[5])
                #setup for k_points
            elseif occursin("celldm(1)", line)
                alat_bohr = parse(T, split(line)[2])
                prefac_k  = T(2pi / alat_bohr * 1.889725)
                #k_cryst

            elseif occursin("cryst.", line) && length(split(line)) == 2
                out[:k_cryst] = Vector{Vec3{T}}()
                line = readline(f)
                while line != "" && !occursin("--------", line)
                    push!(out[:k_cryst], parse_k_line(line, T))
                    line = readline(f)
                end

                #k_cart
            elseif occursin("cart.", line) && length(split(line)) == 5
                out[:k_cart] = Vector{Vec3{T}}()
                line = readline(f)
                while line != "" && !occursin("--------", line)
                    push!(out[:k_cart], prefac_k * parse_k_line(line, T))
                    line = readline(f)
                end

                #bands
            elseif occursin("k", line) && occursin("PWs)", line)
                tmp = T[]
                readline(f)
                line = readline(f)
                while line != "" && !occursin("--------", line)
                    append!(tmp, parse_line(T, line))
                    line = readline(f)
                end
                push!(k_eigvals, tmp)

                #errors
            elseif occursin("mpirun noticed", line)
                @warn "File ended unexpectedly, returning what info has been gathered so far."
                return out
                break
                #vcrel outputs
            elseif occursin("Begin final coordinates", line)
                line = readline(f)
                while !occursin("End final coordinates", line)
                    if occursin("CELL_PARAMETERS", line)
                        out[:alat]            = occursin("angstrom", line) ? :angstrom : parse(T, split(line)[end][1:end-1])
                        out[:cell_parameters] = reshape(T[parse.(T, split(readline(f))); parse.(T, split(readline(f))); parse.(T, split(readline(f)))], (3,3))'
                    elseif occursin("ATOMIC_POSITIONS", line)
                        out[:pos_option]      = cardoption(line)
                        line  = readline(f)
                        atoms = []
                        while !occursin("End", line)
                            s_line = split(line)
                            key    = Symbol(s_line[1])
                            push!(atoms, key=>Point3{T}(parse.(T, s_line[2:end])...))
                            line = readline(f)
                        end
                        posdict = Dict{Symbol, Vector{Point3{T}}}()
                        for (atsym, pos) in atoms
                            if haskey(posdict, atsym)
                                push!(posdict[atsym], pos)
                            else
                                posdict[atsym] = [pos]
                            end
                        end
                        out[:atomic_positions] = posdict
                        break
                    end
                    line = readline(f)
                end
                pseudo_data = InputData(:atomic_species, :none, out[:pseudos])
                tmp_flags = Dict(:ibrav => 0, :A => (out[:alat] == :angstrom ? 1 : conversions[:bohr2ang] * out[:alat]))
                cell_data = InputData(:cell_parameters, :alat, Mat3(out[:cell_parameters]))
                atoms_data = InputData(:atomic_positions, out[:pos_option], out[:atomic_positions])
                out[:final_structure] = extract_structure!("newstruct", tmp_flags, cell_data, atoms_data, pseudo_data)

            elseif occursin("Total force", line)
                force = parse(T, split(line)[4])
                if force <= lowest_force
                    lowest_force      = force
                    out[:total_force] = force
                end
            elseif occursin("Magnetic moment per site", line)
                key = :colin_mag_moments
                out[key] = T[]
                line = readline(f)
                while !isempty(line)
                    push!(out[key], parse(Float64, split(line)[6]))
                    line = readline(f)
                end
            elseif occursin("estimated scf accuracy", line)
                key = :accuracy
                acc = parse(T, split(line)[end-1])
                if haskey(out, key)
                    push!(out[key], acc)
                else
                    out[key] = [acc]
                end
            end
        end

        #process bands
        if !isempty(k_eigvals)
            out[:bands] = Vector{DFBand{T}}()
            for i1=1:length(k_eigvals[1])
                eig_band = T[]
                for i = 1:length(out[:k_cart])
                    push!(eig_band, k_eigvals[i][i1])
                end
                push!(out[:bands], DFBand(get(out, :k_cart,[zero(Vec3)]), get(out, :k_cryst, [zero(Vec3)]), eig_band))
            end
        end
        return out
    end
end

"""
    read_qe_bands(filename::String, T=Float64)

Reads the output file of a 'bands' calculation in Quantum Espresso.
Returns an array of DFBands each with the same k_points and their respective energies.
"""
read_qe_bands_file(filename::String, T=Float64) = read_qe_output(filename, T)[:bands]

"""
    read_ks_from_qe_output(filename::String, T=Float64)

Read k-points from a Quantum Espresso bands output file in cartesian (2pi/alat in Angstrom^-1!) and crystalline coordinates.
Returns (k_points_cart,k_points_cryst).
"""
function read_ks_from_qe_output(filename::String, T=Float64)
    t = read_qe_output(filename, T)
    return t[:k_cart], t[:k_cryst]
end

"""
    read_fermi_from_qe_output(filename::String,T=Float64)

Reads the Fermi level from a Quantum Espresso scf calculation output file
(if there is one).
"""
read_fermi_from_qe_output(filename::String, T=Float64) = read_qe_output(filename, T)[:fermi]

"""
    read_qe_kpdos(filename::String,column=1;fermi=0)

Reads the k_resolved partial density of states from a Quantum Espresso projwfc output file.
Only use this if the flag kresolveddos=true in the projwfc input file!! The returned matrix can be readily plotted using heatmap() from Plots.jl!
Optional input: column = 1 (column of the output, 1 = first column after ik and E)
fermi  = 0 (possible fermi offset of the read energy values)
Return:         Array{Float64,2}(length(k_points),length(energies)) ,
(ytickvals,yticks)
"""
function read_qe_kpdos(filename::String, column=1; fermi=0)
    read_tmp = readdlm(filename)
    zmat     = zeros(typeof(read_tmp[1]), Int64(read_tmp[end, 1]), size(read_tmp)[1] / Int64(read_tmp[end, 1]))
    for i1 = 1:size(zmat)[1]
        for i2 = 1:size(zmat)[2]
            zmat[i1, i2] = read_tmp[size(zmat)[2] * (i1 - 1) + i2, 2 + column]
        end
    end

    yticks    = collect(Int64(div(read_tmp[1, 2] - fermi, 1)):1:Int64(div(read_tmp[end, 2] - fermi, 1)))
    ytickvals = [findfirst(x -> norm(yticks[1] + fermi - x) <= 0.1, read_tmp[:, 2])]
    for (i, tick) in enumerate(yticks[2:end])
        push!(ytickvals, findnext(x -> norm(tick + fermi - x) <= 0.1, read_tmp[:, 2], ytickvals[i]))
    end

    return  zmat', (ytickvals, yticks)
end

"""
    read_qe_pdos(filename::String, column=1; fermi=0)

Reads partial dos file. One can specify the column of values to read.
"""
function read_qe_pdos(filename::String, column=1; fermi=0)
    read_tmp = readdlm(filename)
    energies = read_tmp[:,1] .- fermi
    values   = read_tmp[:,1+column]

    return energies, values
end

"""
    read_qe_polarization(filename::String, T=Float64)

Returns the polarization and modulus.
"""
function read_qe_polarization(filename::String, T=Float64)
    t = read_qe_output(filename, T)
    return t[:polarization], t[:pol_mod]
end

read_qe_vcrel(filename::String, T=Float64) = read_qe_output(filename, T) do x
                                                return x[:cell_parameters], x[:alat], x[:atomic_positions], x[:pos_option]
                                            end

function alat(flags, pop=false)
    if haskey(flags, :A)
        alat = pop ? pop!(flags, :A) : flags[:A]
    elseif haskey(flags, :celldm_1)
        alat = pop ? pop!(flags, :celldm_1) : flags[:celldm_1]
        alat *= conversions[:bohr2ang]
    else
        error("Cell option 'alat' was found, but no matching flag was set. \n
               The 'alat' has to  be specified through 'A' and 'celldm(1)'.")
    end
    return alat
end

#TODO handle more fancy cells
function extract_cell!(flags, cell_block)
    if cell_block != nothing
        _alat = 1.0
        if cell_block.option == :alat
            @assert pop!(flags, :ibrav) == 0 "Only ibrav = 0 allowed for now."
            _alat = alat(flags)

        elseif cell_block.option == :bohr
            _alat = conversions[:bohr2ang]
        end

        return _alat * cell_block.data
    end
end

function extract_atoms!(control, atom_block, pseudo_block, cell)
    atoms = Atom{Float64}[]

    option = atom_block.option
    if option == :crystal || option == :crystal_sg
        primv = cell
    elseif option == :alat
        primv = alat(control, true) * Mat3(Matrix(1.0I, 3, 3))
    elseif option == :bohr
        primv = conversions[:bohr2ang] * Mat3(Matrix(1.0I, 3, 3))
    else
        primv = Mat3(Matrix(1.0I, 3, 3))
    end

    for (at_sym, positions) in atom_block.data
        pseudo = haskey(pseudo_block.data, at_sym) ? pseudo_block.data[at_sym] : error("Please specify a pseudo potential for atom '$at_sym'.")
        for pos in positions
            push!(atoms, Atom(at_sym, element(at_sym), primv' * pos, pseudo=pseudo))
        end
    end

    return atoms
end

function extract_structure!(name, control, cell_block, atom_block, pseudo_block)
    if atom_block == nothing
        return nothing
    end
    cell = extract_cell!(control, cell_block)
    atoms = extract_atoms!(control, atom_block, pseudo_block, cell)
    return Structure(name, cell, atoms)
end

"""
    read_qe_input(filename, T=Float64; exec="pw.x",  runcommand="", run=true, structure_name="NoName")

Reads a Quantum Espresso input file. The exec get's used to find which flags are allowed in this input file, and convert the read values to the correct Types.
Returns a `DFInput{QE}` and the `Structure` that is found in the input.
"""
function read_qe_input(filename; execs=[Exec("pw.x")], run=true, structure_name="noname")
    @assert ispath(filename) "$filename is not a valid path."
    lines = read(filename) |>
        String |>
        x -> replace(x, ", " => "\n") |>
        x -> replace(x, "," => " ") |>
        x -> split(x, "\n") .|>
        strip |>
        x -> filter(!isempty, x) |>
        x -> filter(y -> !occursin("&", y), x) |>
        x -> filter(y -> !(occursin("/", y) && length(y) == 1), x) |>
        x -> filter(y -> y[1] != '!', x)

    exec = getfirst(x->x.exec ∈ QEEXECS, execs)
    flaglines = strip_split.(filter(x -> occursin("=", x), lines), "=")
    parsed_flags = Dict()
    #easy flags
    for (f, v) in filter(x-> !occursin("(", x[1]), flaglines)
        sym = Symbol(f)
        typ = flagtype(QE, exec, sym)
        if eltype(typ) <: Bool
            v = strip(lowercase(v), '.')
        elseif eltype(typ) <: Number
            v = replace(v, "d" => "e")
        end
        tval = typ != String ? parse.((typ,), split(v)) : v
        parsed_flags[sym] = length(tval) == 1 ? tval[1] : tval
    end
    nat  = parsed_flags[:nat]
    ntyp = parsed_flags[:ntyp]
    #difficult flags
    for (f, v) in filter(x-> occursin("(", x[1]), flaglines)
        _s = split(replace(replace(replace(f, "(" => " "), ")" => " "), "," => " "))

        sym = Symbol(_s[1])
        ids = parse.(Int, _s[2:end])
        typ = flagtype(QE, exec, sym)
        v = replace(v, "d" => "e")
        parsedval = parse.((eltype(typ),), split(v))
        if !haskey(parsed_flags, sym)
            if typ <: AbstractMatrix
                parsed_flags[sym] = length(parsedval) == 1 ? zeros(eltype(typ), ntyp, 10) : fill(zeros(eltype(typ), length(parsedval)), ntyp, 10) #arbitrary limit
            else
                parsed_flags[sym] = length(parsedval) == 1 ? zeros(eltype(typ), ntyp) : fill(zeros(eltype(typ), length(parsedval)), ntyp)
            end
        end
        if length(ids) == 1
            parsed_flags[sym][ids[1]] = length(parsedval) == 1 ? parsedval[1] : parsedval
        else
            parsed_flags[sym][ids[1], ids[2]] = length(parsedval) == 1 ? parsedval[1] : parsedval
        end
    end

    findcard(s) = findfirst(l -> occursin(s, lowercase(l)), lines)
    i = findcard("atomic_species")
    pseudos = InputData(:atomic_species, :none, Dict{Symbol, String}())
    for k=1:ntyp
        sline = strip_split(lines[i+k])
        pseudos.data[Symbol(sline[1])] = sline[end]
    end

    i = findcard("cell_parameters")
    cell_block = InputData(:cell_parameters,
                           cardoption(lines[i]),
                           Mat3([parse(Float64, split(lines[i+k])[j]) for k=1:3, j=1:3]))

    i = findcard("atomic_positions")
    atom_block = InputData(:atomic_positions,
                           cardoption(lines[i]),
                           Dict{Symbol, Vector{Point3{Float64}}}() )
    for k=1:nat
        sline = split(lines[i+k])
        atsym = Symbol(sline[1])
        point = Point3(parse.(Float64, sline[2:4]))
        if !haskey(atom_block.data, atsym)
            atom_block.data[atsym] = [point]
        else
            push!(atom_block.data[atsym], point)
        end
    end

    datablocks = InputData[]
    i = findcard("k_points")
    if i!=nothing
        k_option = cardoption(lines[i])
        if k_option == :automatic
            s_line = split(lines[i+1])
            k_data = parse.(Int, s_line)
        else
            nks    = parse(Int, lines[i+1])
            k_data = Vector{Vector{Float64}}(undef, nks)
            for k = 1:nks
                k_data[k] = parse.(Float64, split(lines[i+1+k]))
            end
        end
        push!(datablocks, InputData(:k_points, k_option, k_data))
    end

    structure = extract_structure!(structure_name, parsed_flags, cell_block, atom_block, pseudos)
    delete!.((parsed_flags,), [:ibrav, :nat, :ntyp, :A, :celldm_1, :celldm])
    dir, file = splitdir(filename)
    return DFInput{QE}(splitext(file)[1], dir, parsed_flags, datablocks, execs, run), structure
end

function qe_writeflag(f, flag, value)
    if isa(value, Matrix)
        for i=1:size(value)[1], j=1:size(value)[2]
            if !iszero(value[i,j])
                write(f, "  $flag($i,$j) = $(value[i, j])\n")
            end
        end
    elseif isa(value, Vector)
        for i=1:length(value)
            if !iszero(value[i])
                if length(value[i]) == 1
                    write(f, "  $flag($i) = $(value[i])\n")
                else
                    write(f, "  $flag($i) =")
                    for v in value[i]
                        write(f, " $v")
                    end
                    write(f, "\n")
                end
            end
        end
    else
        write(f, "  $flag = $value\n")
    end
end


"""
    save(input::DFInput{QE}, structure, filename::String=inpath(input))

Writes a Quantum Espresso input file.
"""
function save(input::DFInput{QE}, structure, filename::String=inpath(input))
    if haskey(flags(input), :calculation)
        input[:calculation] = replace(input[:calculation], "_" => "-")
    end
    open(filename, "w") do f
        writeflag(flag_data) = qe_writeflag(f, flag_data[1], flag_data[2])
        write_dat(data)       = write_data(f, data)

        controls = Dict{Symbol, Dict{Symbol, Any}}()

        for (flag, val) in input.flags
            block, variable = qe_block_variable(input, flag)
            if !haskey(controls, block.name)
                controls[block.name] = Dict{Symbol, Any}()
            end
            controls[block.name][flag] = val
        end

        #Here we try to figure out the correct order of the control blocks
        # first we find the order of the pw.x inputs, the rest should follow.
        blocks2file = []
        for name in [:control, :system, :electrons, :ions, :cell]
            push!(blocks2file, name => pop!(controls, name, nothing))
        end
        for name in keys(controls)
            push!(blocks2file, name => pop!(controls, name, nothing))
        end
        filter!(x->x[2]!=nothing, blocks2file)
        for (name, flags) in blocks2file
            write(f, "&$name\n")
            if name == :system
                nat   = length(structure.atoms)
                ntyp  = length(unique(structure.atoms))
                # A     = 1.0
                ibrav = 0
                write(f,"  ibrav = $ibrav\n")
                # write(f,"  A = $A\n")
                write(f,"  nat = $nat\n")
                write(f,"  ntyp = $ntyp\n")
            end
            map(writeflag, [(flag, data) for (flag, data) in flags])
            write(f, "/\n\n")
        end

        write_structure(f, input, structure)
        for dat in input.data
            if dat.option != :none
                write(f, "$(uppercase(String(dat.name))) ($(dat.option))\n")
            else
                write(f, "$(uppercase(String(dat.name)))\n")
            end
            if dat.name == :k_points && dat.option != :automatic
                write(f, "$(length(dat.data))\n")
                write_dat(dat.data)
            else
                write_dat(dat.data)
            end
            write(f, "\n")
        end
    end
end

function write_structure(f, input::DFInput{QE}, structure)
    unique_at = unique(structure.atoms)
    pseudo_lines = String[]
    atom_lines   = String[]
    for at in unique_at
        push!(pseudo_lines, "$(id(at)) $(element(at).atomic_weight)   $(pseudo(at))\n")
    end
    for at in structure.atoms
        pos = position(at)
        push!(atom_lines, "$(id(at))  $(pos[1]) $(pos[2]) $(pos[3])\n")
    end

    write(f, "ATOMIC_SPECIES\n")
    write.((f, ), pseudo_lines)

    write(f, "\n")
    write(f, "CELL_PARAMETERS (angstrom)\n")
    write_cell(f, structure.cell)
    write(f, "\n")

    write(f, "ATOMIC_POSITIONS (angstrom) \n")
    write.((f, ), atom_lines)
    write(f, "\n")
end

function qe_generate_pw2waninput(input::DFInput{Wannier90}, qeprefix, runexecs)
    flags = Dict()
    flags[:prefix] = qeprefix
    flags[:seedname] = name(input)
    flags[:outdir] = "'$(dir(input))'"
    flags[:wan_mode] = "'standalone'"
    flags[:write_mmn] = true
    flags[:write_amn] = true
    if flag(input, :spin) != nothing
        flags[:spin_component] = flag(input, :spin)
    end
    if flag(input, :wannier_plot) != nothing
        flags[:write_unk] = flag(input, :wannier_plot)
    end
    if any(flag(input, :berry_task) .== ("morb", "'morb'"))
        flags[:write_uHu] = true
    end
    pw2wanexec = Exec("pw2wannier90.x", runexecs[2].dir)
    return DFInput{QE}("pw2wan_$(flags[:seedname])", dir(input), flags, InputData[], [runexecs[1], pw2wanexec], input.run)
end
