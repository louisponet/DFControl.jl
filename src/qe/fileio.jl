const parseable_qe_execs = ["pw.x", "projwfc.x", "pw2wannier90.x", "pp.x"]
#this is all pretty hacky with regards to the new structure and atom api. can for sure be a lot better!
"Quantum espresso card option parser"
function card_option(line)
    if contains(line, "{")
        return Symbol(strip(split(line, "{")[end], '}'))
    elseif contains(line, "(")
        return Symbol(strip(split(line, "(")[end], ')'))
    else
        return Symbol(split(line)[end])
    end
end

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
            if contains(line, "C/m^2")
                s_line = split(line)
                P      = parse(T, s_line[3])
                mod    = parse(T, s_line[5][1:end-1])
                readline(f)
                s_line = parse.(T, split(readline(f))[6:2:10])
                out[:polarization] = Point3{T}(P * s_line[1], P * s_line[2], P * s_line[3])
                out[:pol_mod]      = mod

                #fermi energy
            elseif contains(line, "Fermi")
                out[:fermi]        = parse(T, split(line)[5])
            elseif contains(line, "lowest unoccupied") && contains(line, "highest occupied")
                out[:fermi]        = parse(T, split(line)[7])

            elseif contains(line, "lowest unoccupied") || contains(line, "highest occupied")
                out[:fermi]        = parse(T, split(line)[5])
                #setup for k_points
            elseif contains(line, "celldm(1)")
                alat_bohr = parse(T, split(line)[2])
                prefac_k  = T(2pi / alat_bohr * 1.889725)
                #k_cryst

            elseif contains(line, "cryst.") && length(split(line)) == 2
                out[:k_cryst] = Array{Array{T,1},1}()
                line = readline(f)
                while line != "" && !contains(line, "--------")
                    push!(out[:k_cryst], parse_k_line(line, T))
                    line = readline(f)
                end

                #k_cart
            elseif contains(line, "cart.") && length(split(line)) == 5
                out[:k_cart] = Array{Array{T,1},1}()
                line = readline(f)
                while line != "" && !contains(line, "--------")
                    push!(out[:k_cart], prefac_k * parse_k_line(line, T))
                    line = readline(f)
                end

                #bands
            elseif contains(line, "k") && contains(line, "PWs)")
                tmp = T[]
                readline(f)
                line = readline(f)
                while line != "" && !contains(line, "--------")
                    append!(tmp, parse_line(T, line))
                    line = readline(f)
                end
                push!(k_eigvals, tmp)

                #errors
            elseif contains(line, "mpirun noticed")
                warn("File ended unexpectedly, returning what info has been gathered so far.")
                return out
                break

                #vcrel outputs
            elseif contains(line, "Begin final coordinates")
                line = readline(f)
                while !contains(line, "End final coordinates")

                    if contains(line, "CELL_PARAMETERS")
                        out[:alat]            = parse(T, split(line)[end][1:end-1])
                        out[:cell_parameters] = reshape(T[parse.(T, split(readline(f))); parse.(T, split(readline(f))); parse.(T, split(readline(f)))], (3,3))
                    elseif contains(line, "ATOMIC_POSITIONS")
                        out[:pos_option]      = card_option(line)
                        line  = readline(f)
                        atoms = Dict{Symbol,Array{Point3{T},1}}()
                        while !contains(line, "End")
                            s_line = split(line)
                            key    = Symbol(s_line[1])
                            if key in keys(atoms)
                                push!(atoms[key], Point3{T}(parse.(T, s_line[2:end])...))
                            else
                                atoms[key] = [Point3{T}(parse.(T, s_line[2:end])...)]
                            end
                            line = readline(f)
                        end
                        out[:atomic_positions] = atoms
                        break
                    end
                    line = readline(f)
                end

            elseif contains(line, "Total force")
                force = parse(T, split(line)[4])
                if force <= lowest_force
                    lowest_force      = force
                    out[:total_force] = force
                end
            elseif contains(line, "Magnetic moment per site")
                key = :colin_mag_moments
                out[key] = T[]
                line = readline(f)
                while !isempty(line)
                    push!(out[key], parse(split(line)[6]))
                    line = readline(f)
                end
            end
        end

        #process bands
        if !isempty(k_eigvals)
            out[:bands] = Array{DFBand{T},1}()
            for i1=1:length(k_eigvals[1])
                eig_band = T[]
                for i = 1:length(out[:k_cart])
                    push!(eig_band, k_eigvals[i][i1])
                end
                push!(out[:bands], DFBand(get(out, :k_cart,[[0.0,0.0,0.0]]), get(out, :k_cryst, [[0.0,0.0,0.0]]), eig_band))
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

# function read_errors(filename::String)
#     out = String[]
#     open(filename, "r") do f
#         line = readline(f)
#         while !eof(f)
#             if contains(lowercase(line), "error")
#                 push!(out, line)
#             end
#         end
#         return out
#     end
# end
function alat(control, pop=false)
    sysblock = block(control, :system)
    if sysblock == nothing
        error("Could not resolve the alat!")
    end
    if haskey(sysblock.flags, :A)
        alat = pop ? pop!(sysblock.flags, :A) : sysblock.flags[:A]
    elseif haskey(sysblock.flags, celldm_1)
        alat = pop ? pop!(sysblock.flags, celldm_1) : sysblock.flags[celldm_1]
        alat *= conversions[:bohr2ang]
    else
        error("So far alat can be specified only through A and celldm(1).")
    end
    return alat
end

#TODO handle more fancy cells
function extract_cell!(control, cell_block)
    if cell_block != nothing
        _alat = 1.0
        if cell_block.option == :alat
            @assert pop!(block(control,:system).flags, :ibrav) == 0 "Only ibrav = 0 allowed for now."
            _alat = alat(control)

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
        primv = alat(control, true) * Mat3(eye(3))
    elseif option == :bohr
        primv = conversions[:bohr2ang] * Mat3(eye(3))
    else
        primv = Mat3(eye(3))
    end

    for (at_sym, positions) in atom_block.data
        pseudo = haskey(pseudo_block.data, at_sym) ? pseudo_block.data[at_sym] : error("Please specify a pseudo potential for atom '$at_sym'.")
        for pos in positions
            push!(atoms, Atom(at_sym, element(at_sym), primv' * pos, :pseudo => pseudo))
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
Returns a `QEInput` and the `Structure` that is found in the input.
"""
function read_qe_input(filename, T=Float64::Type; exec=Exec("pw.x"), runcommand=Exec(""), run=true, structure_name="NoName")
    control = Array{QEControlBlock,1}()
    data    = Array{QEDataBlock,1}()
    atom_block     = nothing
    cell_block     = nothing
    pseudo_block   = nothing
    open(filename) do f
        line = readline(f)
        while !eof(f)
            @label start_label
            if contains(line, "&")
                c_block_name    = Symbol(lowercase(strip(strip(line), '&')))
                flag_dict       = Dict{Symbol,Any}()
                def_block_flags = all_qe_block_flags(exec.exec, c_block_name)
                line            = readline(f)
                while strip(line) != "/"
                    if contains(line, "!")
                        line = readline(f)
                        continue
                    end
                    split_line = filter(x -> x != "", strip.(split(line, ",")))
                    for s in split_line
                        key, val = String.(strip.(split(s, "=")))
                        qe_flag  = Symbol(split(key, "(")[1])
                        flag_type = filter(x->x.name == qe_flag, def_block_flags)
                        if length(flag_type) != 0
                            t_val = parse_flag_val(val, flag_type[1].typ)
                            flag_dict[Symbol(key)] = eltype(t_val) == flag_type[1].typ || flag_type[1].typ == String ? t_val : error("Couldn't parse the value of flag '$key' in file '$filename'!")
                        else
                            error("Error reading $filename: flag '$key' not found in QE flag Dictionary for control block $c_block_name !")
                        end
                    end
                    line = readline(f)
                end
                push!(control, QEControlBlock(c_block_name, flag_dict))
                @goto start_label

            elseif contains(line, "CELL_PARAMETERS") || contains(line, "cell_parameters")
                cell_unit    = card_option(line)
                cell         = MMatrix{3, 3, T}()
                cell[1, 1:3] = parse.(T, split(readline(f)))
                cell[2, 1:3] = parse.(T, split(readline(f)))
                cell[3, 1:3] = parse.(T, split(readline(f)))
                line = readline(f)
                cell_block = QEDataBlock(:cell_parameters, cell_unit, cell)
                @goto start_label

            elseif contains(line, "ATOMIC_SPECIES") || contains(line, "atomic_species")
                line    = readline(f)
                pseudos = Dict{Symbol,String}()
                while length(split(line)) == 3
                    pseudos[Symbol(split(line)[1])] = split(line)[end]
                    line = readline(f)
                end
                pseudo_block = QEDataBlock(:atomic_species, :none, pseudos)
                @goto start_label

            elseif contains(line, "ATOMIC_POSITIONS") || contains(line, "atomic_positions")
                option = card_option(line)
                atoms  = Dict{Symbol,Array{Point3{T},1}}()
                line   = readline(f)
                while length(split(line)) == 4
                    s_line   = split(line)
                    atom     = Symbol(s_line[1])
                    position = Point3(parse(T, s_line[2]), parse(T, s_line[3]), parse(T, s_line[4]))
                    if !haskey(atoms, atom)
                        atoms[atom] = [position]
                    else
                        push!(atoms[atom], position)
                    end
                    line = readline(f)
                end
                atom_block = QEDataBlock(:atomic_positions, option, atoms)
                @goto start_label

            elseif contains(line, "K_POINTS") || contains(line, "k_points")
                k_option = card_option(line)
                line     = readline(f)
                if k_option == :automatic
                    s_line = split(line)
                    k_data = parse.(Int, s_line)
                else
                    nks    = parse(Int, line)
                    k_data = Array{Array{T,1},1}(nks)
                    for i = 1:nks
                        k_data[i] = parse.(T, split(readline(f)))
                    end
                end
                push!(data, QEDataBlock(:k_points, k_option, k_data))
                @goto start_label
            end
            line = readline(f)
        end
    end

    structure = extract_structure!(structure_name, control, cell_block, atom_block, pseudo_block)
    sysblock = block(control, :system)
    if sysblock != nothing
        pop!(sysblock.flags, :ibrav, nothing)
        pop!(sysblock.flags, :nat, nothing)
        pop!(sysblock.flags, :ntyp, nothing)
        pop!(sysblock.flags, :A, nothing)
        pop!(sysblock.flags, celldm_1, nothing)
        pop!(sysblock.flags, :celldm, nothing)
    end
    return QEInput(splitdir(filename)[2], control, data, runcommand, exec, run), structure
end

function write_block_data(f, data)
    if typeof(data) <: Vector{Vector{Float64}} || typeof(data) <: Vector{NTuple{4, Float64}} #k_points
        for x in data
            for y in x
                write(f, " $y")
            end
            write(f, "\n")
        end
    elseif typeof(data) <: Vector{Int} || typeof(data) <: NTuple{6, Int}
        for x in data
            write(f, " $x")
        end
        write(f, "\n")
    elseif typeof(data) <: Matrix
        im, jm = size(data)
        for i = 1:im
            for j = 1:jm
                write(f, " $(data[i, j])")
            end
            write(f, "\n")
        end
    end
end

"""
    write_input(input::QEInput, structure, filename::String=input.filename)

Writes a Quantum Espresso input file.
"""
function write_input(input::QEInput, structure, filename::String=input.filename)
    open(filename, "w") do f
        write_flag(flag_data) = write_flag_line(f, flag_data[1], flag_data[2])
        write_block(data)     = write_block_data(f, data)

        blocks2file = QEControlBlock[]
        c_block = getfirst(x -> x.name == :control, input.control)
        s_block = getfirst(x -> x.name == :system, input.control)
        if c_block != nothing
            push!(blocks2file, c_block)
            push!(blocks2file, s_block)
        end
        rest = filter(x -> x.name != :control && x.name != :system, input.control)
        for r in rest
            push!(blocks2file, r)
        end
        for block in blocks2file
            write(f, "&$(block.name)\n")
            if block.name == :system
                nat   = length(structure.atoms)
                ntyp  = length(unique_atoms(structure.atoms))
                # A     = 1.0
                ibrav = 0
                write(f,"  ibrav = $ibrav\n")
                # write(f,"  A = $A\n")
                write(f,"  nat = $nat\n")
                write(f,"  ntyp = $ntyp\n")
            end
            map(write_flag, [(flag, data) for (flag, data) in block.flags])
            write(f, "/\n\n")
        end

        write_structure(f, structure, input)
        for block in input.data
            if block.option != :none
                write(f, "$(uppercase(String(block.name))) ($(block.option))\n")
            else
                write(f, "$(uppercase(String(block.name)))\n")
            end
            if block.name == :k_points && block.option != :automatic
                write(f, "$(length(block.data))\n")
                write_block(block.data)
            else
                write_block(block.data)
            end
            write(f, "\n")
            #write the o_ther cards depending on whether or not they are there
        end
    end
end

function write_structure(f, structure, input::QEInput)
    unique_atoms = Symbol[]
    pseudo_lines = String[]
    atom_lines   = String[]
    for at in structure.atoms
        if !in(at.id, unique_atoms)
            push!(unique_atoms, at.id)
            push!(pseudo_lines, "$(at.id) $(at.element.atomic_weight)   $(at.pseudo)\n")
        end
        push!(atom_lines, "$(at.id)  $(at.position[1]) $(at.position[2]) $(at.position[3])\n")
    end

    write(f, "ATOMIC_SPECIES\n")
    write.(f, pseudo_lines)

    write(f, "\n")
    write(f, "CELL_PARAMETERS (angstrom)\n")
    write_cell(f, structure.cell)
    write(f, "\n")

    write(f, "ATOMIC_POSITIONS (angstrom) \n")
    write.(f, atom_lines)
    write(f, "\n")
end
