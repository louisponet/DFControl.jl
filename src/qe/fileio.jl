const parseable_qe_execs = ["pw.x", "projwfc.x", "pw2wannier90.x", "pp.x"]
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
                out[:k_cryst] = Vector{Vec3{T}}()
                line = readline(f)
                while line != "" && !contains(line, "--------")
                    push!(out[:k_cryst], parse_k_line(line, T))
                    line = readline(f)
                end

                #k_cart
            elseif contains(line, "cart.") && length(split(line)) == 5
                out[:k_cart] = Vector{Vec3{T}}()
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
                        out[:alat]            = contains(line, "angstrom") ? :angstrom : parse(T, split(line)[end][1:end-1])
                        out[:cell_parameters] = reshape(T[parse.(T, split(readline(f))); parse.(T, split(readline(f))); parse.(T, split(readline(f)))], (3,3))
                    elseif contains(line, "ATOMIC_POSITIONS")
                        out[:pos_option]      = cardoption(line)
                        line  = readline(f)
                        atoms = []
                        while !contains(line, "End")
                            s_line = split(line)
                            key    = Symbol(s_line[1])
                            push!(atoms, key=>Point3{T}(parse.(T, s_line[2:end])...))
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
        primv = alat(control, true) * Mat3(eye(3))
    elseif option == :bohr
        primv = conversions[:bohr2ang] * Mat3(eye(3))
    else
        primv = Mat3(eye(3))
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
function read_qe_input(filename, T=Float64::Type; exec=Exec("pw.x"), runcommand=Exec(""), run=true, structure_name="NoName")
    data    = Vector{InputData}()
    flags   = Dict{Symbol, Any}()
    atom_block     = nothing
    cell_block     = nothing
    pseudo_block   = nothing
    open(filename) do f
        line = readline(f)
        while !eof(f)
            @label start_label
            if contains(line, "&")
                line = readline(f)
                while strip(line) != "/"
                    if contains(line, "!")
                        line = readline(f)
                        continue
                    end
                    split_line = filter(x -> x != "", strip.(split(line, ",")))
                    for s in split_line
                        key, val = String.(strip.(split(s, "=")))
                        qe_flag  = Symbol(replace(replace(key, "(", "_"), ")",""))
                        flag_type = flagtype(QE, exec, qe_flag)
                        if flag_type != Void
                            t_val = parse_flag_val(val, flag_type)
                            flags[qe_flag] = eltype(t_val) == flag_type || flag_type == String ? t_val : error("Couldn't parse the value of flag '$key' in file '$filename'!")
                        else
                            error("Error reading $filename: flag '$key' not found in QE flag Dictionary for input $(exec.exec)!")
                        end
                    end
                    line = readline(f)
                end
                @goto start_label

            elseif contains(line, "CELL_PARAMETERS") || contains(line, "cell_parameters")
                cell_unit    = cardoption(line)
                cell_        = Matrix{T}(3, 3)
                cell_[1, 1:3] = parse.(T, split(readline(f)))
                cell_[2, 1:3] = parse.(T, split(readline(f)))
                cell_[3, 1:3] = parse.(T, split(readline(f)))
                cell = Mat3(cell_)
                line = readline(f)
                cell_block = InputData(:cell_parameters, cell_unit, cell)
                @goto start_label

            elseif contains(line, "ATOMIC_SPECIES") || contains(line, "atomic_species")
                line    = readline(f)
                pseudos = Dict{Symbol,String}()
                while length(split(line)) == 3
                    pseudos[Symbol(split(line)[1])] = split(line)[end]
                    line = readline(f)
                end
                pseudo_block = InputData(:atomic_species, :none, pseudos)
                @goto start_label

            elseif contains(line, "ATOMIC_POSITIONS") || contains(line, "atomic_positions")
                option = cardoption(line)
                atoms  = Dict{Symbol, Vector{Point3{T}}}()
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
                atom_block = InputData(:atomic_positions, option, atoms)
                @goto start_label

            elseif contains(line, "K_POINTS") || contains(line, "k_points")
                k_option = cardoption(line)
                line     = readline(f)
                if k_option == :automatic
                    s_line = split(line)
                    k_data = parse.(Int, s_line)
                else
                    nks    = parse(Int, line)
                    k_data = Vector{Vector{T}}(nks)
                    for i = 1:nks
                        k_data[i] = parse.(T, split(readline(f)))
                    end
                end
                push!(data, InputData(:k_points, k_option, k_data))
                @goto start_label
            end
            line = readline(f)
        end
    end

    structure = extract_structure!(structure_name, flags, cell_block, atom_block, pseudo_block)
    pop!.(flags, [:ibrav, :nat, :ntyp, :A, :celldm_1, :celldm], nothing)
    return DFInput{QE}(splitdir(filename)[2], flags, data, [runcommand, exec], run), structure
end

"""
    save(input::DFInput{QE}, structure, filename::String=input.filename)

Writes a Quantum Espresso input file.
"""
function save(input::DFInput{QE}, structure, filename::String=input.filename)
    open(filename, "w") do f
        write_flag(flag_data) = write_flag_line(f, flag_data[1], flag_data[2])
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
            map(write_flag, [(flag, data) for (flag, data) in flags])
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
    write.(f, pseudo_lines)

    write(f, "\n")
    write(f, "CELL_PARAMETERS (angstrom)\n")
    write_cell(f, structure.cell)
    write(f, "\n")

    write(f, "ATOMIC_POSITIONS (angstrom) \n")
    write.(f, atom_lines)
    write(f, "\n")
end
