
function parse_k_line(line, T)
    splt = split(line)
    k1   = parse(T, splt[5])
    k2   = parse(T, splt[6])
    k3   = parse(T, splt[7][1:1:end-2])
    return [k1, k2, k3]
end

function write_flag_line(f, flag, data, seperator="=", i="")
    write(f,"  $flag$i $seperator")

    if typeof(data) <: Array

        if length(data) < 20
            write(f, "  $(data[1])")
            for x in data[2:end]
                write(f, " $x")
            end
            write(f, "\n")
        else
            write(f, "\n")
            for i = 1:3:length(data)
                write(f, "  $(data[i]) $(data[i + 1]) $(data[i + 2])\n")
            end
        end

    else #this should work for anything singular valued data such as bools, ''s and other types
        write(f, "$data\n")
    end

end

function parse_flag_val(val, T=Float64)
    if T == String
        return val
    end

    if contains(val, "d")
        val = replace(val, "d", "e")
    end

    val = strip(val, '.') 
    t = convert.(T, parse.(split(lowercase(val))))
    #deal with abinit constants -> all flags that are read which are not part of the abi[:structure] get cast into the correct atomic units!
    if length(t) > 1 && typeof(t[end]) == Symbol
        t = t[1:end-1] .* abi_const[t[end]]
    end

    length(t) == 1 ? t[1] : typeof(t) <: Array{Real,1} ? convert.(T,t) : t  
end


"Quantum espresso card option parser"
function get_card_option(line)
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

Reads a generic quantum espresso input, returning a dictionary with all found data in the file.
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
                out[:polarization] = Point3D{T}(P * s_line[1], P * s_line[2], P * s_line[3])
                out[:pol_mod]      = mod
                
                #fermi energy
            elseif contains(line, "Fermi")
                out[:fermi]        = parse(T, split(line)[5])
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
                        out[:pos_option]      = get_card_option(line)
                        line  = readline(f)
                        atoms = Dict{Symbol,Array{Point3D{T},1}}()
                        while !contains(line, "End")
                            s_line = split(line)
                            key    = Symbol(s_line[1])
                            if key in keys(atoms)
                                push!(atoms[key], Point3D{T}(parse.(T, s_line[2:end])...))
                            else
                                atoms[key] = [Point3D{T}(parse.(T, s_line[2:end])...)]
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
                push!(out[:bands], DFBand(out[:k_cart], out[:k_cryst], eig_band))
            end
        end
        return out
    end
end

"""
    read_qe_bands_file(filename::String, T=Float64)

Reads the output file of a 'bands' calculation in Quantum Espresso.
Returns an array of DFBands each with the same k_points and their respective energies.
"""
read_qe_bands_file(filename::String, T=Float64) = read_qe_output(filename, T)[:bands]

"""
    read_ks_from_qe_bands_file(filename::String, T=Float64)

Read k-points from a Quantum Espresso bands output file in cartesian (2pi/alat in Angstrom^-1!) and crystalline coordinates.
Returns (k_points_cart,k_points_cryst).
"""
function read_ks_from_qe_bands_file(filename::String, T=Float64)
    t = read_qe_output(filename, T)
    return t[:k_cart], t[:k_cryst]
end

"""
    read_fermi_from_qe_file(filename::String,T=Float64)

Reads the Fermi level from a Quantum Espresso scf calculation output file
(if there is one).
"""
read_fermi_from_qe_file(filename::String, T=Float64) = read_qe_output(filename, T)[:fermi]

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

function read_errors(filename::String)
    out = String[]
    open(filename, "r") do f
        line = readline(f)
        while !eof(f)
            if contains(lowercase(line), "error")
                push!(out, line)
            end
        end
        return out
    end
end

#Incomplete. certain variables are not read well 

#I should use the extra info I got out of the QEControlFlags to parse the correct things!
"""
    read_qe_input(filename, T=Float64)

Reads a Quantum Espresso input file.
Returns a DFInput.
"""
function read_qe_input(filename, T=Float64::Type; run_command="", run=true)
    control_blocks = Array{QEControlBlock,1}()
    data_blocks    = Array{QEDataBlock,1}()
    
    flags_to_discard = ["nat", "ntyp"]
    open(filename) do f
        line = readline(f)
        while !eof(f)
            @label start_label
            if contains(line, "&")
                c_block_name    = Symbol(lowercase(strip(strip(line), '&')))
                flag_dict       = Dict{Symbol,Any}()
                def_block_flags = get_qe_flags(c_block_name)
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
                        if key in flags_to_discard
                            continue
                        elseif haskey(def_block_flags, qe_flag)
                            t_val = parse_flag_val(val, def_block_flags[qe_flag])
                            flag_dict[Symbol(key)] = eltype(t_val) == def_block_flags[qe_flag] || def_block_flags[qe_flag] == String ? t_val : error("Couldn't parse the value of flag '$key' in file '$filename'!") 
                        else
                            error("Error reading $filename: flag '$key' not found in QE flag dictionary for control block $c_block_name !")
                        end
                    end
                    line = readline(f)
                end
                push!(control_blocks, QEControlBlock(c_block_name, flag_dict))
                @goto start_label
                
            elseif contains(line, "CELL_PARAMETERS") || contains(line, "cell_parameters")
                cell_unit    = get_card_option(line)
                cell         = Matrix{T}(3, 3)
                cell[1, 1:3] = parse.(T, split(readline(f)))
                cell[2, 1:3] = parse.(T, split(readline(f)))
                cell[3, 1:3] = parse.(T, split(readline(f)))
                line = readline(f)
                push!(data_blocks, QEDataBlock(:cell_parameters, cell_unit, cell))
                @goto start_label
                
            elseif contains(line, "ATOMIC_SPECIES") || contains(line, "atomic_species")
                line    = readline(f)
                pseudos = Dict{Symbol,String}()
                while length(split(line)) == 3
                    pseudos[Symbol(split(line)[1])] = split(line)[end]
                    line = readline(f)
                end
                push!(data_blocks, QEDataBlock(:atomic_species, :none, pseudos))
                @goto start_label
                
            elseif contains(line, "ATOMIC_POSITIONS") || contains(line, "atomic_positions")
                option = get_card_option(line)
                atoms  = Dict{Symbol,Array{Point3D{T},1}}()
                line   = readline(f)
                while length(split(line)) == 4
                    s_line   = split(line)
                    atom     = Symbol(s_line[1])
                    position = Point3D(parse(T, s_line[2]), parse(T, s_line[3]), parse(T, s_line[4]))
                    if !haskey(atoms, atom)
                        atoms[atom] = [position]
                    else
                        push!(atoms[atom], position)
                    end
                    line = readline(f)
                end
                push!(data_blocks, QEDataBlock(:atomic_positions, option, atoms))
                @goto start_label
                
            elseif contains(line, "K_POINTS") || contains(line, "k_points")
                k_option = get_card_option(line)
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
                push!(data_blocks, QEDataBlock(:k_points, k_option, k_data))
                @goto start_label
            end
            line = readline(f)
        end
    end
    return QEInput(splitdir(filename)[2], control_blocks, data_blocks, run_command, run)
end

#can I use @generated here?
function write_block_data(f, data)
    if typeof(data) <: Array{Vector{Float64},1} || typeof(data) <: Array{Vector{Float64},1} #k_points
        for x in data
            for y in x
                write(f, " $y")
            end
            write(f, "\n")
        end
    elseif typeof(data) <: Array{Int,1}
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
    elseif typeof(data) <: Dict{Symbol,<:Any}
        for (key, value) in data
            if typeof(value) == String
                if length(String(key)) > 2 
                    write(f, "$key $(Float64(ELEMENTS[Symbol(String(key)[1:end - 1])].atomic_weight))   $value\n")
                else
                    write(f, "$key $(Float64(ELEMENTS[key].atomic_weight))   $value\n")
                end
            elseif typeof(value) <: Array{<:Point3D,1}
                for at in value
                    write(f, "$key $(at.x) $(at.y) $(at.z)\n")
                end
            end
        end
    end
end

"""
    write_qe_input(filename::String, df_input::DFInput)

Writes a Quantum Espresso input file.
"""
function write_qe_input(input::QEInput, filename::String=input.filename)
    open(filename, "w") do f
        write_flag(flag_data) = write_flag_line(f, flag_data[1], flag_data[2])
        write_block(data)     = write_block_data(f, data)
        
        for block in input.control_blocks
            write(f, "&$(block.name)\n")
            if block.name == :system
                atoms = get_data(input, :atomic_positions)
                nat   = sum(length.([ps for ps in values(atoms)]))
                ntyp  = length(atoms)
                write_flag((:nat,  nat))
                write_flag((:ntyp, ntyp))
            end
            map(write_flag, [(flag, data) for (flag, data) in block.flags])
            write(f, "/\n\n")
        end
        
        for block in input.data_blocks
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


#---------------------------END QUANTUM ESPRESSO SECTION----------------#
#---------------------------START WANNIER SECTION ----------------------#

"""
    read_wannier_input(filename::String, T=Float64)

Reads a DFInput from a wannier90 input file.
"""
function read_wannier_input(filename::String, T=Float64; run_command="", run=true, preprocess=true)
    flags       = Dict{Symbol,Any}()
    data_blocks = Array{WannierDataBlock,1}()
    open(filename,"r") do f
        line = readline(f)
        while !eof(f)
            @label start_label
            
            if contains(line, "!") || line == "" || contains(lowercase(line), "end")
                line = readline(f)
                continue
            end

            if contains(lowercase(line), "begin")
                block_name = Symbol(split(lowercase(line))[end])
                
                if block_name == :projections
                    proj_dict = Dict{Symbol,Array{Symbol,1}}()
                    line      = readline(f)
                    while !contains(lowercase(line), "end")
                        if contains(line, "!") || line == ""
                            line = readline(f)
                            continue
                        end
                        if contains(line, "random")
                            push!(data_blocks, WannierDataBlock(:projections, :random, nothing))
                            line = readline(f)
                            break
                        else
                            split_line      = strip_split(line, ':')
                            atom            = Symbol(split_line[1])
                            projections     = [Symbol(proj) for proj in strip_split(split_line[2], ';')]
                            proj_dict[atom] = projections
                            line = readline(f)
                        end
                    end
                    push!(data_blocks, WannierDataBlock(:projections, :none, proj_dict))
                    @goto start_label
                
                elseif block_name == :kpoint_path
                    line = readline(f)
                    k_path_array = Array{Tuple{Symbol,Array{T,1}},1}()
                    while !contains(lowercase(line), "end")
                    if contains(line, "!") || line == ""
                        line = readline(f)
                        continue
                    end
                    split_line = split(line)
                    push!(k_path_array, (Symbol(split_line[1]), parse_string_array(T, split_line[2:4])))
                    push!(k_path_array, (Symbol(split_line[5]), parse_string_array(T, split_line[6:8])))
                    line = readline(f)
                end
                push!(data_blocks, WannierDataBlock(:kpoint_path, :none, k_path_array))
                @goto start_label
            
                elseif block_name == :unit_cell_cart
                    line = readline(f)
                    if length(split(line)) == 1
                        option = Symbol(lowercase(line))
                    else
                        option = :ang
                    end
                    cell_param = Matrix{T}(3, 3)
                    for i = 1:3
                        cell_param[i, :] = parse_line(T, readline(f))
                    end
                    push!(data_blocks, WannierDataBlock(:unit_cell_cart, option, cell_param))
                    line = readline(f)
                    @goto start_label
            
                elseif block_name == :atoms_frac || block_name == :atoms_cart
                    line   = readline(f)
                    atoms  = Dict{Symbol,Array{Point3D{T},1}}()
                    option = Symbol(split(String(block_name), "_")[end])
                    while !contains(lowercase(line), "end")
                        if contains(line, "!") || line == ""
                            line = readline(f)
                            continue
                        end
                        atom       = Symbol(split_line[1])
                        position   = Point3D(parse_string_array(T, split_line[2:4]))
                        split_line = strip_split(line)
                        if !haskey(atoms,atom)
                            atoms[atom] = [position]
                        else
                            push!(atoms[atom], position)
                        end
                        line = readline(f)
                    end
                    push!(data_blocks, WannierDataBlock(block_name, option, atoms))
                    @goto start_label
        
                elseif block_name == :kpoints
                    line     = readline(f)
                    k_points = Array{Array{T,1},1}()
                    while !contains(lowercase(line), "end")
                        if line == ""
                            line = readline(f)
                            continue
                        end
                        push!(k_points, parse_line(T, line))
                        line = readline(f)
                    end
                    push!(data_blocks, WannierDataBlock(:kpoints, :none, k_points))
                    @goto start_label
                end

            else
                if contains(line, "mp_grid")
                    flags[:mp_grid] = parse_string_array(Int, split(split(line, '=')[2]))
                else
                    flag       = Symbol(split_line[1])
                    value      = split_line[2]
                    split_line = strip_split(line, '=')
                    if  lowercase(value) == "t" || lowercase(value) == "true"
                        flags[flag] = true
                    elseif lowercase(value) == "f" || lowercase(value) == "false"
                        flags[flag] = false
                    elseif !isnull(tryparse(Int, value))
                        flags[flag] = get(tryparse(Int, value))
                    elseif !isnull(tryparse(T, value))
                        flags[flag] = get(tryparse(T, value))
                    else
                        flags[flag] = value
                    end
                end
            end
            line = readline(f)
        end
    end
    return WannierInput(splitdir(filename)[2], flags, data_blocks, run_command, run, preprocess)
end

"""
    write_wannier_input(filename::String, input::DFInput)

Writes the input of a wannier90 input file.
"""
function write_wannier_input(input::WannierInput, filename::String=input.filename)
    open(filename, "w") do f
        for (flag, value) in input.flags
            write_flag_line(f, flag, value)
        end
        write(f, "\n")
        for block in input.data_blocks
            write(f, "begin $(block.name)\n")
            if block.name == :kpoint_path
                for i = 1:2:length(block.data)
                    letter1, k_points1 = block.data[i]
                    letter2, k_points2 = block.data[i+1]
                    write(f, "$letter1 $(k_points1[1]) $(k_points1[2]) $(k_points1[3]) $letter2 $(k_points2[1]) $(k_points2[2]) $(k_points2[3])\n")
                end
                
            elseif block.name == :projections
                for (atom,symbols) in block.data
                    write(f, "$atom: $(symbols[1])")
                    for sym in symbols[2:end]
                        write(f, ";$sym")
                    end
                    write(f, "\n")
                end
                write(f, "\n")
                
            elseif block.name == :unit_cell_cart
                matrix = block.data
                write(f, "$(block.option)\n")
                write(f, "$(matrix[1, 1]) $(matrix[1, 2]) $(matrix[1, 3])\n")
                write(f, "$(matrix[2, 1]) $(matrix[2, 2]) $(matrix[2, 3])\n")
                write(f, "$(matrix[3, 1]) $(matrix[3, 2]) $(matrix[3, 3])\n")
                
            elseif block.name == :atoms_frac || block.name == :atoms_cart
                for (atom, points) in block.data
                    for point in points
                        write(f, "$atom $(point.x) $(point.y) $(point.z)\n")
                    end
                end
                
            elseif block.name == :kpoints
                for k in block.data
                    write(f, "$(k[1]) $(k[2]) $(k[3])\n")
                end
            end
            write(f, "end $(block.name)\n\n")
        end
    end
end
#---------------------------END WANNIER SECTION ------------------------#

#---------------------------START ABINIT SECTION------------------------#
@pyimport abipy.abio.abivars as abivars
@pyimport abipy.abilab as abilab
#also searches the structure but where not doing anything with it now
"""
    read_abi_input(filename::String, T=Float64)

Returns an ABINIT input.
"""
function read_abi_input(filename::String, T=Float64; run_command="", run=true, pseudos=[""])
    abi_input = PyObject(abivars.AbinitInputFile(filename))
    datasets  = abi_input[:datasets]
    structure = abi_input[:structure]
    cell      = T.(structure[:lattice_vectors]())
    atoms     = Dict{Symbol,Array{Point3D{T},1}}()
    for site in structure[:sites]
        atsym   = Symbol(site[:specie][:symbol])
        coord   = Point3D{T}(T.(site[:frac_coords]))
        if !haskey(atoms, atsym)
            atoms[atsym] = [coord]
        else
            push!(atoms[atsym], coord)
        end
    end
    inputs           = Array{AbinitInput,1}()
    filename         = splitdir(filename)[2]
    flags_to_discard = ["acell","rprim","natom","ntypat","typat","znucl","xred","acell"]
    for (i, data) in enumerate(datasets)
        file, ext= splitext(filename)
        if length(datasets) > 1
            newfile = file * "$i" * ext
        else
            newfile = file * ext
        end
        flags = Dict{Symbol,Any}()
        for (flag, value) in data
            if flag in flags_to_discard
                continue
            else
                flag_type           = get_abi_flag_type(Symbol(flag))
                flags[Symbol(flag)] = flag_type != Void ? parse_flag_val(value,flag_type) : error("Couldn't parse flag '$flag' with value '$value'!")
            end
        end
        push!(inputs, AbinitInput(newfile, flags,
                                  [AbinitDataBlock(:cell_parameters,  :angstrom, cell),
                                   AbinitDataBlock(:atomic_positions, :frac, atoms),
                                   AbinitDataBlock(:pseudos, :pseudos, pseudos)],
                                   structure,
                                   run_command,
                                   run))
    end
    return inputs
end

function write_abi_input(input::AbinitInput, filename::String=input.filename)
    flags = input.flags
    open(filename, "w") do f
        write_flag(flag_data) = write_flag_line(f, flag_data[1], flag_data[2], "")
        write(f, input.structure[:abi_string])
        write(f, "\n")
        write_flag.(collect(flags))
    end
end

#question: Does it matter that we might be writing redundant data such as spinat etc?
"""
    write_abi_datasets(inputs::Array{AbinitInput,1}, directory)

Takes all the inputs, sees which ones have the same structure and constructs input files for each seperate sturcture.
The filename of the first input of a certain structure is used as file for all the datasets.
"""
function write_abi_datasets(inputs::Array{AbinitInput,1}, directory)
    input_groups = Array{Array{AbinitInput,1},1}([[pop!(inputs)]])
    while length(inputs) != 0
        input = pop!(inputs)
        for group in input_groups
            if group[1].structure == input.structure
                push!(group, input)
                break
            end
        end
    end
    
    for group in input_groups
        run_indices = Int[]
        for (i, _input) in enumerate(group)
            if _input.run 
                push!(run_indices, i)
            end
        end
        if length(run_indices) == length(group)
            jdtset = ""
        else
            jdtset = join(["$c" for c in run_indices], " ")
        end
        open(directory * group[1].filename,"w") do f
            write(f, "ndtset $(length(group)) jdtset $jdtset\n")
            write(f, group[1].structure[:abi_string])
            write(f, "\n")
            for t = 1:length(group)
                write_flag(flag_data) = write_flag_line(f, flag_data[1], flag_data[2], "", t)
                write(f, "#=============== BEGIN DATASET $t ===============#\n")
                write_flag.(collect(group[t].flags))
                write(f, "#================ END DATASET $t ================#\n")
            end
        end
    end
    return map(x -> (x[1].filename, get_data(x[1], :pseudos), x[1].run_command), input_groups)
end

function read_abi_output(filename::String)
    result = abilab.abiopen(filename)
    return result
end


#---------------------------END ABINIT SECTION--------------------------#

#---------------------------BEGINNING GENERAL SECTION-------------------# 
"""
    write_input(df_input::DFInput, filename::String=df_input.filename)

Writes the input file for a DFInput.Backend of DFInput decides what writing function is called.
"""
function write_input(df_input::DFInput, filename::String=df_input.filename)
    if typeof(df_input) == QEInput
        write_qe_input(df_input, filename)
    elseif typeof(df_input) == WannierInput
        write_wannier_input(df_input, filename)
    elseif typeof(df_input) == AbinitInput
        write_abi_input(df_input, filename)
    end
end

#Incomplete: only works with SBATCH right now
function write_job_name(job::DFJob, f)
    write(f, "#SBATCH -J $(job.name) \n")
end

function write_job_header(job::DFJob, f)
    job_header = job.header == "" ? get_default_job_header() : job.header
    for line in job_header
        if contains(line, "\n")
            write(f, line)
        else
            write(f, line * "\n")
        end
    end
end

"""
    write_job_files(job::DFJob)

Writes all the input files and job file that are linked to a DFJob.
"""
function write_job_files(job::DFJob)
    files_to_remove = search_dir(job.local_dir, ".in")
    new_filenames   = String[]
    num_abi         = 0 
    open(job.local_dir * "job.tt", "w") do f
        write(f, "#!/bin/bash\n")
        write_job_name(job, f)
        write_job_header(job, f)
        i = 1
        while i < length(job.calculations)
            calculation = job.calculations[i]
            run_command = calculation.run_command
            filename    = calculation.filename
            should_run  = calculation.run
            if typeof(calculation) == WannierInput
                write_input(calculation, job.local_dir * filename)
                if calculation.preprocess
                    run_command *= " -pp"
                end
                if !should_run
                    write(f, "#$run_command $(filename[1:end-4])\n")
                else
                    write(f, "$run_command $(filename[1:end-4])\n")
                end
                i += 1
            elseif typeof(calculation) == AbinitInput
                abinit_inputs     = Array{AbinitInput,1}(filter(x -> typeof(x) == AbinitInput, job.calculations))
                i += length(abinit_inputs)
                abinit_jobfiles   = write_abi_datasets(abinit_inputs, job.local_dir)
                for (filename, pseudos, run_command) in abinit_jobfiles
                    file, ext = splitext(filename)
                    write(f, "$run_command << !EOF\n$filename\n$(file * ".out")\n$(job.name * "_Xi$num_abi")\n$(job.name * "_Xo$num_abi")\n$(job.name * "_Xx$num_abi")\n")
                    for pp in pseudos
                        write(f, "$pp\n")
                    end
                    write(f, "!EOF\n")
                    num_abi += 1 
                end
                
            elseif typeof(calculation) == QEInput
                write_input(calculation, job.local_dir * filename)
                if !should_run
                    write(f, "#$run_command <$filename> $(split(filename,".")[1]).out \n")
                else
                    write(f, "$run_command <$filename> $(split(filename,".")[1]).out \n")
                end
                i += 1 
            end
            push!(new_filenames, filename)
        end
    end
    
    for file in files_to_remove
        if !in(file, new_filenames)
            rm(job.local_dir * file)
        end
    end
end


function read_command_line(line)
    if typeof(line) <: String
        line = split(line)
    end
    i = 0
    for (j, s) in enumerate(line)
        if contains(s, ".x")
            i = j
            break
        elseif contains(s, "abinit")
            i=j
            break
        end
    end
    run_command = prod([s * " " for s in line[1:i]])
    if contains(line[i + 1],"-")
        run_command *= line[i + 1]
        i+=1
    end
    return i, run_command
end


"""
    read_job_file(job_file::String)

Reads and returns the name, input files, run_commands and whether or not they need to be commented out.
All files that are read contain "in".
This reads QE and wannier90 inputs for now.
"""
function read_job_file(job_file::String)
    data = Dict{Symbol,Any}()
    data[:name]         = ""
    data[:header]       = Array{String,1}()
    data[:input_files]  = Array{String,1}() 
    data[:output_files] = Array{String,1}() 
    data[:run_commands] = Array{String,1}() 
    data[:should_run]   = Array{Bool,1}()
    open(job_file, "r") do f
        readline(f)
        while !eof(f)
            line = readline(f)
            if line == ""
                continue
            end
            if contains(line, ".x ")
                if contains(line, "#")
                    push!(data[:should_run], false)
                    line = line[2:end]
                else
                    push!(data[:should_run], true)
                end
                
                s_line        = split(line)
                i,run_command = read_command_line(s_line)      
                push!(data[:run_commands], run_command)
                #handles QE and Wannier.
                in_out = strip_split(prod(s_line[i + 1:end]), ">")
                if length(in_out) == 2
                    push!(data[:input_files],  strip(in_out[1], '<'))
                    push!(data[:output_files], in_out[2])
                else
                    push!(data[:input_files], strip(in_out[1], '<'))
                end 
            elseif contains(line, "abinit ")
                data[:abinit_pseudos] = Array{String,1}()
                s_line         = split(line)
                i, run_command = read_command_line(s_line)
                push!(data[:run_commands], strip(run_command, '#'))
                if contains(line, "!EOF")
                    push!(data[:input_files],  strip(readline(f), '#'))
                    push!(data[:output_files], strip(readline(f), '#'))
                    if contains(line, "#")
                        push!(data[:should_run], false)
                    else
                        push!(data[:should_run], true)
                    end
                    line = readline(f)
                    while !contains(line, "EOF")
                        if contains(line, ".xml")
                            push!(data[:abinit_pseudos], strip(line, '#'))
                        end
                        line = readline(f)
                    end
                end
                
                #this is reading the sbatch lines
            elseif contains(line, "#SBATCH")
                if contains(line, "-J")
                    data[:name] = split(line)[end]
                else
                    push!(data[:header], line)
                end
            else  
                push!(data[:header], line)
            end
        end
    end
    return data
end

#Incomplete: because QE is stupid we have to find a way to distinguish nscf and bands outputs hardcoded.
# function read_output(filename::string, args...)
#   open(filename,"r") do f
#     while !eof(f)
#       line = readline(f)

#       if contains(line,"self-consistent calculation")
#         return
#       elseif contains(line,"band structure calculation") && !contains(filename,"nscf")
#         return read_qe_bands_file(filename,args...)
#       elseif contains(line, "


#---------------------------END GENERAL SECTION-------------------#

#Incomplete: we should probably handle writing an array of expressions as well!
function expr2file(filename::String, expression::Expr)
    eq        = Symbol("=")
    lines     = readlines(filename)
    new_lines = String[]
    found     = false
    
    if expression.head != eq
        error("For now only writing of assignment expressions is possible.")
    end
    
    lhs = expression.args[1]
    rhs = expression.args[2]
    
    for line in lines
        if line == ""
            continue
        end
        
        expr = parse(line)
        if typeof(expr) == Void
            continue
        end
        if expr.head != eq
            continue
        end
        
        lhs_t = expr.args[1]
        rhs_t = expr.args[2]
        
        if lhs_t == lhs
            found = true
            push!(new_lines, "$(:($lhs = $rhs))")
        else
            push!(new_lines, "$expr")
        end
    end
    
    open(filename, "w") do f
        for line in new_lines
            write(f, line * "\n")
        end
        if !found
            write(f, "$expression\n")
        end
    end
end

function rm_expr_lhs(filename, lhs)
    lines       = readlines(filename)
    write_lines = String[]
    ind_2_rm    = 0
    
    for line in lines
        lhs_t = parse(line).args[1]
        if lhs_t == lhs
            continue
        else
            push!(write_lines, line)
        end
    end
    
    open(filename, "w") do f
        for line in write_lines
            write(f,line * "\n")
        end
    end
end

