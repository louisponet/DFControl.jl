#THIS IS THE MOST HORRIBLE FUNCTION I HAVE EVER CREATED!!!
function extract_atoms(atoms_block::T, proj_block::T, cell) where T <: WannierDataBlock
    if atoms_block.name == :atoms_cart
        cell = Mat3(eye(3))
    end
    projections = proj_block.data
    atoms = atoms_block.data
    t_start = 1
    out_ats = Atom{Float64}[]
    if projections != nothing
        t_ats = Atom{Float64}[]
        for (proj_at, projs) in projections
            for proj in projs
                for (pos_at, pos) in atoms
                    if proj_at != pos_at
                        continue
                    end
                    for ps in pos
                        size = orbsize(proj)
                        push!(t_ats, Atom(pos_at, element(pos_at), cell' * ps, :projections => [Projection(Orbital(proj), t_start, size, t_start + size - 1)]))
                        t_start += size
                    end
                end
            end
        end
        for (pos_at, pos) in atoms
            for ps in pos
                same_ats = Atom{Float64}[]
                for at in t_ats
                    if at.position == cell' * ps
                        push!(same_ats, at)
                    end
                end
                if length(same_ats) > 1
                    for at in same_ats[2:end]
                        push!(same_ats[1].projections, at.projections[1])
                    end
                end
                push!(out_ats, same_ats[1])
            end
        end
    else
        for (pos_at, pos) in atoms
            for p in pos
                push!(out_ats, Atom(pos_at, element(pos_at), cell' * p, :projections => :random))
            end
        end
    end
    return out_ats
end

function extract_structure(cell_block::T, atoms_block::T, projections_block::T, structure_name="NoName") where T <: Union{WannierDataBlock, Void}
    if atoms_block == nothing || cell_block == nothing
        return nothing
    end
    if cell_block.option == :ang
        cell = cell_block.data
    elseif cell_block.option == :bohr
        cell = conversions[:bohr2ang] * cell_block.data
    end

    atoms = extract_atoms(atoms_block, projections_block, cell)
    return Structure(structure_name, Mat3(cell), atoms)
end

"""
    read_wannier_input(filename::String, T=Float64)

Reads a DFInput from a wannier90 input file.
"""
function read_wannier_input(filename::String, T=Float64; run_command="", run=true, preprocess=true, structure_name="NoName")
    flags       = Dict{Symbol,Any}()
    data_blocks = Array{WannierDataBlock,1}()
    atoms_block = nothing
    cell_block  = nothing
    proj_block  = nothing
    open(filename,"r") do f
        line = readline(f)
        while !eof(f)
            @label start_label

            if contains(line, "!") || line == "" || contains(lowercase(line), "end") || contains(line, "#")
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
                            proj_block = WannierDataBlock(:projections, :random, nothing)
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
                    proj_block = WannierDataBlock(:projections, :none, proj_dict)
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
                        line = readline(f)
                    else
                        option = :ang
                    end
                    cell_param = Matrix{T}(3, 3)
                    for i = 1:3
                        cell_param[i, :] = parse_line(T, line)
                        line = readline(f)
                    end
                    cell_block = WannierDataBlock(:unit_cell_cart, option, cell_param)
                    # line = readline(f)
                    @goto start_label

                elseif block_name == :atoms_frac || block_name == :atoms_cart
                    line   = readline(f)
                    atoms  = Dict{Symbol,Array{Point3D{T},1}}()
                    option = :ang
                    while !contains(lowercase(line), "end")
                        if contains(line, "!") || line == ""
                            line = readline(f)
                            continue
                        end
                        if length(split(line)) == 1
                            option = parse(line)
                            line = readline(f)
                            continue
                        end
                        split_line = strip_split(line)
                        atom       = Symbol(split_line[1])
                        position   = Point3D(parse_string_array(T, split_line[2:4]))
                        if !haskey(atoms,atom)
                            atoms[atom] = [position]
                        else
                            push!(atoms[atom], position)
                        end
                        line = readline(f)
                    end
                    atoms_block = WannierDataBlock(block_name, option, atoms)
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
                    split_line = strip_split(line, '=')
                    flag       = Symbol(split_line[1])
                    value      = split_line[2]
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
    structure = extract_structure(cell_block, atoms_block, proj_block, structure_name)
    return WannierInput(splitdir(filename)[2], structure, flags, data_blocks, run_command, run, preprocess)
end

"""
    write_wannier_input(filename::String, input::DFInput)"

Writes the input of a wannier90 input file.
"""
function write_wannier_input(input::WannierInput, filename::String=input.filename)
    open(filename, "w") do f
        for (flag, value) in input.flags
            write_flag_line(f, flag, value)
        end
        write(f, "\n")

        structure = input.structure
        if structure != nothing
            write(f,"begin unit_cell_cart\n")
            write_cell(f, structure.cell)
            write(f,"end unit_cell_cart\n")
            write(f, "\n")
        end
        write(f, "begin projections\n")
        for at in unique_atoms(structure.atoms)
            if at.projections == :random || !isdefined(at, :projections)
                write(f, "random\n")
                break
            end
            write(f, "$(at.id): $(at.projections[1].orb)")
            if length(at.projections) > 1
                for proj in at.projections[2:end]
                    write(f, ";$(proj.orb)")
                end
            end
            write(f, "\n")
        end
        write(f, "end projections\n")

        write(f, "\n")
        write(f, "begin atoms_cart\n")
        for at in structure.atoms
            write(f, "$(at.id)  $(at.position[1]) $(at.position[2]) $(at.position[3])\n")
        end
        write(f, "end atoms_cart\n")
        write(f, "\n")

        for block in input.data_blocks
            write(f, "begin $(block.name)\n")
            if block.name == :kpoint_path
                for i = 1:2:length(block.data)
                    letter1, k_points1 = block.data[i]
                    letter2, k_points2 = block.data[i+1]
                    write(f, "$letter1 $(k_points1[1]) $(k_points1[2]) $(k_points1[3]) $letter2 $(k_points2[1]) $(k_points2[2]) $(k_points2[3])\n")
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
