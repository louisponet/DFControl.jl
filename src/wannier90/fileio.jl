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
                if typeof(block.data) <: Dict
                    for (atom,symbols) in block.data
                        write(f, "$atom: $(symbols[1])")
                        for sym in symbols[2:end]
                            write(f, ";$sym")
                        end
                        write(f, "\n")
                    end
                elseif typeof(block.data) == String
                    write(f, block.data * "\n")
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