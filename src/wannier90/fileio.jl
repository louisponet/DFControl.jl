#THIS IS THE MOST HORRIBLE FUNCTION I HAVE EVER CREATED!!!
#extracts only atoms with projections
function extract_atoms(atoms_block::T, proj_block::T, cell, spinors=false) where T <: InputData
    if atoms_block.name == :atoms_cart
        cell = Mat3(Matrix(I, 3, 3))
    end
    projections_ = proj_block.data
    atoms = atoms_block.data
    t_start = 1
    out_ats = Atom{Float64}[]
    if projections_ != nothing
        t_ats = Atom{Float64}[]
        for (pos_at, pos) in atoms
            for ps in pos
                for (proj_at, projs) in projections_
                    if proj_at != pos_at
                        continue
                    end
                    for proj in projs
                        size = spinors ? 2*orbsize(proj) : orbsize(proj)
                        push!(t_ats, Atom(pos_at, element(pos_at), cell' * ps, projections=[Projection(orbital(proj), t_start, t_start + size - 1)]))
                        t_start += size
                    end
                end
            end
        end
        for (pos_at, pos) in atoms
            for ps in pos
                same_ats = Atom{Float64}[]
                for at in t_ats
                    if position(at) == cell' * ps
                        push!(same_ats, at)
                    end
                end
                if isempty(same_ats)
                    continue
                end
                if length(same_ats) > 1
                    for at in same_ats[2:end]
                        push!(projections(same_ats[1]), projections(at)[1])
                    end
                end
                push!(out_ats, same_ats[1])
            end
        end
    else
        for (pos_at, pos) in atoms
            for p in pos
                push!(out_ats, Atom(pos_at, element(pos_at), cell' * p, projections=:random))
            end
        end
    end
    return out_ats
end

function extract_structure(name, cell_block::T, atoms_block::T, projections_block::T, spinors=false) where T <: Union{InputData, Nothing}
    if atoms_block == nothing || cell_block == nothing
        return nothing
    end
    if cell_block.option == :bohr
        cell = conversions[:bohr2ang] * cell_block.data
    else
        cell = cell_block.data
    end

    atoms = extract_atoms(atoms_block, projections_block, cell, spinors)
    return Structure(name, Mat3(cell), atoms)
end

function wan_read_input(::Type{T}, f::IO) where T
    flags = Dict{Symbol,Any}()
    data  = Vector{InputData}()
    atoms_block = nothing
    cell_block  = nothing
    proj_block  = nothing
    line = readline(f)
    while !eof(f)
        @label start_label

        if occursin("!", line) || line == "" || occursin("end", lowercase(line)) || occursin("#", line)
            line = readline(f)
            continue
        end

        if occursin("begin", lowercase(line))
            block_name = Symbol(split(lowercase(line))[end])

            if block_name == :projections
                proj_dict = Dict{Symbol,Array{Symbol,1}}()
                line      = readline(f)
                while !occursin("end", lowercase(line))
                    if occursin("!", line) || line == ""
                        line = readline(f)
                        continue
                    end
                    if occursin("random", line)
                        proj_block = InputData(:projections, :random, nothing)
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
                proj_block = InputData(:projections, :none, proj_dict)
                @goto start_label

            elseif block_name == :kpoint_path
                line = readline(f)
                k_path_array = Array{Tuple{Symbol,Array{T,1}},1}()
                while !occursin("end", lowercase(line))
                if occursin("!", line) || line == ""
                    line = readline(f)
                    continue
                end
                split_line = split(line)
                push!(k_path_array, (Symbol(split_line[1]), parse_string_array(T, split_line[2:4])))
                push!(k_path_array, (Symbol(split_line[5]), parse_string_array(T, split_line[6:8])))
                line = readline(f)
            end
            push!(data, InputData(:kpoint_path, :none, k_path_array))
            @goto start_label

            elseif block_name == :unit_cell_cart
                line = readline(f)
                if length(split(line)) == 1
                    option = Symbol(strip(lowercase(line)))
                    line = readline(f)
                else
                    option = :ang
                end
                cell_param = Matrix{T}(undef, 3, 3)
                for i = 1:3
                    cell_param[i, :] = parse_line(T, line)
                    line = readline(f)
                end
                cell_block = InputData(:unit_cell_cart, option, Mat3(cell_param))
                # line = readline(f)
                @goto start_label

            elseif block_name == :atoms_frac || block_name == :atoms_cart
                line   = readline(f)
                atoms  = Dict{Symbol, Array{Point3{T}, 1}}()
                option = :ang
                while !occursin("end", lowercase(line))
                    if occursin("!", line) || line == ""
                        line = readline(f)
                        continue
                    end
                    if length(split(line)) == 1
                        option = Symbol(strip(line))
                        line = readline(f)
                        continue
                    end
                    split_line = strip_split(line)
                    atom       = Symbol(split_line[1])
                    position   = Point3(parse_string_array(T, split_line[2:4]))
                    if !haskey(atoms,atom)
                        atoms[atom] = [position]
                    else
                        push!(atoms[atom], position)
                    end
                    line = readline(f)
                end
                atoms_block = InputData(block_name, option, atoms)
                @goto start_label

            elseif block_name == :kpoints
                line     = readline(f)
                k_points = Array{Array{T,1},1}()
                while !occursin("end", lowercase(line))
                    if line == ""
                        line = readline(f)
                        continue
                    end
                    push!(k_points, parse_line(T, line))
                    line = readline(f)
                end
                push!(data, InputData(:kpoints, :none, k_points))
                @goto start_label
            end

        else
            if occursin("mp_grid", line)
                flags[:mp_grid] = parse_string_array(Int, split(line)[end-2:end])
            else
                flag, val = wan_parse_flag_line(line)
                flags[flag] = val 
            end
        end
        line = readline(f)
    end
    return flags, data, atoms_block, cell_block, proj_block
end

wan_read_input(f::IO) = wan_read_input(Float64, f)

"""
    wan_read_input(filename::String, T=Float64; runcommand= Exec(""), run=true, exec=Exec("wannier90.x"), structure_name="NoName")

Reads a `DFInput{Wannier90}` and the included `Structure` from a WANNIER90 input file.
"""
function wan_read_input(filename::String, T=Float64; execs=[Exec("wannier90.x")], run=true, structure_name="NoName")
    flags = Dict{Symbol,Any}()
    data  = Vector{InputData}()
    atoms_block = nothing
    cell_block  = nothing
    proj_block  = nothing
    open(filename,"r") do f
    	flags, data, atoms_block, cell_block, proj_block = wan_read_input(T, f)
	end
    structure = extract_structure(structure_name, cell_block, atoms_block, proj_block, get(flags, :spinors, false))
    dir, file = splitdir(filename)
    flags[:preprocess] = hasflag(getfirst(x->x.exec == "wannier90.x", execs), :pp)  ? true : false
    return DFInput{Wannier90}(splitext(file)[1], dir, flags, data, execs, run), structure
end

function wan_parse_flag_line(line::String)
    split_line = strip_split(line, '=')
    flag       = Symbol(split_line[1])
    flagtyp   = flagtype(Wannier90, flag)
    value      = strip(lowercase(split_line[2]), '.')
    if flagtyp != String
	    value = replace(value, "d" => "e")
    end
    spl_val    = split(value)
    if eltype(flagtyp) == Bool
        for i = 1:length(spl_val)
            if spl_val[i] == "t"
                spl_val[i] = "true"
            elseif spl_val[i] == "f"
                spl_val[i] = "false"
            end
        end
    end
    parsed_val = flagtyp == String ? spl_val : parse.(eltype(flagtyp), spl_val)
    if length(parsed_val) == 1
        val = parsed_val[1]
    else
        val = parsed_val
    end
    return flag, val
end

function projections_string(at::AbstractAtom)
	prjs = projections(at)
	str = "$(name(at)): $(prjs[1].orb.name)"
    if length(prjs) > 1
        for proj in prjs[2:end]
            str *= ";$(proj.orb.name)"
        end
    end
    return str
end

function wan_write_projections(f::IO, atoms::Vector{<:AbstractAtom})
    write(f, "begin projections\n")
    uniats = unique(atoms)
    projs = projections.(uniats)
    # projs = projections.(unique(atoms(structure)))
    if all(isempty.(projs))
        write(f, "random\n")
    else
        for (at, prjs) in zip(uniats, projs)
            if isempty(prjs)
                continue
            end
            write(f, projections_string(at))
            write(f, "\n")
        end
    end
    write(f, "end projections\n")
end

"""
    save(input::DFInput{Wannier90}, structure, filename::String=inpath(input))

Writes the `DFInput{Wannier90}` and `structure` to a file, that can be interpreted by WANNIER90.
The atoms in the structure must have projections defined.
"""
function save(input::DFInput{Wannier90}, structure, filename::String=inpath(input))
    open(filename, "w") do f
        for (flag, value) in input.flags
            write_flag_line(f, flag, value)
        end
        write(f, "\n")

        if structure != nothing
            write(f,"begin unit_cell_cart\n")
            write_cell(f, cell(structure))
            write(f,"end unit_cell_cart\n")
            write(f, "\n")
        end
        wan_write_projections(f, atoms(structure))

        write(f, "\n")
        write(f, "begin atoms_cart\n")
        for at in atoms(structure)
            pos = position(at)
            write(f, "$(name(at))  $(pos[1]) $(pos[2]) $(pos[3])\n")
        end
        write(f, "end atoms_cart\n")
        write(f, "\n")

        for block in input.data
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

"""
    wan_read_output(filename)

Reads an outputfile for wannier.
Parsed info:
    :disentanglement,
    :wannierise,
    :final_state,
"""
function wan_read_output(filename)
    DisTuple = NamedTuple{(:Iter, :Ω_i_1, :Ω_i, :δ, :Time), Tuple{Int, Float64, Float64, Float64, Float64}}
    WanTuple = NamedTuple{(:Iter, :δ_spread, :∇RMS, :Spread, :Time), Tuple{Int, Float64, Float64, Float64, Float64}}
    WfcTuple = NamedTuple{(:index, :center, :Spread), Tuple{Int, Point3{Float64}, Float64}}
    data = Dict{Symbol, Any}()
    open(filename, "r") do f
        sreadline() = readline(f) |> strip
        line = sreadline()
        while !eof(f)
            if line == "Extraction of optimally-connected subspace"
                data[:disentanglement] = DisTuple[]
                for i = 1:4 readline(f) end
                line = sreadline()
                while !isempty(line)
                    push!(data[:disentanglement], parse(DisTuple, split(line)[1:end-2]))
                    line = sreadline()
                end
            elseif line == "Initial State"
                data[:wannierise] = WanTuple[]
                line = sreadline()
                while !occursin("Final State", line)
                    if occursin("CONV", line)
                        push!(data[:wannierise], parse(WanTuple, split(line)[1:end-2]))
                    end
                    line = sreadline()
                end
                data[:final_state] = WfcTuple[]
                line = sreadline()
                while !occursin("Sum", line)
                    line = replace(replace(replace(line, "(" => ""), ")" => ""), "," => "")
                    push!(data[:final_state], parse(WfcTuple, split(line)[5:end]))
                    line = sreadline()
                end
            end
            line = sreadline()
        end
    end
    return data
end
