
function readoutput(calculation::Calculation{Wannier90}, files...; kwargs...)
    return wan_parse_output(files[1]; kwargs...)
end

#THIS IS THE MOST HORRIBLE FUNCTION I HAVE EVER CREATED!!!
#extracts only atoms with projections
function extract_atoms(atoms_block, proj_block, cell::Mat3, spinors = false)
    if atoms_block.name == :atoms_cart
        cell = Mat3(Matrix(1.0Ang * I, 3, 3))
    end
    projections_ = proj_block.data
    atoms = atoms_block.data
    t_start = 1
    out_ats = Atom[]
    if projections_ !== nothing
        t_ats = Atom[]
        for (proj_at, projs) in projections_
            for (pos_at, pos) in atoms
                for ps in pos
                    if proj_at != pos_at
                        continue
                    end
                    for proj in projs
                        orb = Structures.orbital(proj)
                        size = spinors ? 2 * length(orb) : length(orb)
                        push!(t_ats,
                              Atom(; name = pos_at, position_cart = cell * ps,
                                   position_cryst = ps,
                                   projections = [Projection(orb,
                                                             t_start, t_start + size - 1)]))
                        t_start += size
                    end
                end
            end
        end
        for (pos_at, pos) in atoms
            for ps in pos
                same_ats = Atom[]
                for at in t_ats
                    if at.position_cart == cell * ps
                        push!(same_ats, at)
                    end
                end
                if isempty(same_ats)
                    continue
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
                push!(out_ats,
                      Atom(pos_at; position_cart = cell * p, position_cryst = p,
                           projections = Projection[]))
            end
        end
    end
    return out_ats
end

function extract_structure(cell_block, atoms_block, projections_block,
                           spinors = false)
    if atoms_block == nothing || cell_block == nothing
        return nothing
    end
    if cell_block.option == :bohr
        cell = cell_block.data' .* 1bohr
    else
        cell = cell_block.data' .* 1Ang
    end

    atoms = extract_atoms(atoms_block, projections_block, cell, spinors)
    return Structure(Mat3(cell), atoms)
end

function wan_parse_calculation(file)
    if file isa IO || !occursin("\n", file)
        contents = readlines(file)
    else
        contents = split(file, "\n")
    end
    lines = map(contents) do l
        id = findfirst(isequal('!'), l)
        if id !== nothing
            l[1:id]
        else
            l
        end
    end |> x -> filter(!isempty, x)
     
    flags       = Dict{Symbol,Any}()
    data        = InputData[]
    atoms_block = nothing
    cell_block  = nothing
    proj_block  = nothing
    i = 1
    while i < length(lines)
        line = lines[i]
        if occursin("begin", lowercase(line))
            block_name = Symbol(split(lowercase(line))[end])

            if block_name == :projections
                proj_dict = Tuple{Symbol,Array{String,1}}[]
                i += 1
                line = lines[i]
                while !occursin("end", lowercase(line))
                    if occursin("random", line)
                        proj_block = InputData(:projections, :random, nothing)
                    else
                        split_line  = strip_split(line, ':')
                        atom        = Symbol(split_line[1])
                        projections = [proj for proj in strip_split(split_line[2], ';')]
                        push!(proj_dict, (atom, projections))
                    end
                    i += 1
                    line = lines[i]
                end
                proj_block = InputData(:projections, :none, proj_dict)

            elseif block_name == :kpoint_path
                i += 1
                line = lines[i]
                k_path_array = Array{Tuple{Symbol,Array{Float64,1}},1}()
                while !occursin("end", lowercase(line))
                    split_line = split(line)
                    push!(k_path_array,
                          (Symbol(split_line[1]), parse.(Float64, split_line[2:4])))
                    push!(k_path_array,
                          (Symbol(split_line[5]), parse.(Float64, split_line[6:8])))
                    i += 1
                    line = lines[i]
                end
                push!(data, InputData(:kpoint_path, :none, k_path_array))

            elseif block_name == :unit_cell_cart
                i += 1
                line = lines[i]
                if length(split(line)) == 1
                    option = Symbol(strip(lowercase(line)))
                    i += 1
                    line = lines[i]
                else
                    option = :ang
                end
                cell_param = Matrix{Float64}(undef, 3, 3)
                for i in 1:3
                    cell_param[i, :] = parse_line(Float64, line)
                    i += 1
                    line = lines[i]
                end
                cell_block = InputData(:unit_cell_cart, option, Mat3(cell_param))
                # line = readline(f)

            elseif block_name == :atoms_frac || block_name == :atoms_cart
                i += 1
                line   = lines[i]
                atoms  = Dict{Symbol,Array{Point3{Float64},1}}()
                option = :ang
                while !occursin("end", lowercase(line))
                    if length(split(line)) == 1
                        option = Symbol(strip(line))
                        i += 1
                        line = lines[i]
                    end
                    split_line = strip_split(line)
                    atom       = Symbol(split_line[1])
                    position   = Point3(parse.(Float64, split_line[2:4]))
                    if !haskey(atoms, atom)
                        atoms[atom] = [position]
                    else
                        push!(atoms[atom], position)
                    end
                    i += 1
                    line = lines[i]
                end
                atoms_block = InputData(block_name, option, atoms)
            elseif block_name == :kpoints
                i += 1
                line = lines[i]
                k_points = Array{NTuple{3,Float64},1}()
                while !occursin("end", lowercase(line))
                    push!(k_points, (parse_line(Float64, line)...,))
                    i += 1
                    line = lines[i]
                end
                push!(data, InputData(:kpoints, :none, k_points))
            end

        else
            if occursin("mp_grid", line)
                flags[:mp_grid] = parse_string_array(Int, split(line)[end-2:end])
            else
                flag, val = wan_parse_flag_line(line)
                flags[flag] = val
            end
        end
        i += 1
        line = lines[i]
    end
    structure = extract_structure(cell_block, atoms_block, proj_block,
                                  get(flags, :spinors, false))
    return (flags=flags, data=data, structure = structure)
end

function wan_parse_array_value(eltyp, value_str)
    value_str = replace(value_str, "," => "")
    spl = split(value_str)
    arr = eltyp[]
    for i in 1:length(spl)
        if spl[i] == "t"
            spl[i] = "true"
        elseif spl[i] == "f"
            spl[i] = "false"
        end
        s = spl[i]
        if occursin("-", s)
            t = split(s, '-')
            append!(arr, collect(parse(eltyp, t[1]):parse(eltyp, t[2])))
        else
            append!(arr, parse(eltyp, s))
        end
    end
    return arr
end

function wan_parse_flag_line(line::String)
    split_line = strip_split(line, '=')
    flag       = Symbol(split_line[1])
    flagtyp    = Calculations.flagtype(Wannier90, flag)
    value      = strip(lowercase(split_line[2]), '.')
    if flagtyp != String
        value = replace(value, "d" => "e")
    end
    if flagtyp <: Vector
        parsed_val = wan_parse_array_value(eltype(flagtyp), value)
    elseif flagtyp == String
        parsed_val = value
    else
        parsed_val = parse.(eltype(flagtyp), split(value))
    end
    if length(parsed_val) == 1
        val = parsed_val[1]
    else
        val = parsed_val
    end
    return flag, val
end

function wan_write_projections(f::IO, atoms::Vector{Atom})
    write(f, "begin projections\n")
    uniats = unique(atoms)
    projs = map(x -> x.projections, uniats)
    if all(isempty.(projs))
        write(f, "random\n")
    else
        for (at, prjs) in zip(uniats, projs)
            if isempty(prjs)
                continue
            end
            write(f, Structures.projections_string(at))
            write(f, "\n")
        end
    end
    return write(f, "end projections\n")
end

"""
    write(f, calculation::Calculation{Wannier90}, structure)

Writes the `Calculation{Wannier90}` and `structure` to f.
The atoms in the structure must have projections defined.
"""
function Base.write(f::IO, calculation::Calculation{Wannier90}, structure)

    cursize = f.size
    
    preprocess = pop!(calculation, :preprocess, false)
    wplot = pop!(calculation, :wannier_plot, false)
    
    projs = vcat(map(structure.atoms) do x
                     ps = x.projections
                     # @assert !isempty(ps) "Please first set projections for all atoms in the Structure."
                     return ps
                 end...)

    Structures.sanitize!(projs, Calculations.issoc(calculation))
    for (flag, value) in calculation.flags
        write_flag_line(f, flag, value)
    end
    write(f, "\n")

    if structure != nothing
        write(f, "begin unit_cell_cart\n")
        writedlm(f, ustrip.(structure.cell'))
        write(f, "end unit_cell_cart\n")
        write(f, "\n")
    end
    wan_write_projections(f, structure.atoms)

    write(f, "\n")
    write(f, "begin atoms_frac\n")
    for at in structure.atoms
        pos = round.(at.position_cryst, digits = 5)
        write(f, "$(at.name)  $(pos[1]) $(pos[2]) $(pos[3])\n")
    end
    write(f, "end atoms_frac\n")
    write(f, "\n")

    for block in calculation.data
        write(f, "begin $(block.name)\n")
        if block.name == :kpoint_path
            for i in 1:2:length(block.data)
                letter1, k_points1 = block.data[i]
                letter2, k_points2 = block.data[i+1]
                write(f,
                      "$letter1 $(k_points1[1]) $(k_points1[2]) $(k_points1[3]) $letter2 $(k_points2[1]) $(k_points2[2]) $(k_points2[3])\n")
            end

        elseif block.name == :kpoints
            for k in block.data
                write(f, "$(k[1]) $(k[2]) $(k[3])\n")
            end
        end
        write(f, "end $(block.name)\n\n")
    end
    
    Calculations.set_flags!(calculation, :preprocess => preprocess, :wannier_plot => wplot; print=false)
    
    return f.size - cursize
end
function Base.write(f::AbstractString, calculation::Calculation{Wannier90}, structure)
    open(f, "w") do file
        write(file, calculation, structure)
    end
end

function write_flag_line(f, flag, data, seperator = "=", i = "")
    flagstr = string(flag)
    if flagstr[end-1] == '_' && tryparse(Int, string(flagstr[end])) != nothing
        flagstr = flagstr[1:end-2] * "($(flagstr[end]))"
    end
    write(f, "  $flagstr$i $seperator ")

    if typeof(data) <: Array
        if length(data) % 3 == 0 && eltype(data) != Int
            write(f, "\n")
            for i in 1:3:length(data)
                write(f, "  $(data[i]) $(data[i + 1]) $(data[i + 2])\n")
            end
        else
            write(f, "  $(data[1])")
            for x in data[2:end]
                write(f, " $x")
            end
            write(f, "\n")
        end

    else #this should work for anything singular valued data such as bools, ''s and other types
        write(f, "$data\n")
    end
end

function wan_parse_disentanglement(out, line, f)
    DisTuple = NamedTuple{(:Iter, :Ω_i_1, :Ω_i, :δ, :Time),
                          Tuple{Int,Float64,Float64,Float64,Float64}}
    out[:disentanglement] = DisTuple[]
    for i in 1:4
        readline(f)
    end
    line = strip(readline(f))
    while !isempty(line)
        push!(out[:disentanglement], parse(DisTuple, split(line)[1:end-2]))
        line = strip(readline(f))
    end
end
function wan_parse_wannierise(out, line, f)
    WanTuple = NamedTuple{(:Iter, :δ_spread, :∇RMS, :Spread, :Time),
                          Tuple{Int,Float64,Float64,Float64,Float64}}
    WfcTuple = NamedTuple{(:index, :center, :Spread),Tuple{Int,Point3{Float64},Float64}}
    out[:wannierise] = WanTuple[]
    line = strip(readline(f))
    while !occursin("Final State", line)
        if occursin("CONV", line)
            push!(out[:wannierise], parse(WanTuple, split(line)[1:end-2]))
        end
        line = strip(readline(f))
    end
    out[:final_state] = WfcTuple[]
    line = strip(readline(f))
    while !occursin("Sum", line)
        line = replace(replace(replace(line, "(" => ""), ")" => ""), "," => "")
        push!(out[:final_state], parse(WfcTuple, split(line)[5:end]))
        line = strip(readline(f))
    end
end

const WAN_PARSE_FUNCS = ["Extraction of optimally-connected subspace" => wan_parse_disentanglement,
                         "Initial State" => wan_parse_wannierise,
                         "All done:" => (x, y, z) -> x[:finished] = true]

"""
    wan_parse_output(filename::AbstractString; parse_funcs::Vector{Pair{String,Function}}=Pair{String,Function}[])

Reads an outputfile for wannier.
Parsed info:
    :disentanglement,
    :wannierise,
    :final_state,
"""
function wan_parse_output(filename;
                         parse_funcs::Vector{<:Pair{String}} = Pair{String}[])
    out = parse_file(filename, WAN_PARSE_FUNCS; extra_parse_funcs = parse_funcs)
    out[:finished] = haskey(out, :finished) ? true : false
    return out
end

function outfiles(c::Calculation{Wannier90})
    files = [c.name]
    if get(c, :wannier_plot, false)
        push!(files, "UNK")
    end
    return unique(files)
end
