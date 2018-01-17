
"""
    change_flags!(input::QEInput, new_flag_data...)

Changes the flags inside the input to the new ones if they are already defined and if the new ones have the same type.
"""
function change_flags!(input::QEInput, new_flag_data...)
    found_keys = Symbol[]
    for block in input.control_blocks
        for (flag, value) in new_flag_data
            if haskey(block.flags, flag)
                old_data = block.flags[flag]
                if !(flag in found_keys) push!(found_keys, flag) end
                if typeof(block.flags[flag]) == typeof(value)
                    block.flags[flag] = value
                    dfprintln("$(input.filename):\n -> $(block.name):\n  -> $flag:\n      $old_data changed to: $value\n")
                else
                    dfprintln("$(input.filename):\n -> $(block.name):\n  -> $flag:\n     type mismatch old:$old_data ($(typeof(old_data))), new: $value) ($(typeof(value)))\n    Change not applied.\n")
                end
            end
        end
    end
    return found_keys
end

"""
    set_flags!(input::DFInput, flags...)

Sets the specified flags in the input. A controlblock will be added if necessary.
"""
function set_flags!(input::QEInput, flags...; print=true)
    found_keys = Symbol[]
    for (flag, value) in flags
        flag_block, flag_info = get_qe_block_variable(input, flag)
        flag_type = flag_info._type
        if flag_type != Void
            if !(flag in found_keys) push!(found_keys, flag) end
            
            try
                value = length(value) > 1 ? convert.(flag_type, value) : convert(flag_type, value)
            catch
                if print dfprintln("Filename '$(input.filename)':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n") end
                continue
            end
            
            input_block = get_block(input, flag_block.name)
            if input_block != nothing
                if haskey(input_block.flags, flag)
                    old_data = input_block.flags[flag]
                    input_block.flags[flag] = value
                    if print dfprintln("$(input.filename):\n -> $(input_block.name):\n  -> $flag:\n      $old_data changed to: $value\n") end
                else
                    input_block.flags[flag] = value
                    if print dfprintln("$(input.filename):\n -> $(input_block.name):\n  -> $flag:\n      set to: $value\n") end
                end
            else
                push!(input.control_blocks, QEControlBlock(flag_block.name, OrderedDict(flag => value)))
                if print dfprintln("$(input.filename):\n -> New ControlBlock $flag_block:\n  -> $flag:\n      set to: $value\n") end
            end
        end
    end
    return found_keys
end

"""
    get_flag(input::QEInput, flag::Symbol)

Returns the value of the flag.
"""
function get_flag(input::QEInput, flag::Symbol)
    for block in input.control_blocks
        if haskey(block.flags, flag)
            return block.flags[flag]
        end
    end
end

"""
    get_block(input::QEInput, block_symbol::Symbol)

Returns the block with name `block_symbol`.
"""
function get_block(input::QEInput, block_symbol::Symbol)
    for block in [input.control_blocks; input.data_blocks]
        if block.name == block_symbol
            return block
        end
    end
    return nothing
end

"""
    remove_flags!(input::QEInput, flags...)

Remove the specified flags.
"""
function remove_flags!(input::QEInput, flags...)
    for (i, block) in enumerate(input.control_blocks)
        for flag in flags
            if haskey(block.flags, flag)
                pop!(block.flags, flag)
                dfprintln("Removed flag '$flag' from block '$(block.name)' in input '$(input.filename)'")
                if isempty(block.flags)
                    input.control_blocks = length(input.control_blocks) > i ? [input.control_blocks[1:i - 1]; input.control_blocks[i + 1,end]] : input.control_blocks[1:i - 1]
                    dfprintln("Removed block '$(block.name)' in input '$(input.filename)', because it was empty.'")
                    break
                end
            end
        end
    end
end

"""
    print_flags(input::QEInput, block_symbol::Symbol)

Prints the flags of the specified block.
"""
function print_flags(input::QEInput, block_symbol::Symbol)
    block = filter(x -> x.name == block_symbol, input.control_blocks)[1]
    dfprintln("  $(block.name):")
    for (flag, value) in block.flags
        dfprintln("    $flag => $value")
    end
    dfprintln("")
end

#This can be done with a named tuple.
"""
    change_kpoints!(input::QEInput, k_grid::Union{NTuple{3, Int}, NTuple{6, Int}})

Changes the data in the k point `DataBlock` inside the specified calculation.
"""
function change_kpoints!(input::QEInput, k_grid::Union{NTuple{3, Int}, NTuple{6, Int}}; print=true)
    if length(k_grid) == 3
        calc = get_flag(input, :calculation) 
        @assert calc == "'nscf'" warn("Expected calculation to be 'nscf', got $calc.")
        k_option = :crystal
        k_points = gen_k_grid(k_grid[1], k_grid[2], k_grid[3], :nscf)
    else
        calc = get_flag(input, :calculation)
        @assert calc == "'scf'" || contains(calc, "relax") warn("Expected calculation to be 'scf', 'vc-relax', 'relax'.\nGot $calc.")
        k_option = :automatic
        k_points = [k_grid...]
    end
    change_data!(input, :k_points, k_points, option = k_option, print=print)
end

"""
    change_kpoints!(input::QEInput, k_grid::Array{Array{<:AbstractFloat, 1}, 1})

Changes the data in the k point `DataBlock` inside the specified calculation.
"""
function change_kpoints!(input::QEInput, k_grid::Array{Array{<:AbstractFloat, 1}, 1}; print=true)
    calc = get_flag(input, :calculation) 
    @assert calc == "'bands'" warn("Expected calculation to be 'bands', got $calc.")
    k_option = :crystal_b
    num_k = 0.0
    for k in k_grid
        num_k += k[4]
    end
    if num_k > 100.
        set_flags!(input, :verbosity => "'high'")
        if print 
            dfprintln("Set verbosity to high because num_kpoints > 100,\n
                       otherwise bands won't get printed.")
        end
    end
    change_data!(input, :k_points, k_grid, option = k_option, print)
end

"Returns the cell_parameters in Angstrom."
function get_cell(input::QEInput)
    block = get_block(input, :cell_parameters)
    if block != nothing
        if block.option == :bohr
            alat = conversions[:bohr2ang]
        elseif block.option == :alat
            alat = get_flag(input,:A) == nothing ? get_flag(input, Symbol("celldm(1)")) * conversions[:bohr2ang] : get_flag(input,:A) 
        else
            alat = 1.0
        end
        return block.data * alat
    end
end

function change_cell!(input::QEInput, cell_parameters::Matrix; option=:angstrom, print=true)
    @assert size(cell_parameters) == (3,3) "Cell parameters has wrong size.\nExpected: (3,3) got ($(size(cell_parameters)[1]), $(size(cell_parameters)[2]))."
    change_data!(input, :cell_parameters, cell_parameters, option=option, print=print)
    remove_flags!(input, [:A, Symbol("celldm(1)")])
end


"""
    get_atoms(input::QEInput)

Returns a list of the atomic positions in Angstrom. We only support A or celldm(1) as alat.
"""
function get_atoms(input::QEInput)
    out_atoms = OrderedDict{Symbol, Array{<:Point3D, 1}}()
    block = get_block(input, :atomic_positions)
    if block != nothing
        if block.option == :alat
            alat = get_flag(input,:A) == nothing ? get_flag(input, Symbol("celldm(1)")) * conversions[:bohr2ang] : get_flag(input,:A)
            cell = eye(3) * alat
        elseif block.option == :bohr
            cell = eye(3) * conversions[:bohr2ang]
        elseif block.option == :crystal || block.option == :crystal_sg
            cell = get_cell(input)
        else
            cell = eye(3)
        end
        for (atom, positions) in block.data
            t_pos = Point3D[]
            for pos in positions
                push!(t_pos, cell' * pos) 
            end
            out_atoms[atom] = t_pos
        end
    end
    return out_atoms
end

"""
    change_atoms!(input::QEInput, atoms::OrderedDict{Symbol,<:Array{<:Point3D,1}}; option=:angstrom, pseudos_set=nothing, pseudo_specifier=nothing)

Changes the atoms in the input to the specified atoms. Also sets up the pseudos if possible.
If 'default_pseudos' is defined it will look for the pseudo set and also write the correct values inside the 'DataBlock' that defines which pseudopotentials to use.
"""
function change_atoms!(input::QEInput, atoms::OrderedDict{Symbol,<:Array{<:Point3D,1}}; 
                       option           = :angstrom,
                       pseudo_set       = nothing, 
                       pseudo_specifier = nothing,
                       print            = true)

    changed = change_data!(input, :atomic_positions, atoms, option=option, print=print)
    if !changed
        return
    else
        nat = sum([length(x) for x in values(atoms)])
        ntyp = length(keys(atoms))
        set_flags!(input, :nat => nat, :ntyp => ntyp, print=print)
    end
    if isdefined(:default_pseudos) && pseudo_set != nothing
        atomic_species_OrderedDict = OrderedDict{Symbol,String}()
        for atom in keys(atoms)
            atomic_species_OrderedDict[atom] = get_default_pseudo(atom, pseudo_set, pseudo_fuzzy=pseudo_fuzzy)
        end
        change_data!(input, :atomic_species, atomic_species_OrderedDict)
        set_flags!(input, :pseudo_dir => "'$(default_pseudo_dirs[pseudo_set])'",print=print)
    end
end

