# All the methods to change the inpût control flags, if you want to implement another kind of calculation add a similar one here!
include("qe/input.jl")
include("wannier90/input.jl")
"""
    change_flags!(input::DFInput, new_flag_data...)

Changes the flags inside the input to the new ones if they are already defined and if the new ones have the same type.
"""
function change_flags!(input::Union{AbinitInput, WannierInput}, new_flag_data...)
    found_keys = Symbol[]
    for (flag, value) in new_flag_data
        if haskey(input.flags, flag)
            old_data = input.flags[flag]
            if !(flag in found_keys) push!(found_keys, flag) end
            if typeof(input.flags[flag]) == typeof(value)
                input.flags[flag] = value
                dfprintln("$(input.filename):\n -> $flag:\n      $old_data changed to: $value\n")
            else
                dfprintln("$(input.filename):\n -> $flag:\n     type mismatch old:$old_data ($(typeof(old_data))), new: $value ($(typeof(value)))\n    Change not applied.\n")
            end
        end
    end
    return found_keys
end

"""
    set_flags!(input::DFInput, flags...)

Sets the specified flags in the input. A controlblock will be added if necessary.
"""
function set_flags!(input::Union{AbinitInput, WannierInput}, flags...)
    found_keys = Symbol[]
    flag_func(flag) = typeof(input) == WannierInput ? get_wan_flag_type(flag) : get_abi_flag_type(flag)
    for (flag, value) in flags
        flag_type = flag_func(flag)
        if flag_type != Void
            if !(flag in found_keys) push!(found_keys, flag) end
            try
                value = length(value) > 1 ? convert.(flag_type, value) : convert(flag_type, value)
            catch
                dfprintln("Filename '$(input.filename)':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            
            if haskey(input.flags, flag)
                old_data = input.flags[flag]
                input.flags[flag] = value
                dfprintln("$(input.filename):\n  -> $flag:\n      $old_data changed to: $value\n")
            else
                input.flags[flag] = value
                dfprintln("$(input.filename):\n  -> $flag:\n      set to: $value\n")
            end
        end
    end
    return found_keys
end


"""
    change_data!(input::DFInput, block_name::Symbol, new_block_data; option=nothing)

Changes the data of the specified 'DataBlock' to the new data. Optionally also changes the 'DataBlock' option.
"""
function change_data!(input::DFInput, block_name::Symbol, new_block_data; option=nothing)
    changed = false
    for data_block in input.data_blocks
        if data_block.name == block_name
            if typeof(data_block.data) != typeof(new_block_data) 
                warn("Overwritten data of type '$(typeof(data_block.data))' with type '$(typeof(new_block_data))'.")
            end
            old_data        = data_block.data
            data_block.data = new_block_data
            dfprintln("Block data '$(data_block.name)' in input  '$(input.filename)' is now:")
            dfprintln(string(data_block.data))
            dfprintln("")
            data_block.option = option == nothing ? data_block.option : option
            dfprintln("option: $(data_block.option)")
            changed = true
        end
    end
    return changed
end

"""
    add_data!(input::DFInput, block_name::Symbol, block_data, block_option=:none)

Adds a block with the given name data and option to the calculation.
"""
function add_data!(input::DFInput, block_name::Symbol, block_data, block_option=:none)
    for block in input.data_blocks
        if block_name == block.name
            change_data!(input, block_name, block_data, option=block_option) 
            return
        end
    end
    if typeof(input)==QEInput
        block = QEDataBlock(block_name, block_option, block_data)
    elseif typeof(input) == WannierInput
        block = WannierDataBlock(block_name, block_option, block_data)
    elseif typeof(input) == AbinitInput
        block = AbinitDataBlock(block_name, block_option, block_data)
    end
    add_block!(input, block)
end

"""
    get_flag(input::DFInput, flag::Symbol)

Returns the value of the flag.
"""
function get_flag(input::Union{AbinitInput, WannierInput}, flag::Symbol)
    if haskey(input.flags, flag)
        return input.flags[flag]
    end
end

"""
    get_data(input::DFInput, block_symbol::Symbol)

Returns the specified 'DataBlock'.
"""
function get_data(input::DFInput, block_symbol::Symbol)
    block = get_block(input, block_symbol)
    if block != nothing && typeof(block) <: DataBlock
        return block.data
    else 
        error("No `DataBlock` with name '$block_symbol' found. ")
    end
end

"""
    get_block(input::DFInput, block_symbol::Symbol)

Returns the block with name `block_symbol`.
"""
function get_block(input::Union{AbinitInput, WannierInput}, block_symbol::Symbol)
    for block in input.data_blocks
        if block.name == block_symbol
            return block
        end
    end
    return nothing
end

#removes an input control flag, if you want to implement another input add a similar function here!
"""
    remove_flags!(input::DFInput, flags...)

Remove the specified flags.
"""
function remove_flags!(input::Union{AbinitInput, WannierInput}, flags...)
    for flag in flags
        if haskey(input.flags, flag)
            pop!(input.flags, flag, false)
            dfprintln("Removed flag '$flag' from input '$(input.filename)'")
        end
    end
end

"""
    get_blocks(input::DFInput) 

Returns all the blocks inside a calculation.
"""
function get_blocks(input::DFInput) 
    out = Block[]
    for block in getfield.(input, filter(x -> contains(String(x), "block"), fieldnames(input)))
        push!(out, block...)
    end
    return out 
end

"""
    add_block(input::DFInput, block::Block)

Adds the given block to the input. Should put it in the correct arrays.
"""
function add_block!(input::DFInput, block::Block)
    if typeof(block) <: DataBlock
        push!(input.data_blocks, block)
    else
        typeof(block) <: ControlBlock
        push!(input.control_blocks, block)
    end
    print_blocks(input)
end

"""
    print_flag(input::DFInput, flag)

Prints information of the flag inside the input.
"""
function print_flag(input::DFInput, flag)
    if (:control_blocks in fieldnames(input))
        for block in input.control_blocks
            if haskey(block.flags, flag)
                s = """ 
                Filename: $(input.filename)  
                Block Name: $(block.name)
                $flag => $(block.flags[flag])\n
                """
                dfprintln(s)
            end
        end
    end

    if (:flags in fieldnames(input))
        if haskey(input.flags, flag)
            s = """
            Filename: $(input.filename)
            $flag => $(input.flags[flag])\n
            """
            dfprintln(s)
        end
    end
end

print_flags(input::DFInput, flags::Array) = print_flag.(input, flags)

"""
    print_block(input::DFInput, block_name::Symbol)

Print the information of a 'Block' inside the input.
"""
function print_block(input::DFInput, block_name::Symbol)
    found        = false
    input_blocks = get_blocks(input)
    for block in input_blocks
        if block.name == block_name
            print_filename(input)
            display(block) 
            found = true
        end
    end
    return found
end

"""
    print_blocks(input::DFInput)

Print the information on all Blocks inside the input.
"""
function print_blocks(input::DFInput)
    print_filename(input)
    get_blocks(input) |> display
    dfprintln("")
end

"""
    print_filename(input::DFInput)

Prints the filename associated with the input.
"""
function print_filename(input::DFInput)
    dfprintln("Input file: $(input.filename)")
end

"""
    print_info(input::DFInput)

Prints general info of the input.
"""
function print_info(input::DFInput)
    dfprintln("Filename: $(input.filename)")
    if (:control_blocks in fieldnames(input))
        dfprintln("  Control Blocks:")
        for (i, block) in enumerate(input.control_blocks)
            dfprintln("    $i: $(block.name)")
        end
    end
    dfprintln("  Data Blocks:")
    for (i, block) in enumerate(input.data_blocks)
        dfprintln("    $i: $(block.name)")
    end
    dfprintln("  Run command: $(input.run_command)")
    dfprintln("  Runs: $(input.run)")
end

"""
    print_flags(input::DFInput)

Prints all the flags of the input.
"""
function print_flags(input::DFInput)
    dfprintln("#----------------#")
    dfprintln("Filename: $(input.filename)")
    if (:control_blocks in fieldnames(input))
        for block in input.control_blocks
            dfprintln("  $(block.name):")
            for (flag, value) in block.flags
                dfprintln("    $flag => $value")
            end
            dfprintln("")
        end
    end
    if (:flags in fieldnames(input))
        for (flag, value) in input.flags
            dfprintln("  $flag => $value")
        end
    end
    dfprintln("#----------------#\n")
end

"""
    change_cell_parameters!(input::DFInput, cell_param::Array{AbstractFloat,2})

Changes the cell parameters `DataBlock`.
"""
function change_cell_parameters!(input::DFInput, cell_param::Array{<:AbstractFloat,2})
    assert(size(cell_param) == (3, 3))
    if typeof(input) == WannierInput
        change_data!(input, :unit_cell_cart, cell_param)
    elseif typeof(input) == QEInput
        change_data!(input, :cell_parameters, cell_param)
    end
end

"""
    change_k_points!(input::DFInput, calc_filename, k_points)

Changes the data in the k point `DataBlock` inside the specified calculation.
"""
function change_k_points!(input::DFInput, k_points)
    if typeof(input) == WannierInput
        change_data!(input, :kpoints, k_points)
    elseif typeof(input) == QEInput
        change_data!(input, :k_points, k_points)
    end
end


"""
    change_data_option!(job::DFJob, block_symbol::Symbol, option::Symbol)

Changes the option of specified data block.
"""
function change_data_option!(input::DFInput, block_symbol::Symbol, option::Symbol)
    for fieldname in fieldnames(input)
        field = getfield(input,fieldname)
        if typeof(field) <: Array{<:DataBlock,1}
            for block in field
                if block.name == block_symbol
                    old_option   = block.option
                    block.option = option
                    dfprintln("Option of DataBlock '$(block.name)' in input '$(input.filename)' changed from '$old_option' to '$option'")
                end
            end
        end
    end
end