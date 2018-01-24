mutable struct QEInput <: DFInput
    filename       ::String
    structure      ::Union{Structure, Void} 
    control_blocks ::Array{QEControlBlock,1}
    data_blocks    ::Array{QEDataBlock,1}
    run_command    ::String  #everything before < in the job file
    exec           ::String
    run            ::Bool
end

function QEInput(template::QEInput, filename, newflags...; run_command=template.run_command, run=true, new_data...)
    newflags = Dict(newflags...) # this should handle both OrderedDicts and pairs of flags

    input             = deepcopy(template)
    input.filename    = filename
    input.run_command = run_command
    set_flags!(input, newflags...)

    for (block_name, block_info) in new_data
        if get_block(input, block_name) != nothing
            block = get_block(input, block_name)
            if length(block_info) == 1
                block.option = :none 
                block.data   = block_info
            elseif length(block_info) == 2
                block.option = block_info[1]
                block.data   = block_info[2]
            end
        else
            if length(block_info) == 1
                add_block!(input, QEDataBlock(block_name, :none, block_info))
            elseif length(block_info) == 2
                add_block!(input, QEDataBlock(block_name, block_info[1], block_info[2]))
            end
        end
    end
    return input
end

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
        flag_type = flag_info.typ
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
                push!(input.control_blocks, QEControlBlock(flag_block.name, Dict(flag => value)))
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
    block = getfirst(x -> x.name == block_symbol, input.control_blocks)
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
