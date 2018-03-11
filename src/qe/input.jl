
mutable struct QEControlBlock <: ControlBlock
    name::Symbol
    flags::Dict{Symbol,Any}
end

mutable struct QEDataBlock <: DataBlock
    name::Symbol
    option::Symbol
    data::Any
end
mutable struct QEInput <: DFInput
    filename       ::String
    control_blocks ::Vector{QEControlBlock}
    data_blocks    ::Vector{QEDataBlock}
    run_command    ::String  #everything before < in the job file
    exec           ::String
    run            ::Bool
end

"""
    QEInput(template::QEInput, filename, newflags...; run_command=template.run_command, run=true, new_data...)

Creates a new `QEInput` from a template `QEInput`, setting the newflags in the new one.
"""
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
    return found_keys, input
end

"""
    set_flags!(input::QEInput, flags...)

Sets the specified flags in the input, if they are allowed. The flag values will be converted to the correct type according to the Documentation provided by QE. A ControlBlock will be added to the input if necessary.
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
    return found_keys, input
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
    return input
end
#This can be done with a named tuple.
"""
    change_kpoints!(input::QEInput, k_grid::Union{NTuple{3, Int}, NTuple{6, Int}})

Changes the data in the k point `DataBlock` inside the specified calculation.
If the specified calculation is `'nscf'` the accepted format is `(nka, nkb, nkc)`, and the k_grid will be generated. If the calculation is `'scf'` the format is `(nka, nkb, nkc, sta, stb, stc)`.
"""
function change_kpoints!(input::QEInput, k_grid::Union{NTuple{3, Int}, NTuple{6, Int}}; print=true)
    if length(k_grid) == 3
        calc = get_flag(input, :calculation)
        if calc != "'nscf'"
            return
        end
        k_option = :crystal
        k_points = gen_k_grid(k_grid[1], k_grid[2], k_grid[3], :nscf)
    else
        calc = get_flag(input, :calculation)
        @assert calc == "'scf'" || contains(calc, "relax") warn("Expected calculation to be 'scf', 'vc-relax', 'relax'.\nGot $calc.")
        k_option = :automatic
        k_points = [k_grid...]
    end
    change_data!(input, :k_points, k_points, option = k_option, print=print)
    return input
end

"""
    change_kpoints!(input::QEInput, k_grid::Vector{NTuple{4, <:AbstractFloat}};
    k_option=:crystal_b)

Changes the data in the k point `DataBlock` inside the specified calculation. The format is `[(ka, kb, kc, nk),...]`. This format is to be used with a `'bands'` calculation.
"""
function change_kpoints!(input::QEInput, k_grid::Vector{NTuple{4, T}}; print=true, k_option=:crystal_b) where T<:AbstractFloat
    calc = get_flag(input, :calculation)
    @assert calc == "'bands'" warn("Expected calculation to be 'bands', got $calc.")
    @assert k_option in [:tpiba_b, :crystal_b, :tpiba_c, :crystal_c] error("Only $([:tpiba_b, :crystal_b, :tpiba_c, :crystal_c]...) are allowed as a k_option, got $k_option.")
    if k_option in [:tpiba_c, :crystal_c]
        @assert length(k_grid) == 3 error("If $([:tpiba_c, :crystal_c]...) is selected the length of the k_points needs to be 3, got length: $(length(k_grid)).")
    end
    k_option = k_option
    num_k = 0.0
    for k in k_grid
        num_k += k[4]
    end
    if num_k > 100.
        set_flags!(input, :verbosity => "'high'")
        if print
            dfprintln("Verbosity is set to high because num_kpoints > 100,\n
                       otherwise bands won't get printed.")
        end
    end
    change_data!(input, :k_points, k_grid, option = k_option, print = print)
    return input
end
