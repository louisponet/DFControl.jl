#these are all the control blocks, they hold the flags that guide the calculation
abstract type Block end
abstract type ControlBlock <: Block end

#these are all the data blocks, they hold the specific data for the calculation
abstract type DataBlock <: Block end

mutable struct WannierDataBlock <: DataBlock
    name::Symbol
    option::Symbol
    data::Any
end

mutable struct AbinitDataBlock <: DataBlock
    name::Symbol
    option::Symbol
    data::Any
end

function block(blocks::Vector{<:Block}, name::Symbol)
    found_blocks = filter(x-> x.name == name, blocks)
    if isempty(found_blocks)
        return nothing
    else
        return found_blocks[1]
    end
end
"""
Represents an input for DFT calculation.
"""
abstract type DFInput end

include("qe/input.jl")

mutable struct WannierInput <: DFInput
    filename    ::String
    flags       ::Dict{Symbol,Any}
    data ::Vector{WannierDataBlock}
    runcommand ::Exec #runcommand, flags
    exec        ::Exec #exec, flags
    run         ::Bool
end

"""
    setkpoints!(input::WannierInput, k_grid::NTuple{3, Int}; print=true)

sets the data in the k point `DataBlock` inside the specified calculation.
"""
function setkpoints!(input::WannierInput, k_grid::NTuple{3, Int}; print=true)
    setflags!(input, :mp_grid => [k_grid...])
    k_points = kgrid(k_grid[1], k_grid[2], k_grid[3], :wan)
    setdata!(input, :kpoints, k_points, print=print)
    return input
end

mutable struct AbinitInput <: DFInput
    filename    ::String
    structure   ::Union{Structure, Void}
    flags       ::Dict{Symbol,Any}
    data        ::Vector{AbinitDataBlock}
    runcommand ::Exec
    exec        ::Exec
    run         ::Bool
end

"""
    flag(input::DFInput, flag::Symbol)

Returns the value of the flag.
"""
function flag(input::Union{AbinitInput, WannierInput}, flag::Symbol)
    if haskey(input.flags, flag)
        return input.flags[flag]
    end
end

"""
    setflags!(input::Union{AbinitInput, WannierInput}, flags...; print=true)

Sets the specified flags in the input.
"""
function setflags!(input::Union{AbinitInput, WannierInput}, flags...; print=true)
    found_keys = Symbol[]
    flag_func(flag) = typeof(input) == WannierInput ? wan_flag_type(flag) : abi_flag_type(flag)
    for (flag, value) in flags
        flag_type = flag_func(flag)
        if flag_type != Void
            if !(flag in found_keys) push!(found_keys, flag) end
            try
                value = length(value) > 1 ? convert.(flag_type, value) : convert(flag_type, value)
            catch
                if print dfprintln("Filename '$(input.filename)':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n") end
                continue
            end

            if haskey(input.flags, flag)
                old_data = input.flags[flag]
                input.flags[flag] = value
                if print dfprintln("$(input.filename):\n  -> $flag:\n      $old_data setd to: $value\n") end
            else
                input.flags[flag] = value
                if print dfprintln("$(input.filename):\n  -> $flag:\n      set to: $value\n") end
            end
        end
    end
    return found_keys, input
end

"""
    rmflags!(input::Union{AbinitInput, WannierInput}, flags...)

Remove the specified flags.
"""
function rmflags!(input::Union{AbinitInput, WannierInput}, flags...)
    for flag in flags
        if haskey(input.flags, flag)
            pop!(input.flags, flag, false)
            dfprintln("Removed flag '$flag' from input '$(input.filename)'")
        end
    end
    return input
end


"""
    setdata!(input::DFInput, block_name::Symbol, new_block_data; option=nothing, print=true)

sets the data of the specified 'DataBlock' to the new data. Optionally also sets the 'DataBlock' option.
"""
function setdata!(input::DFInput, block_name::Symbol, new_block_data; option=nothing, print=true)
    setd = false
    for data_block in input.data
        if data_block.name == block_name
            if typeof(data_block.data) != typeof(new_block_data)
                if print warn("Overwritten data of type '$(typeof(data_block.data))' with type '$(typeof(new_block_data))'.") end
            end
            old_data        = data_block.data
            data_block.data = new_block_data
            data_block.option = option == nothing ? data_block.option : option
            if print
                dfprintln("Block data '$(data_block.name)' in input  '$(input.filename)' is now:")
                dfprintln(string(data_block.data))
                dfprintln("option: $(data_block.option)")
                dfprintln("")
            end
            setd = true
        end
    end
    if !setd
        if typeof(input)==QEInput
            block = QEDataBlock(block_name, block_option, block_data)
        elseif typeof(input) == WannierInput
            block = WannierDataBlock(block_name, block_option, block_data)
        elseif typeof(input) == AbinitInput
            block = AbinitDataBlock(block_name, block_option, block_data)
        end
        setblock!(input, block)
        setd = true
    end
    return setd, input
end

"""
    data(input::DFInput, block_symbol::Symbol)

Returns the specified 'DataBlock'.
"""
function data(input::DFInput, block_symbol::Symbol)
    block_ = block(input, block_symbol)
    if block_ != nothing && typeof(block_) <: DataBlock
        return block_.data
    else
        error("No `DataBlock` with name '$block_symbol' found. ")
    end
end

"""
    block(input::DFInput, block_symbol::Symbol)

Returns the block with name `block_symbol`.
"""
function block(input::Union{AbinitInput, WannierInput}, block_symbol::Symbol)
    for block_ in input.data
        if block_.name == block_symbol
            return block_
        end
    end
    return nothing
end

"""
    blocks(input::DFInput)

Returns all the blocks inside a calculation.
"""
function blocks(input::DFInput)
    out = Block[]
    for block in getfield.(input, filter(x -> contains(String(x), "block"), fieldnames(input)))
        push!(out, block...)
    end
    return out
end

function setoradd!(blocks::Vector{<:Block}, block::Block)
    found = false
    for (i, d) in enumerate(blocks)
        if d.name == block.name
            blocks[i] = block
            found = true
            break
        end
    end
    if !found
        push!(blocks, block)
    end
end

"""
    setblock!(input::DFInput, block::Block)

Adds the given block to the input. Should put it in the correct arrays.
"""
function setblock!(input::DFInput, block::Block)
    if typeof(block) <: DataBlock
        setoradd!(input.data, block)
    elseif typeof(block) <: ControlBlock
        setoradd!(input.control, block)
    end
    return input
end

"""
    setoption!(input::DFInput, block_symbol::Symbol, option::Symbol;; print=true)

sets the option of specified data block.
"""
function setoption!(input::DFInput, block_symbol::Symbol, option::Symbol; print=true)
    for fieldname in fieldnames(input)
        field = getfield(input,fieldname)
        if typeof(field) <: Array{<:DataBlock,1}
            for block in field
                if block.name == block_symbol
                    old_option   = block.option
                    block.option = option
                    if print dfprintln("Option of DataBlock '$(block.name)' in input '$(input.filename)' setd from '$old_option' to '$option'") end
                end
            end
        end
    end
    return input
end
