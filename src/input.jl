#these are all the control blocks, they hold the flags that guide the calculation
struct InputData
    name   ::Symbol
    option ::Symbol
    data   ::Any
end

mutable struct DFInput{P <: Package}
    filename    ::String
    flags       ::Dict{Symbol, Any}
    data        ::Vector{InputData}
    runcommand  ::Exec
    exec        ::Exec
    run         ::Bool
end

"""
    DFInput(template::DFInput, filename, newflags...; runcommand=template.runcommand, run=true)

Creates a new `DFInput` from a template `DFInput`, setting the newflags in the new one.
"""
function DFInput(template::DFInput, filename, newflags...; runcommand=template.runcommand, run=true)
    newflags = Dict(newflags...) # this should handle both OrderedDicts and pairs of flags

    input             = deepcopy(template)
    input.filename    = filename
    input.runcommand  = runcommand
    setflags!(input, newflags...)

    return input
end

inputdata(input, name) = getfirst(x-> x.name == name, input.data)
data(input, name)      = inputdata(input, name).data

"""
    setkpoints!(input::DFInput{:Wannier90}, k_grid::NTuple{3, Int}; print=true)

sets the data in the k point `InputData` inside the specified calculation.
"""
function setkpoints!(input::DFInput{Wannier90}, k_grid::NTuple{3, Int}; print=true)
    setflags!(input, :mp_grid => [k_grid...])
    setdata!(input, :kpoints, kgrid(k_grid..., :wan), print=print)
    return input
end

"""
    setkpoints!(input::DFInput{QE}, k_grid::Union{NTuple{3, Int}, NTuple{6, Int}})

sets the data in the k point `InputData` inside the specified calculation.
If the specified calculation is `'nscf'` the accepted format is `(nka, nkb, nkc)`, and the k_grid will be generated. If the calculation is `'scf'` the format is `(nka, nkb, nkc, sta, stb, stc)`.
"""
function setkpoints!(input::DFInput{QE}, k_grid::NTuple{3, Int}; print=true) #nscf

    calc = flag(input, :calculation)
    @assert calc == "'nscf'" warn("Expected calculation to be 'nscf'.\nGot $calc.")
    setdata!(input, :k_points, kgrid(k_grid..., :nscf), option = :crystal, print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::NTuple{6, Int}; print=true) #scf
    calc = flag(input, :calculation)
    @assert calc == "'scf'" || contains(calc, "relax") warn("Expected calculation to be 'scf', 'vc-relax', 'relax'.\nGot $calc.")

    setdata!(input, :k_points, [k_grid...], option = :automatic, print=print)
    return input
end

"""
    setkpoints!(input::DFInput{QE}, k_grid::Vector{NTuple{4, <:AbstractFloat}};
    k_option=:crystal_b)

sets the data in the k point `DataBlock` inside the specified calculation. The format is `[(ka, kb, kc, nk),...]`. This format is to be used with a `'bands'` calculation.
"""
function setkpoints!(input::DFInput{QE}, k_grid::Vector{NTuple{4, T}}; print=true, k_option=:crystal_b) where T<:AbstractFloat
    calc = flag(input, :calculation)
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
        setflags!(input, :verbosity => "'high'")
        if print
            dfprintln("Verbosity is set to high because num_kpoints > 100,\n
                       otherwise bands won't get printed.")
        end
    end
    setdata!(input, :k_points, k_grid, option = k_option, print = print)
    return input
end

# mutable struct AbinitInput <: DFInput
#     filename    ::String
#     structure   ::Union{Structure, Void}
#     flags       ::Dict{Symbol,Any}
#     data        ::Vector{AbinitDataBlock}
#     runcommand ::Exec
#     exec        ::Exec
#     run         ::Bool
# end
#
"""
    flag(input::DFInput, flag::Symbol)

Returns the value of the flag.
"""
function flag(input::DFInput, flag::Symbol)
    if haskey(input.flags, flag)
        return input.flags[flag]
    end
end

"""
    setflags!(input::DFInput, flags...; print=true)

Sets the specified flags in the input.
"""
function setflags!(input::DFInput{T}, flags...; print=true) where T
    found_keys = Symbol[]
    for (flag, value) in flags
        flag_type = flagtype(input, flag)
        if flag_type != Void
            !(flag in found_keys) && push!(found_keys, flag)
            try
                value = length(value) > 1 ? convert.(flag_type, value) : convert(flag_type, value)
            catch
                print && dfprintln("Filename '$(input.filename)':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            old_data = haskey(input.flags, flag) ? input.flags[flag] : ""
            print && dfprintln("$(input.filename):\n  -> $flag:\n      $old_data set to: $value\n")
        end
    end
    return found_keys, input
end

"""
    rmflags!(input::DFInput, flags...)

Remove the specified flags.
"""
function rmflags!(input::DFInput, flags...)
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

sets the data of the specified 'InputData' to the new data. Optionally also sets the 'InputData' option.
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
        setinputdata!(input, InputData(block_name, block_option, block_data))
        setd = true
    end
    return setd, input
end

#todo rename this
blocks(input) = input.data

function setoradd!(datas::Vector{InputData}, data::InputData)
    found = false
    for (i, d) in enumerate(datas)
        if d.name == data.name
            datas[i] = data
            found = true
            break
        end
    end
    if !found
        push!(datas, data)
    end
end

"""
    setinputdata!(input::DFInput, data::InputData)

Adds the given data to the input. Should put it in the correct arrays.
"""
function setinputdata!(input::DFInput, data::InputData)
    setoradd!(input.data, data)
    return input
end

"""
    setdataoption!(input::DFInput, name::Symbol, option::Symbol;; print=true)

Sets the option of specified inputdata.
"""
function setdataoption!(input::DFInput, name::Symbol, option::Symbol; print=true)
    for data in input.data
        if data.name == name
            old_option  = data.option
            data.option = option
            if print dfprintln("Option of InputData '$(data.name)' in input '$(input.filename)' set from '$old_option' to '$option'") end
        end
    end
    return input
end
