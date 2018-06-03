#these are all the control data, they hold the flags that guide the calculation
mutable struct InputData
    name   ::Symbol
    option ::Symbol
    data   ::Any
end

mutable struct DFInput{P <: Package}
    filename ::String
    flags    ::Dict{Symbol, Any}
    data     ::Vector{InputData}
    execs    ::Vector{Exec}
    run      ::Bool
end

"""
    DFInput(template::DFInput, filename, newflags...; runcommand=template.runcommand, run=true)

Creates a new `DFInput` from a template `DFInput`, setting the newflags in the new one.
"""
function DFInput(template::DFInput, filename, newflags...; excs=execs(template), run=true, data=nothing)
    newflags = Dict(newflags...)

    input          = deepcopy(template)
    input.filename = filename
    input.execs    = excs
    setflags!(input, newflags..., print=false)

    if data != nothing
        for (name, (option, data)) in data
            setdata!(input, name, data, option=option, print=false)
        end
    end
    return input
end

filename(input::DFInput) = input.filename
data(input::DFInput)     = input.data
flags(input::DFInput)    = input.flags
function flag(input::DFInput, flag::Symbol)
    if haskey(input.flags, flag)
        return input.flags[flag]
    end
end

Base.eltype(::DFInput{P}) where P = P
package(::DFInput{P}) where P = P
data(input::DFInput, name) = getfirst(x-> x.name == name, input.data)
data(input::Vector{InputData}, name) = getfirst(x-> x.name == name, input.data)
data(input, name)      = data(input, name).data

exec(input::DFInput, exec::String) = getfirst(x -> contains(x.exec, exec), input.execs)
execs(input::DFInput) = input.execs
execs(input::DFInput, exec::String) = filter(x -> contains(x.exec, exec), input.execs)
execflags(input::DFInput, exec::String) = [x.flags for x in execs(input, exec)]
setexecflags!(input::DFInput, exec::String, flags...) = setflags!.(execs(input, exec), flags...)
setexecdir!(input::DFInput, exec, dir) = setexecdir!.(execs(input, exec), dir)

runcommand(input::DFInput) = input.execs[1]

outfile(input::DFInput{QE})        = splitext(input.filename)[1]*".out"
outfile(input::DFInput{Wannier90}) = splitext(input.filename)[1]*".wout"

setflow!(input::DFInput, run) = input.run = run


"""
    setkpoints!(input::DFInput, k_grid)

Sets the kpoints of the input. Will automatically generate the kgrid values if necessary.
"""
function setkpoints!(input::DFInput{Wannier90}, k_grid::NTuple{3, Int}; print=true)
    setflags!(input, :mp_grid => [k_grid...], print=print)
    setdata!(input, :kpoints, kgrid(k_grid..., :wan), print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::NTuple{3, Int}; print=true) #nscf

    calc = flag(input, :calculation)
    print && calc != "'nscf'" && warn("Expected calculation to be 'nscf'.\nGot $calc.")
    setdata!(input, :k_points, kgrid(k_grid..., :nscf), option = :crystal, print=print)
    prod(k_grid) > 100 && setflags!(input, :verbosity => "'high'", print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::NTuple{6, Int}; print=true) #scf
    calc = flag(input, :calculation)
    print && (calc != "'scf'" || !contains(calc, "relax")) && warn("Expected calculation to be 'scf', 'vc-relax', 'relax'.\nGot $calc.")
    setdata!(input, :k_points, [k_grid...], option = :automatic, print=print)
    prod(k_grid[1:3]) > 100 && setflags!(input, :verbosity => "'high'", print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::Vector{NTuple{4, T}}; print=true, k_option=:crystal_b) where T<:AbstractFloat
    calc = flag(input, :calculation)
    print && calc != "'bands'" && warn("Expected calculation to be 'bands', got $calc.")
    @assert in(k_option, [:tpiba_b, :crystal_b, :tpiba_c, :crystal_c]) error("Only $([:tpiba_b, :crystal_b, :tpiba_c, :crystal_c]...) are allowed as a k_option, got $k_option.")
    if k_option in [:tpiba_c, :crystal_c]
        @assert length(k_grid) == 3 error("If $([:tpiba_c, :crystal_c]...) is selected the length of the k_points needs to be 3, got length: $(length(k_grid)).")
    end
    num_k = 0.0
    for k in k_grid
        num_k += k[4]
    end
    if num_k > 100.
        setflags!(input, :verbosity => "'high'", print=print)
        if print
            info("Verbosity is set to high because num_kpoints > 100,\n
                       otherwise bands won't get printed.")
        end
    end
    setdata!(input, :k_points, k_grid, option=k_option, print=print)
    return input
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
                value = convertflag(flag_type, value)
            catch
                print && warn("Filename '$(input.filename)':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            old_data = haskey(input.flags, flag) ? input.flags[flag] : ""
            input.flags[flag] = value
            print && warn("$(input.filename):\n  -> $flag:\n      $old_data set to: $value\n")
        end
    end
    return found_keys, input
end

"Runs through all the set flags and checks if they are allowed and set to the correct value"
function sanitizeflags!(input::DFInput)
    for (flag, value) in flags(input)
        flagtype_ = flagtype(input, flag)
        if flagtype_ == Void
            info("Flag $flag was not found in allowed flags for exec $(execs(input)[2]). Removing flag.")
            rmflags!(input, flag)
            continue
        end
        if !(typeof(value) <: flagtype_)
            try
                flags(input)[flag] = convertflag(flagtype_, value)
            catch
                error("Input $(filename(input)): Could not convert :$flag of value $value to the correct type ($flagtype_), please set it to the correct type.")
            end
        end
    end
end

"""
    rmflags!(input::DFInput, flags...)

Remove the specified flags.
"""
function rmflags!(input::DFInput, flags...; print=true)
    for flag in flags
        if haskey(input.flags, flag)
            pop!(input.flags, flag, false)
            print && info("Removed flag '$flag' from input '$(input.filename)'")
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
                info("Block data '$(data_block.name)' in input  '$(input.filename)' is now:\n\t$(string(data_block.data)) \n\toption: $(data_block.option)\n")
            end
            setd = true
        end
    end
    if !setd
        setdata!(input, InputData(block_name, option, new_block_data))
        setd = true
    end
    return setd, input
end

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
    setdata!(input::DFInput, data::InputData)

Adds the given data to the input. Should put it in the correct arrays.
"""
function setdata!(input::DFInput, data::InputData)
    setoradd!(input.data, data)
    return input
end

"""
    setdataoption!(input::DFInput, name::Symbol, option::Symbol;; print=true)

Sets the option of specified data.
"""
function setdataoption!(input::DFInput, name::Symbol, option::Symbol; print=true)
    for data in input.data
        if data.name == name
            old_option  = data.option
            data.option = option
            if print info("Option of InputData '$(data.name)' in input '$(input.filename)' set from '$old_option' to '$option'") end
        end
    end
    return input
end

spincalc(input::DFInput) = all(flag(input, :nspin) .!= [nothing, 1])

readoutput(input::DFInput{QE}, filename) = read_qe_output(filename)
readoutput(input::DFInput{Wannier90}, filename) = SymAnyDict()
