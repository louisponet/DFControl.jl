#these are all the control data, they hold the flags that guide the calculation
mutable struct InputData
    name   ::Symbol
    option ::Symbol
    data   ::Any
end
name(data::InputData) = data.name

@with_kw mutable struct DFInput{P <: Package}
    name     ::String
    dir      ::String
    flags    ::SymAnyDict
    data     ::Vector{InputData}
    execs    ::Vector{Exec}
    run      ::Bool
    outdata  ::SymAnyDict=SymAnyDict()
end
DFInput{P}(name, dir, flags, data, execs, run) where P<:Package = DFInput{P}(name, dir, flags, data, execs, run, SymAnyDict())

"""
    DFInput(template::DFInput, name, newflags...; runcommand=template.runcommand, run=true)

Creates a new `DFInput` from a template `DFInput`, setting the newflags in the new one.
"""
function DFInput(template::DFInput, name, newflags...; excs=execs(template), run=true, data=nothing, dir=template.dir)
    newflags = Dict(newflags...)

    input          = deepcopy(template)
    input.name     = name
    input.execs    = excs
    input.run      = run
    input.dir      = dir
    setflags!(input, newflags..., print=false)

    if data != nothing
        for (name, (option, data)) in data
            setdata!(input, name, data, option=option, print=false)
        end
    end
    return input
end

name(input::DFInput) = input.name
dir(input::DFInput)  = input.dir
setdir!(input::DFInput, dir) = (input.dir = dir)
namewext(input::DFInput, ext)      = name(input) * ext
infile(input::DFInput{QE})         = namewext(input, ".in")
infile(input::DFInput{Wannier90})  = namewext(input, ".win")
outfile(input::DFInput{QE})        = namewext(input, ".out")
outfile(input::DFInput{Wannier90}) = namewext(input, ".wout")
inpath(input::DFInput)             = joinpath(dir(input), infile(input))
outpath(input::DFInput)            = joinpath(dir(input), outfile(input))


flags(input::DFInput)    = input.flags
function flag(input::DFInput, flag::Symbol)
    if haskey(input.flags, flag)
        return input.flags[flag]
    end
end

function Base.getindex(input::DFInput, n::Symbol)
    if haskey(input.flags, n)
        return input.flags[n]
    else
        tmp = getfirst(x->name(x) == n, data(input))
        if tmp != nothing
            return tmp
        else
            error("Id :$n \n\t Not found in flags or data of Input $(name(input))")
        end
    end
end


Base.eltype(::DFInput{P}) where P = P
package(::DFInput{P}) where P = P

data(input::DFInput)     = input.data
data(input::DFInput, name) = getfirst(x-> x.name == name, input.data)
data(input::Vector{InputData}, name) = getfirst(x-> x.name == name, input.data)

exec(input::DFInput, exec::String) = getfirst(x -> occursin(exec, x.exec), input.execs)
execs(input::DFInput) = input.execs
execs(input::DFInput, exec::String) = filter(x -> occursin(exec, x.exec), input.execs)

execflags(input::DFInput, exec::String) = [x.flags for x in execs(input, exec)]
setexecflags!(input::DFInput, exec::String, flags...) = setflags!.(execs(input, exec), (flags,)...)
setexecdir!(input::DFInput, exec, dir) = setexecdir!.(execs(input, exec), dir)
rmexecflags!(input::DFInput, exec::String, flags...) = rmflags!.(execs(input, exec), flags...)

runcommand(input::DFInput) = input.execs[1]

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
    print && calc != "'nscf'" && (@warn "Expected calculation to be 'nscf'.\nGot $calc.")
    setdata!(input, :k_points, kgrid(k_grid..., :nscf), option = :crystal, print=print)
    prod(k_grid) > 100 && setflags!(input, :verbosity => "'high'", print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::NTuple{6, Int}; print=true) #scf
    calc = flag(input, :calculation)
    print && (calc != "'scf'" || !occursin("relax", calc)) && (@warn "Expected calculation to be 'scf', 'vc-relax', 'relax'.\nGot $calc.")
    setdata!(input, :k_points, [k_grid...], option = :automatic, print=print)
    prod(k_grid[1:3]) > 100 && setflags!(input, :verbosity => "'high'", print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::Vector{NTuple{4, T}}; print=true, k_option=:crystal_b) where T<:AbstractFloat
    calc = flag(input, :calculation)
    print && calc != "'bands'" && (@warn "Expected calculation to be 'bands', got $calc.")
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
            @info "Verbosity is set to high because num_kpoints > 100,\n
                       otherwise bands won't get printed."
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
        if flag_type != Nothing
            !(flag in found_keys) && push!(found_keys, flag)
            try
                value = convertflag(flag_type, value)
            catch
                print && ( @warn "Filename '$(name(input))':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            old_data = haskey(input.flags, flag) ? input.flags[flag] : ""
            input.flags[flag] = value
            print && (@info "$(name(input)):\n  -> $flag:\n      $old_data set to: $value\n")
        end
    end
    return found_keys, input
end

Base.setindex!(input::DFInput, dat, key) = setflags!(input, key => dat; print=false)
function Base.setindex!(input::DFInput, dat::InputData, id)
    index = findfirst(x -> name(x)==id, data(input))
    data(input)[index] = dat
end

"Runs through all the set flags and checks if they are allowed and set to the correct value"
function cleanflags!(input::DFInput)
    for (flag, value) in flags(input)
        flagtype_ = flagtype(input, flag)
        if flagtype_ == Nothing
            @warn "Flag $flag was not found in allowed flags for exec $(execs(input)[2]). Removing flag."
            rmflags!(input, flag)
            continue
        end
        if !(eltype(value) <: flagtype_)
            try
                flags(input)[flag] = convertflag(flagtype_, value)
            catch
                error("Input $(name(input)): Could not convert :$flag of value $value to the correct type ($flagtype_), please set it to the correct type.")
            end
        end
    end
end

#TODO implement abinit and wannier90
"""
    sanitizeflags!(input::DFInput)

Cleans up flags, i.e. remove flags that are not allowed and convert all
flags to the correct types.
Tries to correct common errors for different input types.
"""
function sanitizeflags!(input::DFInput)
    cleanflags!(input)
end
function sanitizeflags!(input::DFInput{QE})
    cleanflags!(input)
    flag(input, :outdir) != fortstring(dir(input)) && setflags!(input, :outdir => fortstring(dir(input)), print=false)
    #TODO add all the required flags
end

"""
    rmflags!(input::DFInput, flags...)

Remove the specified flags.
"""
function rmflags!(input::DFInput, flags...; print=true)
    for flag in flags
        if haskey(input.flags, flag)
            pop!(input.flags, flag, false)
            print && (@info "Removed flag '$flag' from input '$(name(input))'")
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
                if print @warn "Overwritten data of type '$(typeof(data_block.data))' with type '$(typeof(new_block_data))'." end
            end
            old_data        = data_block.data
            data_block.data = new_block_data
            data_block.option = option == nothing ? data_block.option : option
            if print
                @info "Block data '$(data_block.name)' in input  '$(name(input))' is now:\n\t$(string(data_block.data)) \n\toption: $(data_block.option)\n"
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
    setdataoption!(input::DFInput, name::Symbol, option::Symbol; print=true)

Sets the option of specified data.
"""
function setdataoption!(input::DFInput, name::Symbol, option::Symbol; print=true)
    for data in input.data
        if data.name == name
            old_option  = data.option
            data.option = option
            print && (@info "Option of InputData '$(data.name)' in input '$(name(input))' set from '$old_option' to '$option'")
        end
    end
    return input
end

isbandscalc(input::DFInput{QE}) = flag(input, :calculation) == "'bands'"
isnscfcalc(input::DFInput{QE}) = flag(input, :calculation) == "'nscf'"
isscfcalc(input::DFInput{QE}) = flag(input, :calculation) == "'scf'"
isspincalc(input::DFInput{QE}) = all(flag(input, :nspin) .!= [nothing, 1])

outdata(input::DFInput) = input.outdata
hasoutput(input::DFInput) = !isempty(outdata(input))

hasoutfile(input::DFInput) = ispath(outpath(input))
hasnewout(input::DFInput, time) = mtime(outpath(input)) > time

"Returns the outputdata for the input."
function outputdata(input::DFInput; print=true, overwrite=true)
    if hasoutput(input) && !overwrite
        return outdata(input)
    end
    if hasoutfile(input)
        input.outdata = readoutput(input)
        return input.outdata
    end
    print && (@warn "No output data or output file found for input: $(name(input)).")
    return SymAnyDict()
end

readoutput(input::DFInput{QE}) = read_qe_output(outpath(input))
readoutput(input::DFInput{Wannier90}) = SymAnyDict()

function readbands(input::DFInput)
    to = readoutput(input)
    if haskey(to, :bands)
        return to[:bands]
    else
        error("No bands found in $(name(input)).")
    end
end
