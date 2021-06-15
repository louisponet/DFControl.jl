
#----------- Basic Interaction --------------#

"Returns the flag with given symbol."
function Base.getindex(input::DFInput, n::Symbol)
    if haskey(input.flags, n)
        return input.flags[n]
    else
        error("Flag:$n Not found in flags of Input $(name(input))")
    end
end
"Sets the flag"
Base.setindex!(input::DFInput, dat, key) = set_flags!(input, key => dat)

"""
    set_flags!(input::DFInput, flags...; print=true)

Sets the specified flags in the input.
"""
function set_flags!(input::DFInput{T}, flags...; print=true) where T
    found_keys = Symbol[]
    for (flag, value) in flags
        flag_type = flagtype(input, flag)
        if flag_type !== Nothing
            !(flag in found_keys) && push!(found_keys, flag)
            try
                if isa(value, AbstractVector{<:AbstractVector}) && flag_type <: AbstractVector
                    value = [convert.(eltype(flag_type), v) for v in value]
                else
                    value = convert(flag_type, value)
                end
            catch
                print && ( @warn "Filename '$(name(input))':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            old_data = haskey(input.flags, flag) ? input.flags[flag] : ""
            input.flags[flag] = value
            print && (@info "$(name(input)):\n  -> $flag:\n      $old_data set to: $value\n")
        else
            print && @warn "Flag $flag was ignored since it could not be found in the allowed flags for input $(name(input))."
        end
    end
    return found_keys, input
end

"""
    rm_flags!(input::DFInput, flags...)

Remove the specified flags.
"""
function rm_flags!(input::DFInput, flags...; print=true)
    for flag in flags
        if haskey(input.flags, flag)
            pop!(input.flags, flag, false)
            print && (@info "Removed flag '$flag' from input '$(name(input))'")
        end
    end
    return input
end

data(input::DFInput, n) = getfirst(x-> name(x) == n, data(input))

"""
    set_data!(input::DFInput, block_name::Symbol, new_block_data; option=nothing, print=true)

sets the data of the specified 'InputData' to the new data. Optionally also sets the 'InputData' option.
"""
function set_data!(input::DFInput, block_name::Symbol, new_block_data; option=nothing, print=true)
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
        set_data!(input, InputData(block_name, option, new_block_data))
        setd = true
    end
    return setd, input
end

"""
    set_dataoption!(input::DFInput, name::Symbol, option::Symbol; print=true)

Sets the option of specified data.
"""
function set_dataoption!(input::DFInput, name::Symbol, option::Symbol; print=true)
    for data in input.data
        if data.name == name
            old_option  = data.option
            data.option = option
            print && (@info "Option of InputData '$(data.name)' in input '$(name(input))' set from '$old_option' to '$option'")
        end
    end
    return input
end

execs(input::DFInput, exec::String) =
	filter(x -> occursin(exec, x.exec), input.execs)

exec(input::DFInput, exec::String) =
	getfirst(x -> occursin(exec, x.exec), input.execs)

execflags(input::DFInput, exec::String) =
	[x.exec => x.flags for x in execs(input, exec)]

function set_execflags!(input::DFInput, exec::String, flags...)
	for e in execs(input, exec)
		set_flags!(e, flags...)
	end
end

set_execdir!(input::DFInput, exec, dir) =
	set_execdir!.(execs(input, exec), dir)

rmexecflags!(input::DFInput, exec::String, flags...) =
	rm_flags!.(execs(input, exec), flags...)

runcommand(input::DFInput) = input.execs[1]


"Returns the outputdata for the input."
function outputdata(input::DFInput; print=true, overwrite=true)
    if hasoutput(input)
        if !overwrite && !isempty(outdata(input))
            return outdata(input)
        else
            input.outdata = readoutput(input)
            return input.outdata
        end
    end
    print && (@warn "No output data or output file found for input: $(name(input)).")
    return SymAnyDict()
end
#---------- Extended Interaction ------------#

"""
    readbands(input::DFInput)

Tries to read bands from the outputfile of the input.
"""
function readbands(input::DFInput)
    hasoutput_assert(input)
    to = readoutput(input)
    if haskey(to, :bands)
        return to[:bands]
    elseif haskey(to, :bands_up)
        return (up=to[:bands_up], down=to[:bands_down])
    else
        error("No bands found in $(name(input)).")
    end
end

"""
    readfermi(input::DFInput)

Tries to read the fermi level from the outputfile of the input.
"""
function readfermi(input::DFInput)
    hasoutput_assert(input)
    to = readoutput(input)
    if haskey(to, :fermi)
        return to[:fermi]
    else
        error("No fermi found in $(name(input)).")
    end
end

set_kpoints!(::DFInput{P}, args...; kwargs...) where {P} =
    @error "set_kpoints! not implemented for package $P."
    
#-------- Generating new DFInputs ---------- #

function input_from_kpoints(template::DFInput, newname, kpoints, newflags...)
    newcalc = DFInput(template, newname, newflags...)
    set_name!(newcalc, newname)
    set_kpoints!(newcalc, kpoints, print=false)
    return newcalc
end

gencalc_nscf(::DFInput{P}, args...) where {P} =
    @error "gencalc_nscf is not implemented for package $P."
    
gencalc_scf(::DFInput{P}, args...) where {P} =
    @error "gencalc_scf is not implemented for package $P."
    
gencalc_projwfc(::DFInput{P}, args...) where {P} =
    @error "gencalc_projwfc is not implemented for package $P."
    
gencalc_wan(::DFInput{P}, args...) where {P} =
    @error "gencalc_wan is not implemented for package $P."
    
gencalc_vcrelax(::DFInput{P}, args...) where {P} =
    @error "gencalc_vcrelax is not implemented for package $P."
    
gencalc_bands(::DFInput{P}, args...) where {P} =
    @error "gencalc_bands is not implemented for package $P."

isconverged(::DFInput{P}) where {P} =
    @error "isconverged is not implemented for package $P."

"""
    set_name!(input::DFInput, name::AbstractString)

Sets the name of `input` to `name` and updates `input.infile` and `input.outfile` to conform
with the new name.
"""
function set_name!(input::DFInput, name::AbstractString; print=true)
    input.name = name
    input.infile = name * splitext(infilename(input))[2]
    input.outfile = name * splitext(outfilename(input))[2]
    print && @info "\ninput.name = $name\ninput.infile = $(infilename(input))\ninput.outfile = $(outfilename(input))"
    name
end
