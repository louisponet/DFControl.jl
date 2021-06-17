
#----------- Basic Interaction --------------#

"Returns the flag with given symbol."
function Base.getindex(calculation::DFCalculation, n::Symbol)
    if haskey(calculation.flags, n)
        return calculation.flags[n]
    else
        error("Flag:$n Not found in flags of Input $(name(calculation))")
    end
end
"Sets the flag"
Base.setindex!(calculation::DFCalculation, dat, key) = set_flags!(calculation, key => dat)

"""
    set_flags!(calculation::DFCalculation, flags...; print=true)

Sets the specified flags in the calculation.
"""
function set_flags!(calculation::DFCalculation{T}, flags...; print=true) where T
    found_keys = Symbol[]
    for (flag, value) in flags
        flag_type = flagtype(calculation, flag)
        if flag_type !== Nothing
            !(flag in found_keys) && push!(found_keys, flag)
            try
                if isa(value, AbstractVector{<:AbstractVector}) && flag_type <: AbstractVector
                    value = [convert.(eltype(flag_type), v) for v in value]
                else
                    value = convert(flag_type, value)
                end
            catch
                print && ( @warn "Filename '$(name(calculation))':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            old_data = haskey(calculation.flags, flag) ? calculation.flags[flag] : ""
            calculation.flags[flag] = value
            print && (@info "$(name(calculation)):\n  -> $flag:\n      $old_data set to: $value\n")
        else
            print && @warn "Flag $flag was ignored since it could not be found in the allowed flags for calculation $(name(calculation))."
        end
    end
    return found_keys, calculation
end

"""
    rm_flags!(calculation::DFCalculation, flags...)

Remove the specified flags.
"""
function rm_flags!(calculation::DFCalculation, flags...; print=true)
    for flag in flags
        if haskey(calculation.flags, flag)
            pop!(calculation.flags, flag, false)
            print && (@info "Removed flag '$flag' from calculation '$(name(calculation))'")
        end
    end
    return calculation
end

data(calculation::DFCalculation, n) = getfirst(x-> name(x) == n, data(calculation))

"""
    set_data!(calculation::DFCalculation, block_name::Symbol, new_block_data; option=nothing, print=true)

sets the data of the specified 'InputData' to the new data. Optionally also sets the 'InputData' option.
"""
function set_data!(calculation::DFCalculation, block_name::Symbol, new_block_data; option=nothing, print=true)
    setd = false
    for data_block in calculation.data
        if data_block.name == block_name
            if typeof(data_block.data) != typeof(new_block_data)
                if print @warn "Overwritten data of type '$(typeof(data_block.data))' with type '$(typeof(new_block_data))'." end
            end
            old_data        = data_block.data
            data_block.data = new_block_data
            data_block.option = option == nothing ? data_block.option : option
            if print
                @info "Block data '$(data_block.name)' in calculation  '$(name(calculation))' is now:\n\t$(string(data_block.data)) \n\toption: $(data_block.option)\n"
            end
            setd = true
        end
    end
    if !setd
        set_data!(calculation, InputData(block_name, option, new_block_data))
        setd = true
    end
    return setd, calculation
end

"""
    set_dataoption!(calculation::DFCalculation, name::Symbol, option::Symbol; print=true)

Sets the option of specified data.
"""
function set_dataoption!(calculation::DFCalculation, name::Symbol, option::Symbol; print=true)
    for data in calculation.data
        if data.name == name
            old_option  = data.option
            data.option = option
            print && (@info "Option of InputData '$(data.name)' in calculation '$(name(calculation))' set from '$old_option' to '$option'")
        end
    end
    return calculation
end

execs(calculation::DFCalculation, exec::String) =
	filter(x -> occursin(exec, x.exec), calculation.execs)

exec(calculation::DFCalculation, exec::String) =
	getfirst(x -> occursin(exec, x.exec), calculation.execs)

execflags(calculation::DFCalculation, exec::String) =
	[x.exec => x.flags for x in execs(calculation, exec)]

function set_execflags!(calculation::DFCalculation, exec::String, flags...)
	for e in execs(calculation, exec)
		set_flags!(e, flags...)
	end
end

set_execdir!(calculation::DFCalculation, exec, dir) =
	set_execdir!.(execs(calculation, exec), dir)

rmexecflags!(calculation::DFCalculation, exec::String, flags...) =
	rm_flags!.(execs(calculation, exec), flags...)

runcommand(calculation::DFCalculation) = calculation.execs[1]


"Returns the outputdata for the calculation."
function outputdata(calculation::DFCalculation; print=true, overwrite=true)
    if hasoutput(calculation)
        if !overwrite && !isempty(outdata(calculation))
            return outdata(calculation)
        else
            calculation.outdata = readoutput(calculation)
            return calculation.outdata
        end
    end
    print && (@warn "No output data or output file found for calculation: $(name(calculation)).")
    return SymAnyDict()
end
#---------- Extended Interaction ------------#

"""
    readbands(calculation::DFCalculation)

Tries to read bands from the outputfile of the calculation.
"""
function readbands(calculation::DFCalculation)
    hasoutput_assert(calculation)
    to = readoutput(calculation)
    if haskey(to, :bands)
        return to[:bands]
    elseif haskey(to, :bands_up)
        return (up=to[:bands_up], down=to[:bands_down])
    else
        error("No bands found in $(name(calculation)).")
    end
end

"""
    readfermi(calculation::DFCalculation)

Tries to read the fermi level from the outputfile of the calculation.
"""
function readfermi(calculation::DFCalculation)
    hasoutput_assert(calculation)
    to = readoutput(calculation)
    if haskey(to, :fermi)
        return to[:fermi]
    else
        error("No fermi found in $(name(calculation)).")
    end
end

set_kpoints!(::DFCalculation{P}, args...; kwargs...) where {P} =
    @error "set_kpoints! not implemented for package $P."
    
#-------- Generating new DFCalculations ---------- #

function calculation_from_kpoints(template::DFCalculation, newname, kpoints, newflags...)
    newcalc = DFCalculation(template, newname, newflags...)
    set_name!(newcalc, newname)
    set_kpoints!(newcalc, kpoints, print=false)
    return newcalc
end

gencalc_nscf(::DFCalculation{P}, args...) where {P} =
    @error "gencalc_nscf is not implemented for package $P."
    
gencalc_scf(::DFCalculation{P}, args...) where {P} =
    @error "gencalc_scf is not implemented for package $P."
    
gencalc_projwfc(::DFCalculation{P}, args...) where {P} =
    @error "gencalc_projwfc is not implemented for package $P."
    
gencalc_wan(::DFCalculation{P}, args...) where {P} =
    @error "gencalc_wan is not implemented for package $P."
    
gencalc_vcrelax(::DFCalculation{P}, args...) where {P} =
    @error "gencalc_vcrelax is not implemented for package $P."
    
gencalc_bands(::DFCalculation{P}, args...) where {P} =
    @error "gencalc_bands is not implemented for package $P."

isconverged(::DFCalculation{P}) where {P} =
    @error "isconverged is not implemented for package $P."

"""
    set_name!(calculation::DFCalculation, name::AbstractString)

Sets the name of `calculation` to `name` and updates `calculation.infile` and `calculation.outfile` to conform
with the new name.
"""
function set_name!(calculation::DFCalculation, name::AbstractString; print=true)
    calculation.name = name
    calculation.infile = name * splitext(infilename(calculation))[2]
    calculation.outfile = name * splitext(outfilename(calculation))[2]
    print && @info "\ncalculation.name = $name\ncalculation.infile = $(infilename(calculation))\ncalculation.outfile = $(outfilename(calculation))"
    name
end
