
#----------- Basic Interaction --------------#
"""
    set_name!(calculation::DFCalculation, name::AbstractString)

Sets `calculation.name`, and `calculation.infile` and `calculation.outfile` to conform
with the new `name`.
"""
function set_name!(calculation::DFCalculation, name::AbstractString; print = true)
    calculation.name = name
    calculation.infile = name * splitext(infilename(calculation))[2]
    calculation.outfile = name * splitext(outfilename(calculation))[2]
    print &&
        @info "\ncalculation.name = $name\ncalculation.infile = $(infilename(calculation))\ncalculation.outfile = $(outfilename(calculation))"
    return name
end

"""
    getindex(calculation::DFCalculation, n::Symbol)

Returns the flag with given symbol.

    getindex(job::DFJob, flag::Symbol)
    
Searches through the job's calculations for the requested flag.
A `Dict` will be returned with calculation names as the keys and the flags as values.
"""
function Base.getindex(calculation::DFCalculation, n::Symbol)
    if haskey(calculation.flags, n)
        return calculation.flags[n]
    else
        error("Flag:$n Not found in flags of Input $(name(calculation))")
    end
end

function Base.getindex(job::DFJob, flg::Symbol)
    outdict = Dict()
    for i in calculations(job)
        tfl = flag(i, flg)
        if tfl != nothing
            outdict[name(i)] = tfl
        end
    end
    return outdict
end

"""
    setindex!(calculation::DFCalculation, value, flag::Symbol)

Sets flags.

    setindex!(job::DFJob, value, flag::Symbol)

Set `flag` in all the appropriate calculations to the `value`.
"""
Base.setindex!(calculation::DFCalculation, dat, key) = set_flags!(calculation, key => dat)

function Base.setindex!(job::DFJob, value, key::Symbol)
    for calculation in calculations(job)
        calculation[key] = value
    end
end

"""
    set_flags!(calculation::DFCalculation, flags::Pair{Symbol, Any}...; print=true)

Sets multiple flags in one go. Flag validity and type are verified.

    set_flags!(job::DFJob, calculations::Vector{<:DFCalculation}, flags::Pair{Symbol,<:Any}...; print=true)
    set_flags!(job::DFJob, flags::Pair{Symbol,<:Any}...; print=true)

Sets the flags in the names to the flags specified.
This only happens if the specified flags are valid for the names.

The values that are supplied will be checked whether they are valid.
"""
function set_flags!(calculation::DFCalculation{T}, flags...; print = true) where {T}
    found_keys = Symbol[]
    for (flag, value) in flags
        flag_type = flagtype(calculation, flag)
        if flag_type !== Nothing
            !(flag in found_keys) && push!(found_keys, flag)
            try
                if isa(value, AbstractVector{<:AbstractVector}) &&
                   flag_type <: AbstractVector
                    value = [convert.(eltype(flag_type), v) for v in value]
                else
                    value = convert(flag_type, value)
                end
            catch
                print &&
                    (@warn "Filename '$(name(calculation))':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            old_data = haskey(calculation.flags, flag) ? calculation.flags[flag] : ""
            calculation.flags[flag] = value
            print &&
                (@info "$(name(calculation)):\n  -> $flag:\n      $old_data set to: $value\n")
        else
            print &&
                @warn "Flag $flag was ignored since it could not be found in the allowed flags for calculation $(name(calculation))."
        end
    end
    return found_keys, calculation
end

function set_flags!(job::DFJob, calculations::Vector{<:DFCalculation}, flags...;
                    print = true)
    found_keys = Symbol[]

    for calc in calculations
        t_, = set_flags!(calc, flags...; print = print)
        push!(found_keys, t_...)
    end
    nfound = setdiff([k for (k, v) in flags], found_keys)
    if print && length(nfound) > 0
        f = length(nfound) == 1 ? "flag" : "flags"
        dfprintln("$f '$(join(":" .* String.(setdiff(flagkeys, found_keys)),", "))' were not found in the allowed calculation variables of the specified calculations!")
    end
    return job
end
function set_flags!(job::DFJob, flags...; kwargs...)
    return set_flags!(job, calculations(job), flags...; kwargs...)
end
function set_flags!(job::DFJob, name::String, flags...; fuzzy = true, kwargs...)
    return set_flags!(job, calculations(job, name, fuzzy), flags...; kwargs...)
end

"""
    rm_flags!(calculation::DFCalculation, flags::Symbol...)
    rm_flags!(job::DFJob, flags::Symbol...)

Remove the specified flags.
"""
function rm_flags!(calculation::DFCalculation, flags::Symbol...; print = true)
    for flag in flags
        if haskey(calculation.flags, flag)
            pop!(calculation.flags, flag, false)
            print && (@info "Removed flag '$flag' from calculation '$(name(calculation))'")
        end
    end
    return calculation
end

function rm_flags!(job::DFJob, flags::Symbol...; kwargs...)
    for c in calculations(job)
        rm_flags!(c, flags...; kwargs...)
    end
end

"""
    data(calculation::DFCalculation)
    data(calculation::DFCalculation, n::Symbol)

The former returns `calculation.data`, the later -- the `InputData` with name `n`.
"""
data(calculation::DFCalculation, n::Symbol) = getfirst(x -> name(x) == n, data(calculation))

"""
    set_data!(calculation::DFCalculation, block_name::Symbol, new_block_data; option::Symbol=nothing, print::Bool=true)

Searches for an `InputData` for which `InputData.name == block_name`, and sets `DFCalculation.data = new_block_data`.
If `option` is specified it is set, i.e. `InputData.option = option`.
"""
function set_data!(calculation::DFCalculation, block_name::Symbol, new_block_data;
                   option::Symbol = nothing, print::Bool = true)
    setd = false
    for data_block in calculation.data
        if data_block.name == block_name
            if typeof(data_block.data) != typeof(new_block_data)
                if print
                    @warn "Overwritten data of type '$(typeof(data_block.data))' with type '$(typeof(new_block_data))'."
                end
            end
            old_data = data_block.data
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
    set_data_option!(calculation::DFCalculation, name::Symbol, option::Symbol; print=true)

Searches for an `InputData` for which `InputData.name == block_name`, and sets `InputData.option = option`.
"""
function set_data_option!(calculation::DFCalculation, name::Symbol, option::Symbol;
                          print = true)
    for data in calculation.data
        if data.name == name
            old_option  = data.option
            data.option = option
            print &&
                (@info "Option of InputData '$(data.name)' in calculation '$(name(calculation))' set from '$old_option' to '$option'")
        end
    end
    return calculation
end

"""
    outputdata(calculation::DFCalculation; extra_parse_funcs=[], print=true, overwrite=true)

If an output file exists for `calculation` this will parse it and return a `Dict` with the parsed data.
If `overwrite=false` and `calculation.outputdata` is not empty, this will be returned instead of reparsing the
output file.

`extra_parse_funcs` should be `Vector{Pair{String,Function}}`, where the string will be used to match `occursin(str, line)`
for each line of the output file. If this returns `true`, the function will be called as `func(results_dict, line, file)`.
The purpose is to allow for additional parsing that is not implemented, or for temporary non-standard values that are printed
while working on the DFT code, e.g. for debugging.

!!! note
    This only works for files where not all information is pulled out for yet,
    e.g. projwfc.x outputs are fully parsed already.

Example (from src/qe/fileio.jl):
```julia

function qe_parse_nat(results, line, f)
    results[:nat] = parse(Int, split(line)[end])
end

outputdata(job["scf"], extra_parse_funcs = ["number of atoms/cell" => qe_parse_nat])
```
"""
function outputdata(calculation::DFCalculation;
                    extra_parse_funcs::Vector{<:Pair{String}} = Pair{String}[],
                    print = true, overwrite = true)
    if hasoutput(calculation)
        if !overwrite && !isempty(outdata(calculation))
            return outdata(calculation)
        else
            t = readoutput(calculation; parse_funcs = extra_parse_funcs)
            calculation.outdata = t === nothing ?
                                  parse_file(outpath(calculation), extra_parse_funcs) : t
            return calculation.outdata
        end
    end
    print &&
        (@warn "No output data or output file found for calculation: $(name(calculation)).")
    return SymAnyDict()
end
#---------- Extended Interaction ------------#

"""
    readbands(calculation::DFCalculation)

Parses the `outputdata` associated with `calculation` and returns the bands if present.

    readbands(job::DFJob)

Will try to find the first bandstructure calculation present in the `job` with a valid
output file, and return the `bands` if present.
If no normal bandstructure calculation is present, it will return the bands produced
by a possible `nscf` calculation if one exists with a valid output file.
"""
function readbands(calculation::DFCalculation)
    hasoutput_assert(calculation)
    to = readoutput(calculation)
    if haskey(to, :bands)
        return to[:bands]
    elseif haskey(to, :bands_up)
        return (up = to[:bands_up], down = to[:bands_down])
    else
        error("No bands found in $(name(calculation)).")
    end
end

"""
    readfermi(calculation::DFCalculation)

Parses the `outputdata` associated with `calculation` and returns the fermi level if present.

    readfermi(job::DFJob)

Finds the first `scf` calculation present in the `job`,
reads its `outputdata` and returns the fermi level if present.
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

"""
    set_kpoints!(calculation::DFCalculation{QE}, k_grid::NTuple{3, Int}; print=true)
    set_kpoints!(calculation::DFCalculation{QE}, k_grid::NTuple{6, Int}; print=true)
    set_kpoints!(calculation::DFCalculation{QE}, k_grid::Vector{<:NTuple{4}}; print=true, k_option=:crystal_b)

Convenience function to set the `:k_points` data block of `calculation`.
The three different methods are targeted at `nscf`, `scf` or `vcrelax`,
and `bands` calculations, respectively.
For the `nscf` version an explicit list of `k_points` will be generated.

    set_kpoints!(calculation::DFCalculation{Wannier90}, k_grid::NTuple{3, Int})

Similar to the `nscf` targeted function in the sense that it will generate
an explicit list of `k_points`, adhering to the same rules as for the `nscf`.
The `mp_grid` flag will also automatically be set.
"""
function set_kpoints!(::DFCalculation{P}, args...; kwargs...) where {P}
    @error "set_kpoints! not implemented for package $P."
end

#-------- Generating new DFCalculations ---------- #

function calculation_from_kpoints(template::DFCalculation, newname, kpoints, newflags...)
    newcalc = DFCalculation(template, newname, newflags...)
    set_name!(newcalc, newname)
    set_kpoints!(newcalc, kpoints; print = false)
    return newcalc
end

function gencalc_nscf(::DFCalculation{P}, args...) where {P}
    @error "gencalc_nscf is not implemented for package $P."
end

function gencalc_scf(::DFCalculation{P}, args...) where {P}
    @error "gencalc_scf is not implemented for package $P."
end

function gencalc_projwfc(::DFCalculation{P}, args...) where {P}
    @error "gencalc_projwfc is not implemented for package $P."
end

function gencalc_wan(::DFCalculation{P}, args...) where {P}
    @error "gencalc_wan is not implemented for package $P."
end

function gencalc_vcrelax(::DFCalculation{P}, args...) where {P}
    @error "gencalc_vcrelax is not implemented for package $P."
end

function gencalc_bands(::DFCalculation{P}, args...) where {P}
    @error "gencalc_bands is not implemented for package $P."
end

function isconverged(::DFCalculation{P}) where {P}
    @error "isconverged is not implemented for package $P."
end
