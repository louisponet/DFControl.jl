
"""
    documentation(::Type{QE}, searchstring::AbstractString)

Searches through the description of all the known flags that can be set in QuantumEspresso,
returns the flags where the description contains the `searchstring`.

    documentation(::Type{Elk}, flag::Symbol)

Returns the documentation for a given flag.
"""
function documentation(::Type{QE}, searchstring::AbstractString)
    found = Pair{String,Vector{QEFlagInfo}}[]
    for calculationinfo in QECalculationInfos
        foundflags = QEFlagInfo[]
        for fi in allflags(calculationinfo)
            if occursin(searchstring, fi.description)
                push!(foundflags, fi)
            end
        end
        if !isempty(foundflags)
            push!(found, calculationinfo.exec => foundflags)
        end
    end
    return found
end

function documentation(::Type{QE}, flagsymbol::Symbol)
    flag = qe_flaginfo(flagsymbol)
    if flag == nothing
        error("No documentation found for flag $flagsymbol, are you sure it is a valid flag for any QE executables?")
    end
    return flag
end

"""
    documentation(::Type{Elk}, searchstring::AbstractString)

Searches through the description of all the known flags that can be set in Elk,
returns the flags where the description contains the `searchstring`.

    documentation(::Type{Elk}, flag::Symbol)

Returns the documentation for a given flag.
"""
function documentation(::Type{Elk}, searchstring::AbstractString)
    found = Pair{ElkControlBlockInfo,Vector{ElkFlagInfo}}[]
    for b in ELK_CONTROLBLOCKS
        found_flags = ElkFlagInfo[]
        for f in b.flags
            if occursin(searchstring, f.description)
                push!(found_flags, f)
            end
        end
        if !isempty(found_flags) || occursin(searchstring, b.description)
            push!(found, b => found_flags)
        end
    end
    return found
end

function documentation(::Type{Elk}, flagsymbol::Symbol)
    flag = elk_flaginfo(flagsymbol)
    if flag == nothing
        block = getfirst(x -> x.name == flagsymbol, ELK_CONTROLBLOCKS)
        if block == nothing
            error("No documentation found for flag $flagsymbol, are you sure it is a valid flag for Elk?")
        else
            return block
        end
    else
        return flag
    end
end
