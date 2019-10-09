
"""
    documentation(::Type{QE}, searchstring::AbstractString)

Searches through the description of all the known flags that can be set in QuantumEspresso,
returns the flags where the description contains the `searchstring`.

    documentation(::Type{Elk}, flag::Symbol)

Returns the documentation for a given flag.
"""
function documentation(::Type{QE}, searchstring::AbstractString)
    found = Pair{String, Vector{QEFlagInfo}}[]
    for inputinfo in QEInputInfos
        foundflags = QEFlagInfo[]
        for fi in allflags(inputinfo)
            if occursin(searchstring, fi.description)
                push!(foundflags, fi)
            end
        end
        if !isempty(foundflags)
            push!(found, inputinfo.exec => foundflags)
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
    found = Pair{Symbol, Vector{ElkFlagInfo}}[]
    for b in ELK_CONTROLBLOCKS
        found_flags = ElkFlagInfo[]
        for f in b.flags
            if occursin(searchstring, f.description)
                push!(found_flags, f)
            end
        end
        !isempty(found_flags) && push!(found, b.name => found_flags)
    end
    return found
end

function documentation(::Type{Elk}, flagsymbol::Symbol)
    flag = elk_flaginfo(flagsymbol)
    if flag == nothing
        error("No documentation found for flag $flagsymbol, are you sure it is a valid flag for Elk?")
    end
    return flag
end
