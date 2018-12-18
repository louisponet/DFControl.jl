
"""
    search_documentation(::Type{QE}, searchstring)

Searches through the description of all the known flags that can be set in QuantumEspresso,
returns the flags where the description contains the `searchstring`.
"""
function search_documentation(::Type{QE}, searchstring)
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

"""
    documentation(::Type{QE}, flagsymbol::Symbol)

Returns the documentation for a given flag.
"""
documentation(::Type{QE}, flagsymbol::Symbol) = qe_flaginfo(flagsymbol)
