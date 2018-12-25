using DFControl: ABINIT, SymAnyDict, DFInput, InputData, separate!

#REVIEW: Not sure if 0.1 eV 0.1 eV 0.1 eV == 0.1 0.1 0.1 eV
function abi_convert_flagvalue(flag::String, value::String)
    flagsym = Symbol(flag)
    FT = flagtype(ABINIT, flagsym)
    @assert FT != Nothing "Flag $flagsym was not found in Abinit's input flags. If it should exist (according to the Abinit Documentation), please file an issue."
    flagvals = split(replace(lowercase(value), "d"=>"e"))

    conversion = FT(1.0)
    if length(flagvals) > 1
        units, flagvals = separate!(x -> x ∈ ABI_UNIT_NAMES, flagvals)
        if !isempty(units)
            #REVIEW: Only one unit allowed
            conversion = FT(abi_conversions[Symbol(units[end])])
        end
    end
    if length(flagvals) > 1
        return flagsym, conversion .* parse.(FT, flagvals)
    else
        return flagsym, conversion * parse(FT, flagvals[1])
    end
end

function abi_extract_kpoints!(flagdict::Dict)
    if !haskey(flagdict, :kptopt)
        return
    end
    kptopt = pop!(flagdict, :kptopt)
    option = abi_kptopt_to_option(kptopt)
    req_flags = abi_required_flags(option)

    # @assert all(haskey.((flagdict,), req_flags)) "Not all required flags for kptopt $kptopt were found. \n required flags are: $req_flags"

    kflags, flagdict = separate!(x -> first(x) ∈ req_flags, flagdict)
    return InputData(:kpoints, option, kflags)
end

function abi_read_input(filename; execs=[Exec("abinit")], run=true, structure_name="noname")
    t_ = abilab_opt[:abiopen](filename)
    datasets = [filter(x -> first(x) ∉ ABI_STRUCTURE_FLAGS, t) for t in t_[:datasets]]
    structure = Structure(t_[:structure], name=structure_name)

    flagdicts = [SymAnyDict() for i=1:length(datasets)]
    for (data, fdict) in zip(datasets, flagdicts)
        for (flag, val) in data
            flagsym, flagval = abi_convert_flagvalue(flag, val)
            fdict[flagsym] = flagval
        end
    end

    t_keys = keys(flagdicts[1])
    #TODO: Handle  more input data
    kpointdata = [abi_extract_kpoints!(fl) for fl in flagdicts]

    dir, file = splitdir(splitext(filename)[1])

    abi_inputs = DFInput{ABINIT}[]
    for (i, flags) in enumerate(flagdicts)
        push!(abi_inputs, DFInput{ABINIT}(file*"$i", dir, flags, [kpointdata[i]], execs, run))
    end

    return abi_inputs, structure
end
