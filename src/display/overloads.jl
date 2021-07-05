#printing that is not needed in Atom
import Base: show
df_show_type(io::IO, ::Type{T}) where {T} = dfprintln(io, crayon"red", "$T:", crayon"reset")
df_show_type(io::IO, x) = df_show_type(io, typeof(x))

function show(io::IO, block::InputData)
    s = """Block name: $(block.name)
    Block option: $(block.option)
    Block data:
    """
    dfprintln(io, s)
    dfprintln(io, string(block.data) * "\n\n")
    return
end

show(io::IO, data::Vector{InputData}) = map(x -> show(io, x), data)

function show(io::IO, band::DFBand{T}) where {T<:AbstractFloat}
    df_show_type(io, band)
    string = """
    k_points of length $(length(band.k_points_cart)):
    cart:    $(band.k_points_cart[1]) -> $(band.k_points_cart[end])
    """
    !isempty(band.k_points_cryst) &&
        (string *= "cryst:   $(band.k_points_cryst[1]) -> $(band.k_points_cryst[end])\n")
    string *= """eigvals: $(band.eigvals[1]) -> $(band.eigvals[end])
    extra:   $(band.extra)
    """
    dfprintln(io, string)
    return
end

show(io::IO, bands::Vector{<:DFBand}) = map(x -> show(io, x), bands)

function show(io::IO, job::DFJob)
    reset = crayon"reset"
    fieldns = [:name, :version]
    fs = string.(filter(x -> !isempty(x), getfield.((job,), fieldns)))
    fns = string.(fieldns)
    insert!(fns, 2, "local_dir")
    insert!(fs, 2, main_job_dir(job))
    push!(fns, "available versions")
    push!(fs, join(string.(versions(job)), ", "))
    if haskey(job.metadata, :timestamp)
        push!(fns, "last submission")
        push!(fs, string(round(job.metadata[:timestamp], Dates.Second)))
    end
    push!(fns, "running")
    is_running = isrunning(job)
    push!(fs, string(is_running))
    lfns = maximum(length.(fns)) + maximum(length.(fs)) + 4
    line = "+"
    for i in 1:div(lfns, 2)+1
        line *= "-"
    end
    dfprint(io, crayon"cyan", line[1:end-2])
    dfprint(io, crayon"cyan", "DFJOB")
    for i in 1:div(lfns, 2)
        line *= "-"
    end
    line *= "+"
    totlen = length(line)
    dfprintln(io, crayon"cyan", line[div(lfns, 2)+6:end], reset)
    lfns = maximum(length.(string.(fns)))
    for (fn, f) in zip(fns, fs)
        isname = fn == "name"
        l = length(fn)
        fn *= ":"
        for i in 1:lfns-l
            fn *= " "
        end
        dfprint(io, crayon"cyan", "|", reset)
        isname ? dfprint(io, " $fn ", crayon"magenta", "$f", reset) :
        dfprint(io, " $fn $f", reset)
        for i in 1:totlen-(length(fn)+4+length(f))
            dfprint(io, " ")
        end
        dfprintln(io, crayon"cyan", "|", reset)
    end
    is = calculations(job)
    last = last_running_calculation(job)
    if !isempty(is)
        dfprintln(io, crayon"cyan", line, reset)
        dfprintln(io, reset, "(", crayon"green", "scheduled", reset, ", ", crayon"red",
                  "not scheduled", reset, ")")
        ln = maximum(length.(string.(name.(is))))
        for (si, i) in enumerate(is)
            n = name(i)
            cr = i.run ? crayon"green" : crayon"red"
            dfprint(io, cr,
                    i == last ? (is_running ? "\t$n <- running\n" : "\t$n <- ran last\n") :
                    "\t$n\n")
        end
    end
    dfprint(io, reset)
    return
end

function write_flags(io, fls, prestr = "")
    fl = maximum(length.(string.(keys(fls))))
    for (f, v) in fls
        fs = string(f)
        l = length(fs)
        for i in 1:fl-l
            fs *= " "
        end
        dfprintln(io, crayon"cyan", prestr * "\t$fs", crayon"yellow", " => ",
                  crayon"magenta", "$v")
    end
end

function show(io::IO, c::DFCalculation)
    df_show_type(io, c)
    s = """name  = $(c.name)
    dir   = $(c.dir)
    execs = $(join([e.exec for e in c.execs],", "))
    run   = $(c.run)
    data  = $([e.name for e in data(c)])
    flags:"""
    dfprintln(io, s)

    if !isempty(flags(c))
        if package(c) == QE
            namelist_flags = Dict{Symbol,Dict{Symbol,Any}}()
            info = qe_calculation_info(c)
            if info !== nothing
                flag_keys = keys(flags(c))
                for i in info.control
                    for f in flags(i)
                        if f.name ∈ flag_keys
                            if haskey(namelist_flags, i.name)
                                namelist_flags[i.name][f.name] = c[f.name]
                            else
                                namelist_flags[i.name] = Dict{Symbol,Any}()
                                namelist_flags[i.name][f.name] = c[f.name]
                            end
                        end
                    end
                end
                for (inf, flgs) in namelist_flags
                    dfprintln(io, crayon"green", "\t&$inf", crayon"reset")
                    write_flags(io, flgs, "\t\t")
                end
            else
                write_flags(io, flags(c))
            end
        else
            write_flags(io, flags(c))
        end
    end
    dfprint(io, crayon"reset")
    return
end

function show(io::IO, flag_info::QEFlagInfo{T}) where {T}
    df_show_type(io, flag_info)
    dfprintln(io, crayon"yellow", "name : $(flag_info.name)")
    dfprintln(io, crayon"cyan", "description:", crayon"reset")
    dfprint(io, "\t" * replace(flag_info.description, "\n" => "\n\t"))
    return
end

function show(io::IO, flag_info::ElkFlagInfo{T}) where {T}
    df_show_type(io, flag_info)
    dfprintln(io, crayon"yellow", "name : $(flag_info.name)")
    dfprintln(io, crayon"green", "default:", crayon"reset")
    if flag_info.default != nothing
        dfprintln(io, "\t" * "$(flag_info.default)")
    else
        dfprintln(io, "\t" * "nothing")
    end
    dfprintln(io, crayon"cyan", "description:", crayon"reset")
    dfprintln(io, "\t" * replace(flag_info.description, "\n" => "\n\t"))
    return
end

function show(io::IO, info::ElkControlBlockInfo)
    dfprintln(io, crayon"yellow", "name = $(info.name)")
    dfprintln(io, crayon"cyan", "flags:", crayon"reset")
    for flag in info.flags
        dfprintln(io, "\t$(flag.name)")
    end
    dfprintln(io, crayon"cyan", "description:", crayon"reset")
    dfprintln(io, "\t" * replace(info.description, "\n" => "\n\t"))
    dfprint(io, crayon"reset")
    return
end

function show(io::IO, el::Element)
    for f in fieldnames(typeof(el))
        dfprintln(io, crayon"red", "$f: ", crayon"reset", "$(getfield(el, f))")
    end
    return dfprint(io, crayon"reset")
end

function show(io::IO, str::AbstractStructure)
    dfprintln(io, crayon"cyan", "Structure", crayon"reset")
    dfprintln(io, crayon"red", "    cell parameters:")
    dfprint(io, crayon"reset",
            "\t a = $((str.cell[:,1]...,))\n\t b = $((str.cell[:,2]...,))\n\t c = $((str.cell[:,3]...,))\n")
    dfprintln(io, crayon"red", "    nat:", crayon"reset", " $(length(str.atoms))")
    dfprintln(io, crayon"red", "    ntyp:", crayon"reset", " $(length(unique(str.atoms)))")
    for a in atoms(str)
        show(io, a)
    end
    return dfprintln(io, crayon"reset")
end

function show(io::IO, at::AbstractAtom{T,LT}) where {T,LT<:Length{T}}
    dfprintln(io)
    dfprintln(io, crayon"cyan", "Atom")
    dfprintln(io, crayon"red", "    name: ", crayon"reset", "$(name(at))")
    for f in fieldnames(typeof(at))[3:end-1]
        fld = getfield(at, f)
        if f in (:position_cart, :position_cryst)
            dfprintln(io, crayon"red", "    $f: ", crayon"reset", "$((fld...,))")
        else
            dfprintln(io, crayon"red", "    $f: ", crayon"reset", "$fld")
        end
    end
    dfprint(io, crayon"red", "    dftu: ", crayon"reset")
    for f in fieldnames(DFTU)
        val = getfield(dftu(at), f)
        if !isdefault(val)
            dfprint(io, "$f: $val, ")
        end
    end
    dfprintln(io, crayon"reset")
    return
end
