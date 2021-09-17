
df_show(args...) = Base.show(args...)

df_show_type(io::IO, ::Type{T}) where {T} = dfprintln(io, crayon"red", "$T:", crayon"reset")
df_show_type(io::IO, x) = df_show_type(io, typeof(x))

function df_show(io::IO, block::InputData)
    s = """Block name: $(block.name)
    Block option: $(block.option)
    Block data:
    """
    dfprintln(io, s)
    dfprintln(io, string(block.data) * "\n\n")
    return
end

function df_show(io::IO, band::Band)
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

function df_show(io::IO, job::Job)
    reset = crayon"reset"
    fieldns = [:name, :version]
    fs = string.(filter(x -> !isempty(x), getfield.((job,), fieldns)))
    fns = string.(fieldns)
    insert!(fns, 2, "dir")
    insert!(fs, 2, Jobs.main_job_dir(job))
    push!(fns, "versions")
    versions = Client.versions(job)
    push!(fs, join(string.(versions), ", "))
    if haskey(job.metadata, :timestamp)
        push!(fns, "last submission")
        push!(fs, string(round(job.metadata[:timestamp], Second)))
    end
    push!(fns, "running")
    is_running = Jobs.main_job_dir(job) != abspath(job) ? false : Client.isrunning(job)
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
    is = job.calculations
    last = !isempty(versions) ? Client.last_running_calculation(job) : -1
    if !isempty(is)
        dfprintln(io, crayon"cyan", line, reset)
        dfprintln(io, reset, "(", crayon"green", "scheduled", reset, ", ", crayon"red",
                  "not scheduled", reset, ")")
        ln = maximum(length.(string.(map(x -> x.name, is))))
        for i in is
            n = i.name
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

function df_show(io::IO, c::Calculation)
    df_show_type(io, c)
    s = """name  = $(c.name)
    dir   = $(c.dir)
    exec = $(c.exec.exec)
    run   = $(c.run)
    data  = $([e.name for e in c.data])
    flags:"""
    dfprintln(io, s)

    if !isempty(c.flags)
        if eltype(c) == QE
            namelist_flags = Dict{Symbol,Dict{Symbol,Any}}()
            info = Calculations.qe_calculation_info(c)
            if info !== nothing
                flag_keys = keys(c.flags)
                for i in info.control
                    for f in i.flags
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
                write_flags(io, c.flags)
            end
        else
            write_flags(io, c.flags)
        end
    end
    dfprint(io, crayon"reset")
    return
end

function df_show(io::IO, flag_info::Calculations.QEFlagInfo{T}) where {T}
    df_show_type(io, flag_info)
    dfprintln(io, crayon"yellow", "name : $(flag_info.name)")
    dfprintln(io, crayon"cyan", "description:", crayon"reset")
    dfprint(io, "\t" * replace(flag_info.description, "\n" => "\n\t"))
    return
end

function df_show(io::IO, flag_info::Calculations.ElkFlagInfo{T}) where {T}
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

function df_show(io::IO, info::Calculations.ElkControlBlockInfo)
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

function df_show(io::IO, el::Element)
    for f in fieldnames(typeof(el))
        dfprintln(io, crayon"red", "$f: ", crayon"reset", "$(getfield(el, f))")
    end
    return dfprint(io, crayon"reset")
end

function df_show(io::IO, str::Structure)
    dfprintln(io, crayon"cyan", "Structure", crayon"reset")
    dfprintln(io, crayon"red", "    cell parameters:")
    dfprint(io, crayon"reset",
            "\t a = $((str.cell[:,1]...,))\n\t b = $((str.cell[:,2]...,))\n\t c = $((str.cell[:,3]...,))\n")
    dfprintln(io, crayon"red", "    nat:", crayon"reset", " $(length(str.atoms))")
    dfprintln(io, crayon"red", "    ntyp:", crayon"reset", " $(length(unique(str.atoms)))")
    show(io, str.atoms)
    return dfprintln(io, crayon"reset")
end

function df_show(io::IO, at::Atom)
    dfprintln(io)
    dfprintln(io, crayon"cyan", "Atom")
    dfprintln(io, crayon"red", "    name: ", crayon"reset", "$(at.name)")
    for f in fieldnames(typeof(at))[3:end-1]
        f == :pseudo && continue
        fld = getfield(at, f)
        if f in (:position_cart, :position_cryst)
            dfprintln(io, crayon"red", "    $f: ", crayon"reset", "$((fld...,))")
        else
            dfprintln(io, crayon"red", "    $f: ", crayon"reset", "$fld")
        end
    end
    dfprintln(io, crayon"red", "    dftu:", crayon"reset")
    df_show(io, at.dftu)
    dfprintln(io, crayon"reset")
    return
end

function df_show(io::IO, proj::Projection)
    println(io, crayon"cyan", "Orbital: ", crayon"reset", "$(proj.orbital.name)")
    println(io, crayon"red", "start index: ", crayon"reset", "$(proj.start)")
    return println(io, crayon"red", "last index: ", crayon"reset", "$(proj.last)")
end
