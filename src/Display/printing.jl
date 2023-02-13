show_type(io::IO, ::MIME"text/plain", ::Type{T}) where {T} = println(io, crayon"red", "$T:", crayon"reset")
show_type(io::IO, ::MIME"text/plain", x) = show_type(io, typeof(x))

function Base.show(io::IO, ::MIME"text/plain", block::InputData)
    s = """Block name: $(block.name)
    Block option: $(block.option)
    Block data:
    """
    println(io, s)
    println(io, string(block.data) * "\n\n")
    return
end

function Base.show(io::IO, ::MIME"text/plain", band::Band)
    show_type(io, band)
    string = """
    k_points of length $(length(band.k_points_cart)):
    cart:    $(band.k_points_cart[1]) -> $(band.k_points_cart[end])
    """
    !isempty(band.k_points_cryst) &&
        (string *= "cryst:   $(band.k_points_cryst[1]) -> $(band.k_points_cryst[end])\n")
    string *= """eigvals: $(band.eigvals[1]) -> $(band.eigvals[end])
    extra:   $(band.extra)
    """
    println(io, string)
    return
end

Base.show(io::IO, job::Job) = print(io, "Job(name: $(job.name), dir: $(job.dir), server: $(job.server), calcs: $(length(job.calculations)); $(map(x->x.name, job.calculations)))")

function Base.show(io::IO, ::MIME"text/plain", job::Job)
    server = Server(job.server)
    reset = crayon"reset"
    fieldns = [:name, :dir, :version]
    fs = string.(getfield.((job,), fieldns))
    fs[2] = Jobs.main_job_dir(job)
    fns = string.(fieldns)
    push!(fns, "server")
    if isalive(server) 
        push!(fs, "$(job.server), alive")
    else
        push!(fs, "$(job.server), not alive")
    end
        
    
    push!(fns, "versions")
    if isalive(server)
        versions = Client.versions(job)
        push!(fs, join(string.(versions), ", "))
        timestamp = Client.submission_time(job)
        if timestamp != 0
            push!(fns, "last submission")
            push!(fs, string(round(Dates.unix2datetime(timestamp), Second)))
        end
        push!(fns, "state")
        state = Client.state(job)
        push!(fs, string(state))
    else
        push!(fs, "unknown")
        push!(fns, "state")
        state = Jobs.Unknown
        push!(fs, string(Jobs.Unknown))
    end
        
    lfns = maximum(length.(fns)) + maximum(length.(fs)) + 4
    line = "+"
    for i in 1:div(lfns, 2) + 1
        line *= "-"
    end
    print(io, crayon"cyan", line[1:end-2])
    print(io, crayon"cyan", "JOB")
    for i in 1:div(lfns, 2)
        line *= "-"
    end
    line *= "+"
    totlen = length(line)
    println(io, crayon"cyan", line[div(lfns, 2)+4:end], reset)
    lfns = maximum(length.(string.(fns)))
    for (fn, f) in zip(fns, fs)
        isname = fn == "name"
        l = length(fn)
        fn *= ":"
        for i in 1:lfns-l
            fn *= " "
        end
        print(io, crayon"cyan", "|", reset)
        isname ? print(io, " $fn ", crayon"magenta", "$f", reset) :
        print(io, " $fn $f", reset)
        for i in 1:totlen-(length(fn)+4+length(f))
            print(io, " ")
        end
        println(io, crayon"cyan", "|", reset)
    end
    is = job.calculations
    last = isalive(server) && ispath(server, joinpath(Jobs.main_job_dir(job), "job.sh")) ? Client.last_running_calculation(job) : -1
    if !isempty(is)
        println(io, crayon"cyan", line, reset)
        println(io, reset, "(", crayon"green", "scheduled", reset, ", ", crayon"red",
                  "not scheduled", reset, ")")
        ln = maximum(length.(string.(map(x -> x.name, is))))
        for i in is
            n = i.name
            cr = i.run ? crayon"green" : crayon"red"
            if i == last
                if state == Jobs.Running
                    print(io, cr, "\t$n <- running\n")
                elseif state == Jobs.Completed || state == Jobs.Failed || state == Jobs.Cancelled
                    print(io, cr, "\t$n <- ran last\n")
                else
                    print(io, cr, "\t$n\n")
                end
            else
                print(io, cr, "\t$n\n")
            end
        end
    end
    print(io, reset)
    return
end

function write_flags(io, fls, prestr = "")
    flags = String[]
    vals  = String[]
    for (f, v) in fls
        if v isa AbstractArray
            for i in findall(!iszero, v)
                push!(flags, string(f) * "($(join(Tuple(i), ",")))")
                push!(vals, string(v[i]))
            end
        else
            push!(flags, string(f))
            push!(vals, string(v))
        end
    end
    fl = maximum(length, flags)
    sp = sortperm(flags, by=length)
    for s in sp
        fs, v = flags[s], vals[s]
        l = length(fs)
        for i in 1:fl-l
            fs *= " "
        end
        println(io, crayon"cyan", prestr * "  $fs", crayon"yellow", " => ",
                  crayon"magenta", v, crayon"reset")
    end
    println(io,"")
end

function Base.show(io::IO, c::Calculation)
print(io, typeof(c), "(name: $(c.name))")
end

function Base.show(io::IO, ::MIME"text/plain", c::Calculation)
    show_type(io, c)
    s = """name  = $(c.name)
    exec = $(c.exec.exec)
    run   = $(c.run)
    data  = $([e.name for e in c.data])
    flags:"""
    println(io, s)

    if !isempty(c.flags)
        write_flags_separately = false
        for (b, bflags) in c.flags
            if bflags isa Dict
                if !isempty(bflags)
                    println(io, crayon"green", "  &$b", crayon"reset")
                    write_flags(io, bflags, "  ")
                    
                end
            else
                write_flags_separately = true
            end 
        end
        if write_flags_separately
            write_flags(io, c.flags)
        end
    end
    print(io, crayon"reset")
    return
end

function Base.show(io::IO, ::MIME"text/plain", flag_info::Calculations.QEFlagInfo{T}) where {T}
    show_type(io, flag_info)
    println(io, crayon"yellow", "name : $(flag_info.name)")
    println(io, crayon"cyan", "description:", crayon"reset")
    print(io, "\t" * replace(flag_info.description, "\n" => "\n\t"))
    return
end

function Base.show(io::IO, ::MIME"text/plain", flag_info::Calculations.ElkFlagInfo{T}) where {T}
    show_type(io, flag_info)
    println(io, crayon"yellow", "name : $(flag_info.name)")
    println(io, crayon"green", "default:", crayon"reset")
    if flag_info.default != nothing
        println(io, "\t" * "$(flag_info.default)")
    else
        println(io, "\t" * "nothing")
    end
    println(io, crayon"cyan", "description:", crayon"reset")
    println(io, "\t" * replace(flag_info.description, "\n" => "\n\t"))
    return
end

function Base.show(io::IO, ::MIME"text/plain", info::Calculations.ElkControlBlockInfo)
    println(io, crayon"yellow", "name = $(info.name)")
    println(io, crayon"cyan", "flags:", crayon"reset")
    for flag in info.flags
        println(io, "\t$(flag.name)")
    end
    println(io, crayon"cyan", "description:", crayon"reset")
    println(io, "\t" * replace(info.description, "\n" => "\n\t"))
    print(io, crayon"reset")
    return
end

function Base.show(io::IO, ::MIME"text/plain", el::Element)
    for f in fieldnames(typeof(el))
        println(io, crayon"red", "$f: ", crayon"reset", "$(getfield(el, f))")
    end
    return print(io, crayon"reset")
end

function Base.show(io::IO, str::Structure)
    names = map(x->x.name, unique(str.atoms))
    n = map(x->length(filter(y->y.name == x, str.atoms)), names)
    ats = ""
    for (num, name) in zip(n, names)
        ats *= " $num $name"
    end
    print(io, "Structure(volume: $(Structures.volume(str)), atoms:$ats)")
end
function Base.show(io::IO, ::MIME"text/plain", str::Structure)
    println(io, crayon"cyan", "Structure", crayon"reset")
    
    println(io, crayon"red", "    volume: ", crayon"reset", Structures.volume(str))
    println(io, crayon"red", "    cell parameters:")
    print(io, crayon"reset",
            "\t a = $((str.cell[:,1]...,))\n\t b = $((str.cell[:,2]...,))\n\t c = $((str.cell[:,3]...,))\n")
    println(io, crayon"red", "    nat:", crayon"reset", " $(length(str.atoms))")
    println(io, crayon"red", "    ntyp:", crayon"reset", " $(length(unique(str.atoms)))")
    println(io, crayon"red", "    Atoms:", crayon"reset")
    for a in str.atoms
        print("\t")
        show(io, a)
        println(io)
    end
    # show(io, MIME"text/plain"(), str.atoms)
    return print(io, crayon"reset")
end

Base.show(io::IO, at::Atom) = print(io, crayon"cyan", string(at.name), crayon"reset",  " " * "$(at.position_cryst)")

function Base.show(io::IO, ::MIME"text/plain",  at::Atom)
    println(io)
    println(io, crayon"cyan", "Atom")
    println(io, crayon"red", "    name: ", crayon"reset", "$(at.name)")
    for f in fieldnames(typeof(at))[3:end-1]
        f == :pseudo && continue
        fld = getfield(at, f)
        if f in (:position_cart, :position_cryst)
            println(io, crayon"red", "    $f: ", crayon"reset", "$((fld...,))")
        else
            println(io, crayon"red", "    $f: ", crayon"reset", "$fld")
        end
    end
    print(io, crayon"red", "    dftu: ", crayon"reset")
    show(io, at.dftu)
    println(io, crayon"reset")
    return
end

function Base.show(io::IO, ::MIME"text/plain", proj::Projection)
    println(io, crayon"cyan", "Orbital: ", crayon"reset", "$(proj.orbital.name)")
    println(io, crayon"red", "start index: ", crayon"reset", "$(proj.start)")
    return println(io, crayon"red", "last index: ", crayon"reset", "$(proj.last)")
end

Base.show(io::IO, dftu::DFTU) = print(io, "l=$(dftu.l), U=$(dftu.U), J0=$(dftu.J0), α=$(dftu.α), β=$(dftu.β), J=$(dftu.J)")

