#printing that is not needed in Atom
import Base: show
df_show_type(io::IO, ::Type{T}) where T = dfprintln(io, crayon"red", "$T:", crayon"reset")
df_show_type(io::IO, x) = df_show_type(io, typeof(x))


function show(io::IO, block::InputData)
    s = """Block name: $(block.name)
    Block option: $(block.option)
    Block data:
    """
    dfprintln(io, s)
    dfprintln(io, string(block.data) * "\n\n")
end

show(io::IO, data::Vector{InputData}) = map(x-> show(io, x), data)

function show(io::IO, band::DFBand{T}) where T <: AbstractFloat
    df_show_type(io, band)
    string = """
    k_points of length $(length(band.k_points_cryst)):
    cart:    $(band.k_points_cart[1]) -> $(band.k_points_cart[end])
    cryst:   $(band.k_points_cryst[1]) -> $(band.k_points_cryst[end])
    eigvals: $(band.eigvals[1]) -> $(band.eigvals[end])
    extra:   $(band.extra)
    """
    dfprintln(io, string)
end

show(io::IO, bands::Vector{<:DFBand}) = map(x->show(io,x),bands)

function show(io::IO, job::DFJob)
    reset = crayon"reset"
    fieldns  = [:name, :local_dir, :server, :server_dir]
    fs   = string.(filter(x->!isempty(x), getfield.((job,), fieldns)))
    fns = string.(fieldns)
    lfns = maximum(length.(fns)) + maximum(length.(fs)) + 4
    line = "+"
    for i=1:div(lfns,2) + 1
        line *= "-"
    end
    dfprint(io, crayon"cyan", line[1:end-2])
    dfprint(io, crayon"cyan", "DFJOB")
    for i=1:div(lfns,2)
        line *= "-"
    end
    line *= "+"
    totlen = length(line)
    dfprintln(io, crayon"cyan", line[div(lfns, 2)+6:end],  reset)
    lfns = maximum(length.(string.(fns)))
    for (fn, f) in zip(fns, fs)
        isname = fn == "name"
        l = length(fn)
        fn *= ":"
        for i=1:lfns-l
            fn *= " "
        end
        dfprint(io, crayon"cyan", "|", reset)
        isname ? dfprint(io, " $fn ", crayon"magenta", "$f", reset) : dfprint(io, " $fn $f", reset)
        for i=1:totlen - (length(fn) + 4 + length(f))
            dfprint(io, " ")
        end
        dfprintln(io, crayon"cyan", "|", reset)
    end
    is = inputs(job)
    dfprintln(io, crayon"cyan", line, reset)
    dfprintln(io, reset,"(", crayon"green", "scheduled", reset, ", ", crayon"red", "not scheduled", reset, ")")
    ln = maximum(length.(string.(name.(is))))
    for (si, i) in enumerate(is)
        n = name(i)
        l = length(n)
        for j=1:ln-l
            n *= " "
        end
        cr = i.run ? crayon"green" : crayon"red"
        dfprint(io, cr, "\t\t$n\n")
    end
    dfprint(io, reset)
end

show(io::IO, at::AbstractAtom) = dfprintln(io, "$(id(at)): $(position(at)[1]) $(position(at)[2]) $(position(at)[3])")

function show(io::IO, in::DFInput)
    df_show_type(io, in)
    s = """name  = $(in.name)
    dir   = $(in.dir)
    execs = $(join([e.exec for e in in.execs],", "))
    run   = $(in.run)
    data  = $([e.name for e in data(in)])
    $(crayon"cyan")flags$(crayon"reset"):"""
    dfprintln(io, s)
    fl = maximum(length.(string.(keys(flags(in)))))
    for (f, v) in flags(in)
        fs = string(f)
        l = length(fs)
        for i=1:fl - l
            fs *= " "
        end
        dfprint(io, crayon"cyan", "\t$fs", crayon"yellow"," => ", crayon"magenta", "$v\n")
    end
    dfprint(io, crayon"reset")
end

function show(io::IO, flag_info::QEFlagInfo{T}) where T
    df_show_type(io, flag_info)
    dfprintln(io, "name = $(flag_info.name)")
    dfprintln(io, crayon"cyan", "description:", crayon"reset")
    dfprint(io, "\t"*replace(flag_info.description, "\n" => "\n\t"))
end
