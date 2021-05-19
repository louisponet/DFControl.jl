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
    k_points of length $(length(band.k_points_cart)):
    cart:    $(band.k_points_cart[1]) -> $(band.k_points_cart[end])
    """
    !isempty(band.k_points_cryst) && (string *= "cryst:   $(band.k_points_cryst[1]) -> $(band.k_points_cryst[end])\n")
    string *= """eigvals: $(band.eigvals[1]) -> $(band.eigvals[end])
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
    if !isempty(is)
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
    end
    dfprint(io, reset)
end


function show(io::IO, in::DFInput)
    df_show_type(io, in)
    s = """name  = $(in.name)
    dir   = $(in.dir)
    execs = $(join([e.exec for e in in.execs],", "))
    run   = $(in.run)
    data  = $([e.name for e in data(in)])
    $(crayon"cyan")flags$(crayon"reset"):"""
    dfprintln(io, s)
    if !isempty(flags(in))
        fl = maximum(length.(string.(keys(flags(in)))))
        for (f, v) in flags(in)
            fs = string(f)
            l = length(fs)
            for i=1:fl - l
                fs *= " "
            end
            dfprint(io, crayon"cyan", "\t$fs", crayon"yellow"," => ", crayon"magenta", "$v\n")
        end
    end
    dfprint(io, crayon"reset")
end

function show(io::IO, flag_info::QEFlagInfo{T}) where T
    df_show_type(io, flag_info)
    dfprintln(io, crayon"yellow", "name : $(flag_info.name)")
    dfprintln(io, crayon"cyan", "description:", crayon"reset")
    dfprint(io, "\t"*replace(flag_info.description, "\n" => "\n\t"))
end

function show(io::IO, flag_info::ElkFlagInfo{T}) where T
    df_show_type(io, flag_info)
    dfprintln(io, crayon"yellow", "name : $(flag_info.name)")
    dfprintln(io, crayon"green", "default:", crayon"reset")
    if flag_info.default != nothing
	    dfprintln(io, "\t" * "$(flag_info.default)")
    else
	    dfprintln(io, "\t" * "nothing")
    end
    dfprintln(io, crayon"cyan", "description:", crayon"reset")
    dfprintln(io, "\t"*replace(flag_info.description, "\n" => "\n\t"))
end

function show(io::IO, info::ElkControlBlockInfo)
    dfprintln(io, crayon"yellow", "name = $(info.name)")
    dfprintln(io, crayon"cyan", "flags:", crayon"reset")
    for flag in info.flags
	    dfprintln(io, "\t$(flag.name)")
    end
    dfprintln(io, crayon"cyan", "description:", crayon"reset")
    dfprintln(io, "\t"*replace(info.description, "\n" => "\n\t"))
end

function show(io::IO, el::Element)
	for f in fieldnames(typeof(el))
		dfprintln(io, crayon"red", "$f: ", crayon"reset", "$(getfield(el, f))")
	end
end

function show(io::IO, str::AbstractStructure)
    dfprintln(io, crayon"cyan", "Structure")
    dfprintln(io, crayon"red","    cell parameters:")
    dfprint(io, crayon"reset", "\t a = $((str.cell[:,1]...,))\n\t b = $((str.cell[:,2]...,))\n\t c = $((str.cell[:,3]...,))\n")
    dfprintln(io, crayon"red","    nat:", crayon"reset", " $(length(str.atoms))")
    dfprintln(io, crayon"red","    ntyp:", crayon"reset", " $(length(unique(str.atoms)))")
    for a in atoms(str)
        show(io, a)
    end
end

function show(io::IO, at::AbstractAtom{T, LT}) where {T,LT<:Length{T}}
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
	dfprintln(io)
end
