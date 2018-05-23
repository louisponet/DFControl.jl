#printing that is not needed in Atom

function Base.display(block::ControlBlock)
    dfprintln("Block name: $(block.name)\n  flags:")
    for (flag, value) in block.flags
        dfprintln("    $flag => $value")
    end
    dfprintln("")
end

function Base.show(io::IO, block::ControlBlock)
    println("Block name: $(block.name)\n  flags:")
    for (flag, value) in block.flags
        println("    $flag => $value")
    end
    println("")
end

function Base.display(block::DataBlock)
    s = """Block name: $(block.name)
    Block option: $(block.option)
    Block data:
    """
    dfprintln(s)
    dfprintln(string(block.data) * "\n\n")
end

function Base.display(data::Array{<:Block})
    map(display, data)
end

function Base.display(input::DFInput)
    print_info(input)
end

function Base.display(band::DFBand{T}) where T <: AbstractFloat
    string = """DFBand{$T}:
    k_points of length $(length(band.k_points_cryst)):
    cart:    $(band.k_points_cart[1]) -> $(band.k_points_cart[end])
    cryst:   $(band.k_points_cryst[1]) -> $(band.k_points_cryst[end])
    eigvals: $(band.eigvals[1]) -> $(band.eigvals[end])
    extra:   $(band.extra)
    """
    dfprintln(string)
end

function Base.display(bands::Array{<:DFBand})
    map(display,bands)
end

function Base.show(job::DFJob)
    try
        print_info(job)
    end
end
