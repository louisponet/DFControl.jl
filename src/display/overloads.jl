#printing that is not needed in Atom


function Base.show(block::InputData)
    s = """Block name: $(block.name)
    Block option: $(block.option)
    Block data:
    """
    dfprintln(s)
    dfprintln(string(block.data) * "\n\n")
end
function Base.display(block::InputData)
    s = """Block name: $(block.name)
    Block option: $(block.option)
    Block data:
    """
    dfprintln(s)
    dfprintln(string(block.data) * "\n\n")
end

function Base.display(data::Vector{InputData})
    map(display, data)
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

function Base.display(bands::Vector{<:DFBand})
    map(display,bands)
end

function Base.show(io::IO, job::DFJob)
    s = """--------------------
    DFJob:      $(job.name)
    Local_dir:  $(job.local_dir)
    Server:     $(job.server)
    Server_dir: $(job.server_dir)
    $(length(job.inputs)) calculations
    --------------------
    """
    dfprintln(io, s)
end
