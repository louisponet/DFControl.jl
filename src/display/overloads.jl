#printing that is not needed in Atom


function Base.show(io::IO, block::InputData)
    s = """Block name: $(block.name)
    Block option: $(block.option)
    Block data:
    """
    dfprintln(io, s)
    dfprintln(io, string(block.data) * "\n\n")
end

Base.show(io::IO, data::Vector{InputData}) = map(x-> show(io, x), data)

function Base.show(io::IO, band::DFBand{T}) where T <: AbstractFloat
    string = """DFBand{$T}:
    k_points of length $(length(band.k_points_cryst)):
    cart:    $(band.k_points_cart[1]) -> $(band.k_points_cart[end])
    cryst:   $(band.k_points_cryst[1]) -> $(band.k_points_cryst[end])
    eigvals: $(band.eigvals[1]) -> $(band.eigvals[end])
    extra:   $(band.extra)
    """
    dfprintln(io, string)
end

Base.show(io::IO, bands::Vector{<:DFBand}) = map(x->show(io,x),bands)

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
    for i in inputs(job)
        dfprintln(io, "\t$(i.filename) => runs=$(i.run)")
    end
end

Base.show(io::IO, at::AbstractAtom) = dfprintln(io, "$(id(at)): $(position(at)[1]) $(position(at)[2]) $(position(at)[3])")
