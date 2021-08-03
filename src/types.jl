abstract type Band end

"""
Energy band from DFT calculation.
"""
mutable struct DFBand <: Band
    k_points_cart  :: Vector{Vec3}
    k_points_cryst :: Vector{Vec3{Float64}}
    eigvals        :: Vector{Float64}
    extra          :: SymAnyDict
end
function DFBand(k_points_cart, k_points_cryst, eigvals)
    return DFBand(k_points_cart, k_points_cryst, eigvals, Dict{Symbol,Any}())
end

function kpoints(band::DFBand, kind = :cryst)
    return kind == :cart ? band.k_points_cart : band.k_points_cryst
end
eigvals(band::DFBand) = band.eigvals
StructTypes.StructType(::Type{<:Band}) = StructTypes.Mutable()

"""
    bandgap(bands::AbstractVector{DFBand}, fermi=0.0)

Calculates the bandgap (possibly indirect) around the fermi level.
"""
function bandgap(bands::Union{Iterators.Flatten,AbstractVector{<:Band}}, fermi = 0.0)
    max_valence = -Inf
    min_conduction = Inf
    for b in bands
        max = maximum(eigvals(b) .- fermi)
        min = minimum(eigvals(b) .- fermi)
        if max_valence <= max <= 0.0
            max_valence = max
        end
        if 0.0 <= min <= min_conduction
            min_conduction = min
        end
    end
    return min_conduction - max_valence
end
function bandgap(u_d_bands::Union{NamedTuple,Tuple}, args...)
    return bandgap(Iterators.flatten(u_d_bands), args...)
end

mutable struct TimingData
    name::String
    cpu::Dates.AbstractTime
    wall::Dates.AbstractTime
    calls::Int
    children::Vector{TimingData}
end

function Base.hash(data::T, h::UInt) where {T<:Union{<:InputData,Projection,Exec}}
    for f in fieldnames(T)
        h = hash(getfield(data, f), h)
    end
    return h
end

@enum Scheduler Slurm=1 Bash=2

@with_kw mutable struct Server
    name::String = "nothing"
    username::String = "nothing"
    domain::String = "nothing"
    port::Int = 8080
    scheduler::Scheduler = Bash
    mountpoint::String = ""
    julia_exec::String = "julia"
    default_jobdir::String = homedir()
end

function Server(s::String)
    server = Client.load_server(s) #First check if previous server exists
    if server !== nothing
        return server
    end
    #TODO load existing config 
    # Create new server 
    if occursin("@", s)
        username, domain = split(s, "@")
        name = ""
    else
        @info "Server with name $s not found."
        name = s
        username, domain = "", ""
    end
    @info "Creating new Server configuration..."
    while isempty(name)
        print("Please specify the Server's identifying name:")
        name = readline()
    end
    if load_server(name) !== nothing
        @warn "A server with $name was already configured and will be overwritten."
    end
    while isempty(username)
        print("Username:")
        username = readline()
    end
    while isempty(domain)
        print("Domain:")
        domain = readline()
    end
    print("Port (default: 8080):")
    port_str = readline()
    port = isempty(port_str) ? 8080 : parse(Int, port_str)
    
    print("Julia Exec (default: julia):")
    julia_str = readline()
    julia = isempty(julia_str) ? "julia" : julia_str

    dir = remotecall_fetch(homedir, Distributed.addprocs([("$username"*"@"*"$domain", 1)])[1])
    print("Jobs top dir (default: $dir):")
    dir = abspath(readline())
    
    scheduler_choice = request("Please select scheduler:", RadioMenu([string.(instances(Scheduler))...]))
    scheduler_choice == -1 && return 
    scheduler = Scheduler(scheduler_choice)
    mounted_choice = request("Has the server been mounted?", RadioMenu(["yes", "no"]))
    mounted_choice == -1 && return
    if mounted_choice == 1
        print("Please specify mounting point:")
        mountpoint = readline()
    else
        mountpoint = ""
    end
    server = Server(name, username, domain, port, scheduler, mountpoint, julia_str, dir)
    println("Server configured as:")
    println(server)
    save(server)
    return server
end
Server(j::Job) = Server(j.server)

StructTypes.StructType(::Type{Server}) = StructTypes.Mutable()

Base.joinpath(s::Server, p...) = joinpath(s.default_jobdir, p...)
