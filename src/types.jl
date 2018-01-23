"Point in 3D space in cartesian coordinates with specified float type"
struct Point3D{T<:AbstractFloat}
    x::T
    y::T
    z::T
end
Point3D()                                   = Point3D(0.0)
Point3D(x::T) where T<:AbstractFloat        = Point3D{T}(x, x, x)
Point3D(::Type{T},x) where T<:AbstractFloat = Point3D{T}(x, x, x)
Point3D(x::Array{<:AbstractFloat,1})        = Point3D(x[1], x[2], x[3])
Point3D{T}() where T<:AbstractFloat         = Point3D{T}(0)

import Base: +, -, *, /, convert, promote_rule, show, zero, norm
+(x::Point3D, y::Point3D) = Point3D(x.x + y.x, x.y + y.y, x.z + y.z)
-(x::Point3D, y::Point3D) = Point3D(x.x - y.x, x.y - y.y, x.z - y.z)
*(x::Point3D, y::Point3D) = Point3D(x.x * y.x, x.y * y.y, x.z * y.z)
*(x::Point3D, y::Number)  = Point3D(x.x * y, x.y * y, x.z * y)
*(y::Number, x::Point3D)  = Point3D(x.x * y, x.y * y, x.z * y)
*(x::Matrix, y::Point3D)  = Point3D(x * Array(y))
/(a::Point3D, b::Point3D) = Point3D(a.x / b.x, a.y / b.y, a.z / b.z)
/(a::Point3D, b::Number)  = Point3D(a.x / b, a.y / b, a.z / b)
@inline norm(a::Point3D)  = sqrt(a.x^2 + a.y^2 + a.z^2)
zero(::Type{Point3D{T}}) where T = Point3D(zero(T), zero(T), zero(T))
zero(::Type{Point3D})     = Point3D(zero(Float64), zero(Float64), zero(Float64))

convert(::Type{Point3D}, x::Point3D) = x
convert(::Type{Array}, x::Point3D)   = [x.x, x.y, x.z]
convert(::Type{Point3D}, x::T) where T<:AbstractFloat             = Point3D{T}(x, x, x)
convert(::Type{Point3D{T}}, x::Real) where T<:AbstractFloat       = Point3D{T}(x, x, x)
convert(::Type{Point3D{T}}, x::Point3D) where T<:AbstractFloat    = Point3D{T}(x.x, x.y, x.z)
convert(::Type{Point3D{T}}, x::Array{T,1}) where T<:AbstractFloat = Point3D{T}(x[1], x[2], x[3])
promote_rule(::Type{Point3D{S}}, ::Type{T}) where {S<:AbstractFloat,T<:Real}                   = Point3D{promote_type(S,T)}
promote_rule(::Type{Point3D{S}}, ::Type{Point3D{T}}) where {S<:AbstractFloat,T<:AbstractFloat} = Point3D{promote_type(S,T)}

show(io::IO, x::Point3D)=print(io, "x = $(x.x), y = $(x.y), z = $(x.z)")
Base.write(f::IO, x::Point3D) = write(f, "$(x.x) $(x.y) $(x.z)")


include("atom.jl")
include("structure.jl")
abstract type Band end

"""
Energy band from DFT calculation.
"""
mutable struct DFBand{T<:AbstractFloat} <: Band
    k_points_cart::Array{Array{T,1},1}
    k_points_cryst::Array{Array{T,1},1}
    eigvals::Array{T,1}
    extra::Dict{Symbol,Any}
end
DFBand(k_points_cart::Array{Array{T,1},1}, k_points_cryst::Array{Array{T,1},1}, eigvals::Array{T,1}) where T <: AbstractFloat = DFBand{T}(k_points_cart, k_points_cryst, eigvals, Dict{Symbol,Any}())


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

#these are all the control blocks, they hold the flags that guide the calculation
abstract type Block end
abstract type ControlBlock <: Block end

mutable struct QEControlBlock <: ControlBlock
    name::Symbol
    flags::Dict{Symbol,Any}
end

function Base.display(block::ControlBlock)
    dfprintln("Block name: $(block.name)\n  flags:")
    for (flag, value) in block.flags
        dfprintln("    $flag => $value")
    end
    dfprintln("")
end

#these are all the data blocks, they hold the specific data for the calculation
abstract type DataBlock <: Block end

mutable struct QEDataBlock <: DataBlock
    name::Symbol
    option::Symbol
    data::Any
end

mutable struct WannierDataBlock <: DataBlock
    name::Symbol
    option::Symbol
    data::Any
end

mutable struct AbinitDataBlock <: DataBlock
    name::Symbol
    option::Symbol
    data::Any
end

function block(blocks::Array{<:Block,1}, name::Symbol)
    found_blocks = filter(x-> x.name == name, blocks)
    if isempty(found_blocks)
        return nothing
    else
        return found_blocks[1]
    end
end

function Base.display(block::DataBlock)
    s = """Block name: $(block.name)
    Block option: $(block.option)
    Block data:
    """
    dfprintln(s)
    dfprintln(string(block.data) * "\n\n")
end


function Base.display(blocks::Array{<:Block})
    map(display, blocks)
end
#here all the different input structures for the different calculations go
"""
Represents an input for DFT calculation.

Fieldnames: backend::Symbol -> the DFT package that reads this input.
control_blocks::Dict{Symbol,Dict{Symbol,Any}} -> maps different control blocks to their Dict of flags and values.
pseudos::Dict{Symbol,String} -> maps atom symbol to pseudo input file.
cell_param::Dict{Symbol,Any} -> maps the option of cell_parameters to the cell parameters.
atoms::Dict{Symbol,Any} -> maps atom symbol to position.
k_points::Dict{Symbol,Any} -> maps option of k_points to k_points.
"""
abstract type DFInput end

mutable struct QEInput <: DFInput
    filename       ::String
    structure      ::Union{Structure, Void} 
    control_blocks ::Array{QEControlBlock,1}
    data_blocks    ::Array{QEDataBlock,1}
    run_command    ::String  #everything before < in the job file
    exec           ::String
    run            ::Bool
end

function QEInput(template::QEInput, filename, newflags...; run_command=template.run_command, run=true, new_data...)
    newflags = Dict(newflags...) # this should handle both OrderedDicts and pairs of flags

    input             = deepcopy(template)
    input.filename    = filename
    input.run_command = run_command
    set_flags!(input, newflags...)

    for (block_name, block_info) in new_data
        if get_block(input, block_name) != nothing
            block = get_block(input, block_name)
            if length(block_info) == 1
                block.option = :none 
                block.data   = block_info
            elseif length(block_info) == 2
                block.option = block_info[1]
                block.data   = block_info[2]
            end
        else
            if length(block_info) == 1
                add_block!(input, QEDataBlock(block_name, :none, block_info))
            elseif length(block_info) == 2
                add_block!(input, QEDataBlock(block_name, block_info[1], block_info[2]))
            end
        end
    end
    return input
end

mutable struct WannierInput <: DFInput
    filename    ::String
    structure   ::Structure
    flags       ::Dict{Symbol,Any}
    data_blocks ::Array{WannierDataBlock,1}
    run_command ::String
    run         ::Bool
    preprocess  ::Bool
end

mutable struct AbinitInput <: DFInput
    filename    ::String
    structure   ::Union{Structure, Void}
    flags       ::Dict{Symbol,Any}
    data_blocks ::Array{AbinitDataBlock,1}
    run_command ::String
    run         ::Bool
end 

function Base.display(input::DFInput)
    print_info(input)
end

"""
Represents a full DFT job with multiple input files and calculations.

Fieldnames: name::String
calculations::Dict{String,DFInput} -> calculation type to DFInput
flow::Array{Tuple{String,String},1} -> flow chart of calculations. The tuple is (calculation type, input file).
local_dir::String -> directory on local machine.
server::String -> server in full host@server t.
server_dir::String -> directory on server.
"""
mutable struct DFJob
    id::Int
    name::String
    structure::Structure
    calculations::Array{DFInput,1}
    local_dir::String
    server::String
    server_dir::String
    header::Array{String,1}
    function DFJob(name, structure, calculations, local_dir, server,server_dir, header = get_default_job_header())
        if local_dir != ""
            local_dir = form_directory(local_dir)
        end

        if server_dir != ""
            server_dir = form_directory(server_dir)
        end

        test = filter(x -> x.name == name,UNDO_JOBS)
        if length(test) == 1
            job = new(test[1].id, name, structure, calculations, local_dir, server, server_dir, header)
            UNDO_JOBS[test[1].id] = deepcopy(job)
        elseif length(test) == 0
            job = new(length(UNDO_JOBS) + 1, name, structure, calculations, local_dir, server, server_dir, header)
            push!(UNDO_JOBS, deepcopy(job))
        end
        job
    end
end

"""
    DFJob(job_name, local_dir, args...; server=get_default_server(),server_dir="")

Creates a new DFJob, possibly passing in calculations in args... 
When inputs (args) are passed in, the structure of the job will be set to the first found in the inputs. 
"""
function DFJob(job_name, local_dir, args...; server=get_default_server(), server_dir="")
    local_dir = form_directory(local_dir)
    inputs    = DFInput[]
    structure = nothing
    for arg in args
        push!(inputs,arg)
        if arg.structure != nothing && structure != nothing
            structure = arg.structure
        end
    end
    return DFJob(job_name, structure, inputs, local_dir, server, server_dir)
end

#TODO implement abinit
# function DFJob(job_name, local_dir, calculations::Array{Pair{Union{Symbol, String}, Dict},1}, atoms, cell_parameters=eye(3);
function DFJob(job_name, local_dir, calculations::Array, atoms, cell_parameters=eye(3);
                    server=get_default_server(), 
                    server_dir="", 
                    package=:qe,
                    bin_dir="~/bin/",
                    run_command="mpirun -np 24",
                    common_flags=Dict{Symbol,Any}(),
                    pseudo_set=:default,
                    pseudo_specifier="",
                    header=get_default_job_header())

    @assert package==:qe "Only implemented for Quantum Espresso!"
    local_dir = form_directory(local_dir)
    job_atoms = convert_2atoms(atoms,pseudo_set=pseudo_set, pseudo_specifier=pseudo_specifier)
    job_calcs = DFInput[]
    structure = Structure()
    if typeof(common_flags) != Dict
        common_flags = Dict(common_flags) 
    end

    req_flags = Dict(:prefix  => "'$job_name'",
                     :outdir => "'$server_dir'",
                     :ecutwfc => 25.)
    merge!(req_flags, common_flags)    
    for (calc, data) in calculations
        calc_ = typeof(calc) == String ? Symbol(calc) : calc
        if in(calc_, [Symbol("vc-relax"), :relax, :scf])
            k_points = pop!(data, :k_points, [1, 1, 1, 0, 0, 0])
            k_option = :automatic
        elseif calc_ == :nscf
            k_points = pop!(data, :k_points, (1, 1, 1))
            k_grid   = gen_k_grid(k_points..., :nscf)
            k_option = :crystal
        elseif calc_ == :bands
            k_points = pop!(data, :k_points, [0., 0., 0., 1.])
            num_k = 0.0
            for point in k_points
                num_k += point[4]
            end
            if num_k > 100.
                push!(data[:flags], :verbosity => "'high'")
            end
            k_option = :crystal_b
        end
        flags  = pop!(data, :flags, Dict{Symbol, Any}())
        push!(flags, :calculation => "'$(string(calc_))'")
        input_ = QEInput(string(calc_) * ".in",
                         structure,
                         QEControlBlock[], 
                         [QEDataBlock(:k_points, k_option, k_points)],
                         run_command,
                         bin_dir * "pw.x",
                         true)
        set_flags!(input_, req_flags..., print=false)
        set_flags!(input_, flags..., print=false)
        change_atoms!(input_, job_atoms, pseudo_set = pseudo_set, pseudo_specifier = pseudo_specifier, print=false)
        push!(job_calcs, input_)
    end
    return DFJob(job_name, job_calcs, local_dir, server, server_dir, header)
end

function Base.display(job::DFJob)
    try
        print_info(job)
    end
end
