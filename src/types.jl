"Point in 3D space in cartesian coordinates with specified float type"
struct Point3D{T<:AbstractFloat}
  x::T
  y::T
  z::T
end
(Point3D(x::T) where T<:AbstractFloat) = Point3D(T(x),T(x),T(x))
Point3D(::Type{T},x) where T<:AbstractFloat = Point3D(T(x),T(x),T(x))
Point3D(x::Array{<:AbstractFloat,1}) = Point3D(x[1],x[2],x[3])

abstract type Band{T<:AbstractFloat} end

"""
Energy band from DFT calculation.
"""
mutable struct DFBand{T<:AbstractFloat} <: Band{T}
  k_points_cart::Array{Array{T,1},1}
  k_points_cryst::Array{Array{T,1},1}
  eigvals::Array{T,1}
end

"""
Represents an input for DFT calculation.

Fieldnames: backend::Symbol -> the DFT package that reads this input.
            control_blocks::Dict{Symbol,Dict{Symbol,Any}} -> maps different control blocks to their dict of flags and values.
            pseudos::Dict{Symbol,String} -> maps atom symbol to pseudo input file.
            cell_param::Dict{Symbol,Any} -> maps the option of cell_parameters to the cell parameters.
            atoms::Dict{Symbol,Any} -> maps atom symbol to position.
            k_points::Dict{Symbol,Any} -> maps option of k_points to k_points.
"""
#these are all the control blocks, they hold the flags that guide the calculation
abstract type ControlBlock end

mutable struct QEControlBlock<:ControlBlock
  name::Symbol
  flags::Dict{Symbol,Any}
end

function Base.display(block::ControlBlock)
  println("Block name: $(block.name)")
  println("Block flags:")
  display(block.flags)
end

#these are all the data blocks, they hold the specific data for the calculation
abstract type DataBlock end

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

function Base.display(block::DataBlock)
  println("Block name: $(block.name)")
  println("Block option: $(block.option)")
  println("Block data:")
  display(block.data)
end

#here all the different input structures for the different calculations go
abstract type DFInput end

mutable struct QEInput<:DFInput
  filename::String
  control_blocks::Array{QEControlBlock,1}
  data_blocks::Array{QEDataBlock,1}
  run_command::String  #everything before < in the job file
  run::Bool
end

mutable struct WannierInput<:DFInput
  filename::String
  flags::Dict{Symbol,Any}
  data_blocks::Array{WannierDataBlock,1}
  run_command::String
  run::Bool
end

# mutable struct DFInput
#   backend::Symbol
#   control_blocks::Dict{Symbol,Any}
#   pseudos::Dict{Symbol,String}
#   cell_param::Dict{Symbol,Any}
#   atoms::Dict{Symbol,Union{Symbol,Array{<:Point3D,1}}}
#   k_points::Dict{Symbol,Any}
# end

#having both the flow and array tuple might be overkill
"""
Represents a full DFT job with multiple input files and calculations.

Fieldnames: job_name::String
            calculations::Dict{String,DFInput} -> calculation type to DFInput
            flow::Array{Tuple{String,String},1} -> flow chart of calculations. The tuple is (calculation type, input file).
            local_dir::String -> directory on local machine.
            server::String -> server in full host@server t.
            server_dir::String -> directory on server.
"""
mutable struct DFJob
  job_name::String
  calculations::Array{DFInput,1}
  # flow::Array{Tuple{String,String},1}
  local_dir::String
  server::String
  server_dir::String
  function DFJob(job_name,calculations,local_dir,server,server_dir)
  # function DFJob(job_name,calculations,flow,local_dir,server,server_dir)
    if local_dir != ""
      local_dir = form_directory(local_dir)
    end
    if server_dir != ""
      server_dir = form_directory(server_dir)
    end
    new(job_name,calculations,local_dir,server,server_dir)
    # new(job_name,calculations,flow,local_dir,server,server_dir)
  end
end

"""
Represents an element.
"""
struct Element
  Z::Int64
  Name::String
  atomic_weight::Float64
end

"""
Reads all the elements from the file.
"""
const ELEMENTS = Dict()
open(joinpath(@__DIR__,"../assets/elements.txt"),"r") do f
  while !eof(f)
    line = split(readline(f))
    ELEMENTS[Symbol(line[4])] = Element(parse(Int64,line[1]),line[9],parse(Float64,line[10]))
  end
end
