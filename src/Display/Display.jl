module Display
using Requires, Crayons, Dates, ANSIColoredPrinters
using ..Structures
using ..Calculations
using ..Jobs
using ..Client
using ..DFControl: Band, TimingData
using ..Utils
using RemoteHPC: Server, isalive

include("printing.jl")
include("pluto.jl")
include("base.jl")

const dfprintln = println
const dfprint = print

function __init__()
    @require Glimpse = "f6e19d58-12a4-5927-8606-ac30a9ce9b69" include("glimpse.jl")
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plotting.jl")
end
end
