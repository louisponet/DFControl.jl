module Display
using Requires, Crayons, Dates, ANSIColoredPrinters
using ..Structures
using ..Calculations
using ..Jobs
using ..Client
using ..DFControl: Band, TimingData

include("printing.jl")

const dfprintln = println
const dfprint = print

function __init__()
    @require Glimpse = "f6e19d58-12a4-5927-8606-ac30a9ce9b69" include("glimpse.jl")
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plotting.jl")
    if isdefined(Base, :active_repl)
        Display.include(joinpath(@__DIR__, "base.jl"))
    elseif isdefined(Main, :PlutoRunner)
        Display.include(joinpath(@__DIR__, "pluto.jl"))
        
    else
        @require IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a" include("base.jl")
    end
end
end
