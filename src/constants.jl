const assets_dir = joinpath(@__DIR__, "..", "assets")
const conversions = Dict{Symbol,Float64}(:bohr2ang => 0.52917721092)
conversions[:ang2bohr] = 1 / conversions[:bohr2ang]
# include("abinit/constants.jl")
include("qe/constants.jl")
include("wannier90/constants.jl")
include("elk/constants.jl")
