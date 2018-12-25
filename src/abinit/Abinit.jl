module Abinit
    using PyCall
    const abilab_opt = PyNULL()
    using ..DFControl
    include("constants.jl")
    include("structure.jl")
    include("fileio.jl")
    function __init__()
        copy!(abilab_opt, pyimport("abipy.abilab"))
    end
end
