module Abinit
    using PyCall
    const abilab_opt = PyNULL()
    using ..DFControl

    abi_works_assert() = @assert abilab_opt != PyNULL() "`abipy` cannot be accessed, please install `abipy` before using the `Abinit` backend."

    include("constants.jl")
    include("structure.jl")
    include("fileio.jl")
    function __init__()
        try
            copy!(abilab_opt, pyimport("abipy.abilab"))
        catch
            @warn "`abipy` cannot be accessed, please install `abipy` before using the `Abinit` backend."
        end
    end
end
