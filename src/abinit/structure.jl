using DFControl: Structure, Atom, Mat3, Point3, element
# struct AbiStructure{T <: AbstractFloat} <: AbstractStructure{T}
#     abi  ::PyObject
#     core ::Structure{T}
# end
#
# AbiStructure(abi::PyObject) = AbiStructure(abi, Structure(abi))
#
# structure(abi::AbiStructure) = structure(abi.core)
function Structure(abi_structure::PyObject; name="NoName")
    cell = Mat3(abi_structure[:lattice][:matrix])
    syms_positions = abi_structure[:get_symbol2coords]()
    atoms = Atom{Float64}[]
    for (s, p) in syms_positions
        atsym = Symbol(s)
        el    = element(atsym)
        position = [cell' * Point3(p[i,:]) for i = 1:size(p, 2)]
        for pos in position
            push!(atoms, Atom(atsym, el, pos))
        end
    end
    return Structure(name, cell, atoms)
end
