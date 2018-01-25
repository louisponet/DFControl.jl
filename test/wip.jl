using DFControl

using StaticArrays
using NamedTuples
#We use angstrom everywhere
mutable struct TestAtom{T <: AbstractFloat}
    id       ::Symbol
    element  ::Element
    position ::Point3D{T}
    test     ::T
    pseudo   ::AbstractString
    function TestAtom(id::Symbol, element::Element, position::Point3D{T}, args::KeyValIterable...) where T
        atom = new{T}()
        names = fieldnames(Atom)
        atom.id = id
        atom.element = element
        atom.position=position
        for name in names[4:end] 
            for (field, value) in args
                if field == name
                    setfield!(atom, field, value)
                end
            end
        end
        return atom 
    end
end

Atom(id::Symbol, element::Symbol, position::Point3D, data...)  = Atom(id, ELEMENTS[element], position, data...)
@code_warntype(test(t))