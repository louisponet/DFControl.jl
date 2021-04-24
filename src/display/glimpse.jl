@reexport using Glimpse

function Glimpse.Diorama(str::AbstractStructure)
    dio = Diorama(name=Symbol(str.name))
    Entity(dio, assemble_wire_axis_box(zero(Point3f0), a(str), b(str), c(str)...))
    for a in atoms(str)
        
    Entity(dio, 
end
