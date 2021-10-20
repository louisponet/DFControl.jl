using .Glimpse
const Gl = Glimpse

using ..Structures
using ..DFControl

function Gl.Diorama(str::Structure)
    dio = Diorama(; name = Symbol("DFControl Structure"))
    add_structure!(dio, str)
    dio[Gl.PointLight].data[1] = Gl.PointLight(; diffuse = 1.5f0, ambient = 1.5f0)
    Gl.center_cameras(dio)
    return dio
end
function add_structure!(dio::Diorama, str::Structure; polyhedra::Vector{Symbol} = Symbol[])
    Entity(dio,
           Gl.assemble_wire_axis_box(; position = zero(Point3f0), x = Structures.ustrip(Structures.a(str)),
                                     y = Structures.ustrip(Structures.b(str)), z = Structures.ustrip(Structures.c(str)),
                                     color = Gl.BLACK)...)
    for at in str.atoms
        origin = Point3f0(Structures.ustrip(at.position_cart)...)
        color = RGB{Float32}(at.element.color...)
        Entity(dio,
               Gl.assemble_sphere(origin; color = color,
                                  radius = Float32(0.2 + at.element.atomic_weight / 400.0f0))...,
               Gl.Text(; str = String(at.name), font_size = 2))
        if sum(at.magnetization) != 0
            Entity(dio,
                   Gl.assemble_arrow(origin, origin + Point3f0(at.magnetization...);
                                     color = color, radius_ratio = 1.5f0,
                                     thickness = 0.05f0)...)
        end
    end
    for el in map(x->x.element, polyhedra)
        for at in str[el]
            distances = sort(distance.((at,), str.atoms))
            count = length(filter(i -> distances[i+1] - distances[i] < 1e-1angstrom),
                           1:length(distances)-1)
            poly = polyhedron(str, count)
        end
    end
    return dio
end
