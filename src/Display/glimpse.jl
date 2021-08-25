using ..Glimpse
const Gl = Glimpse

function Gl.Diorama(str::Structure)
    dio = Diorama(; name = Symbol(str.name))
    add_structure!(dio, str)
    dio[Gl.PointLight].data[1] = Gl.PointLight(; diffuse = 1.5f0, ambient = 1.5f0)
    Gl.center_cameras(dio)
    return dio
end
function add_structure!(dio::Diorama, str::Structure; polyhedra::Vector{Symbol} = Symbol[])
    Entity(dio,
           Gl.assemble_wire_axis_box(; position = zero(Point3f0), x = ustrip(a(str)),
                                     y = ustrip(b(str)), z = ustrip(c(str)),
                                     color = Gl.BLACK)...)
    for at in atoms(str)
        origin = Point3f0(ustrip(at.position_cart)...)
        color = RGB{Float32}(at.element.color...)
        Entity(dio,
               Gl.assemble_sphere(origin; color = color,
                                  radius = Float32(0.2 + at.element.atomic_weight / 400.0f0))...,
               Gl.Text(; str = String(at.name), font_size = 2))
        if sum(magnetization(at)) != 0
            Entity(dio,
                   Gl.assemble_arrow(origin, origin + Point3f0(magnetization(at)...);
                                     color = color, radius_ratio = 1.5f0,
                                     thickness = 0.05f0)...)
        end
    end
    for el in element.(polyhedra)
        for at in atoms(str, el)
            distances = sort(distance.((at,), atoms(str)))
            count = length(filter(i -> distances[i+1] - distances[i] < 1e-1angstrom),
                           1:length(distances)-1)
            poly = polyhedron(str, count)
        end
    end
    return dio
end
