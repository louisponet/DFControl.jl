@reexport using Glimpse
const Gl = Glimpse

function Gl.Diorama(str::AbstractStructure)
    dio = Diorama(name=Symbol(str.name))
    Entity(dio, Gl.assemble_wire_axis_box(position=zero(Point3f0), x=ustrip(a(str)), y=ustrip(b(str)), z=ustrip(c(str)), color=Gl.BLACK)...)
    for at in atoms(str)
        origin = Point3f0(ustrip(position_cart(at))...)
        color = RGB{Float32}(element(at).color...)
        Entity(dio, Gl.assemble_sphere(origin, color=color,radius=Float32(0.2 + element(at).atomic_weight/400f0))..., Gl.Text(str=String(name(at)), font_size=2))
        if sum(magnetization(at)) != 0
            
            Entity(dio, Gl.assemble_arrow(origin, origin + Point3f0(magnetization(at)...), color=color, radius_ratio=1.5f0, thickness = 0.05f0)...)
        end
    end
    dio[Gl.PointLight].data[1] = Gl.PointLight(diffuse=1.5f0, ambient=1.5f0)
    expose(dio)
end
