
import Media
import Juno: @render, Tree, SubTree, Inline, Model,  render
import Atom: span, c, UNDEF
getfield_(x, f) = isdefined(x, f) ? getfield(x, f) : UNDEF

@render i::Inline x::ControlBlock begin
    Tree(x.name, [SubTree(span(c(render(i, key), " → ")), val) for (key,val) in x.flags])
end


@render i::Inline x::DataBlock begin
    Tree(x.name, [SubTree(Text("$f → "), getfield(x, f)) for f in fieldnames(x)[2:end]])
end

@render i::Inline x::DFInput begin
    runstring = x.run ? span(".syntax--constant.syntax--other.syntax--symbol",x.filename) : span(".syntax--comment", x.filename)
    fields = filter(x->x != run, fieldnames(x)[2:end])
    Tree(runstring , [SubTree(Text("$f → "), getfield(x, f)) for f in fields])
end

@render i::Inline x::Atom begin
    Tree(Text("$(x.id)  ang: $(round(x.position.x, 5)) $(round(x.position.y, 5)) $(round(x.position.z, 5))"), [SubTree(Text("$f → "), getfield_(x, f)) for f in fieldnames(x)[2:end]])
end
