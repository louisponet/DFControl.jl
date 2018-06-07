
import Media
import Juno: @render, Tree, SubTree, Inline, Model,  render
import Atom: span, c, UNDEF
getfield_(x, f) = isdefined(x, f) ? getfield(x, f) : UNDEF

@render i::Inline x::InputData begin
    Tree(x.name, [SubTree(Text("$f → "), getfield(x, f)) for f in fieldnames(x)[2:end]])
end

@render i::Inline x::DFInput begin
    runstring = x.run ? span(".syntax--constant.syntax--other.syntax--symbol",x.filename) : span(".syntax--comment", x.filename)
    fields = filter(x->x != run, fieldnames(x)[2:end])
    Tree(runstring , [SubTree(Text("$f → "), getfield(x, f)) for f in fields])
end

@render i::Inline x::AbstractAtom begin
    pos = position(x)
    Tree(Text("$(id(x))  ang: $(round(pos[1], 5)) $(round(pos[2], 5)) $(round(pos[3], 5))"), [SubTree(Text("$f → "), getfield_(x, f)) for f in fieldnames(x)[2:end]])
end

@render i::Inline x::DFJob begin
    Tree(x.name, [SubTree(Text("$f → "), getfield(x, f)) for f in fieldnames(x)[3:end]])
end

@render i::Inline x::AbstractStructure begin
    Tree(x.name , [SubTree(Text("$f → "), getfield(x, f)) for f in fields])
end
