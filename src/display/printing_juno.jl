
import Media
import Juno: @render, Tree, SubTree, Inline, Model,  render
import Atom: span, c, UNDEF
getfield_(x, f) = isdefined(x, f) ? getfield(x, f) : UNDEF

@render i::Inline x::InputData begin
    Tree(x.name, [SubTree(Text("$f → "), getfield(x, f)) for f in fieldnames(typeof(x))[2:end]])
end

@render i::Inline x::DFCalculation begin
    runstring = x.run ? span(".syntax--constant.syntax--other.syntax--symbol", name(x)) : span(".syntax--comment", name(x))
    fields = filter(x->x != run, [fieldnames(typeof(x))[2:end]...])
    Tree(runstring , [SubTree(Text("$f → "), getfield(x, f)) for f in fields])
end

@render i::Inline x::AbstractAtom begin
    pos = position_cart(x)
    Tree(Text("$(name(x))  ang: $(round(pos[1], digits=5)) $(round(pos[2], digits=5)) $(round(pos[3], digits=5))"), [SubTree(Text("$f → "), getfield_(x, f)) for f in fieldnames(typeof(x))[2:end]])
end

@render i::Inline x::DFJob begin
    Tree(x.name, [SubTree(Text("$f → "), getfield(x, f)) for f in fieldnames(typeof(x))[3:end]])
end

@render i::Inline x::AbstractStructure begin
    Tree(name(x) , [SubTree(Text("$f → "), getfield(x, f)) for f in fieldnames(typeof(x))])
end
