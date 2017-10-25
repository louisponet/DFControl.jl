
const default_file = joinpath(@__DIR__,"../user_defaults/user_defaults.jl")

function define_def(default,expr1,expr2)
  if !isdefined(default)
    expr2file(default_file,expr1)
  else
    expr2file(default_file,expr2)
  end
  load_defaults(default_file)
end

function add_default_pseudo_dir!(pseudo_symbol::Symbol, dir::String)
  expr_ndef = :(default_pseudo_dirs = Dict{Symbol,String}($(Expr(:quote,pseudo_symbol)) => $dir)) 
  expr_def  = :(default_pseudo_dirs[$(QuoteNode(pseudo_symbol))] = $dir)
  define_def(:default_pseudo_dirs,expr_ndef,expr_def)
end

function remove_default_pseudo_dir!(pseudo_symbol::Symbol)
  if isdefined(:default_pseudo_dirs) && haskey(default_pseudo_dirs,pseudo_symbol)
    pop!(default_pseudo_dirs,pseudo_symbol)
  end
end

function set_default_server!(server::String)
  expr_ndef = :(default_server = $server)
  expr_def  = expr_ndef
  define_def(:default_server,expr_ndef,expr_def)
end

function g_default_server()
  if isdefined(:default_server)
    return default_server
  else
    return ""
  end
end

function g_default_pseudo_dirs()
  if isdefined(:default_pseudo_dirs)
    return default_pseudo_dirs
  else
    return nothing
  end
end

function configure_default_pseudos!(server=g_default_server(),pseudo_dirs = g_default_pseudo_dirs())
  if server == ""
    error("Either supply a valid server string or setup a default server through 'set_default_server!()'.")
  end
  if pseudo_dirs == nothing
    error("Either supply valid pseudo directories or setup a default pseudo dir through 'add_default_pseudo_dir()'.")
  end
  outputs = Dict{Symbol,String}()
  if typeof(pseudo_dirs) == String
    outputs[:default] = readstring(`ssh -t $server ls $pseudo_dirs`)
  elseif typeof(pseudo_dirs) <: Dict
    for (name, directory) in pseudo_dirs
      outputs[name] = readstring(`ssh -t $server ls $directory`)
    end
  end
  
  expr2file(default_file,:(default_pseudos = Dict{Symbol,Dict{Symbol,Array{String,1}}}()))
  
  # atoms = Dict{Symbol,Dict{Symbol,Array{String,1}}}()
  for el in keys(ELEMENTS)
    expr2file(default_file,:(default_pseudos[$(QuoteNode(el))] = Dict{Symbol,Array{String,1}}()))
  end
  for (name,pseudo_string) in outputs
    pseudos = filter(x-> x!= "",split(pseudo_string,"\n"))
    i = 1 
    while i <= length(pseudos)
      pseudo = pseudos[i]
      element = Symbol(split(pseudo,".")[1])
      t_expr = :(String[$pseudo])
      j=1
      while j+i<=length(pseudos) && Symbol(split(pseudos[i+j],".")[1]) == element
        push!(t_expr.args,pseudos[i+j])
        j+=1
      end
      i+=j
      expr2file(default_file,:(default_pseudos[$(QuoteNode(element))][$(QuoteNode(name))] = $t_expr))
    end
  end
  load_defaults(default_file)
end

function get_default_pseudo(atom::Symbol,pseudo_set_name=:default;pseudo_fuzzy=nothing)
  if isdefined(:default_pseudos)
    if pseudo_fuzzy != nothing
      return filter(x->contains(x,pseudo_fuzzy),default_pseudos[atom][pseudo_set_name])[1]
    else
      return default_pseudos[atom][pseudo_set_name][1]
    end
  end
end