"File with all the user defaults inside it"
const default_file = joinpath(@__DIR__,"../user_defaults/user_defaults.jl")


function define_def(default, expr1, expr2)
  if !isdefined(default)
    expr2file(default_file,expr1)
  else
    expr2file(default_file,expr2)
  end
  load_defaults(default_file)
end

"""
    add_default_pseudo_dir!(pseudo_symbol::Symbol, dir::String)

Adds an entry inside the `default_pseudo_dirs` dictionary with flag `pseudo_symbol`.
"""
function add_default_pseudo_dir!(pseudo_symbol::Symbol, dir::String)
  expr_ndef = :(default_pseudo_dirs = Dict{Symbol,String}($(Expr(:quote,pseudo_symbol)) => $dir)) 
  expr_def  = :(default_pseudo_dirs[$(QuoteNode(pseudo_symbol))] = $dir)
  define_def(:default_pseudo_dirs,expr_ndef,expr_def)
end

"""
    remove_default_pseudo_dir!(pseudo_symbol::Symbol)

Removes entry with flag `pseudo_symbol` from the `default_pseudo_dirs` dictionary. 
"""
function remove_default_pseudo_dir!(pseudo_symbol::Symbol)
  if isdefined(:default_pseudo_dirs) && haskey(default_pseudo_dirs,pseudo_symbol)
    pop!(default_pseudo_dirs,pseudo_symbol)
  end
end

"""
    set_default_server!(server::String)

Sets the default server.
"""
function set_default_server!(server::String)
  expr_ndef = :(default_server = $server)
  expr_def  = expr_ndef
  define_def(:default_server,expr_ndef,expr_def)
end

"""
    get_default_server()

Returns the default server if it's defined. If it is not defined return "".
"""
function get_default_server()
  if isdefined(:default_server)
    return default_server
  else
    return ""
  end
end

"""
    get_default_pseudo_dirs()

Returns the default pseudo dirs dictionary if it's defined. If it is not defined return nothing.
"""
function get_default_pseudo_dirs()
  if isdefined(:default_pseudo_dirs)
    return default_pseudo_dirs
  else
    return nothing
  end
end

"""
    configure_default_pseudos!(server = get_default_server(), pseudo_dirs = get_default_pseudo_dirs())

Reads the specified `default_pseudo_dirs` on the `default_server` and sets up the `default_pseudo` dictionary.
"""
function configure_default_pseudos!(server = get_default_server(), pseudo_dirs = get_default_pseudo_dirs())
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

"""
    get_default_pseudo(atom::Symbol, pseudo_set_name=:default; pseudo_fuzzy=nothing)

Returns the pseudo potential string linked to the atom.
"""
function get_default_pseudo(atom::Symbol, pseudo_set_name=:default; pseudo_fuzzy=nothing)
  if isdefined(:default_pseudos)
    if pseudo_fuzzy != nothing
      return filter(x->contains(x,pseudo_fuzzy),default_pseudos[atom][pseudo_set_name])[1]
    else
      return default_pseudos[atom][pseudo_set_name][1]
    end
  end
end