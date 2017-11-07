"File with all the user defaults inside it"
const default_file = joinpath(@__DIR__,"../user_defaults/user_defaults.jl")

function init_defaults(filename::String)
  raw_input =""
  names_to_export = Symbol[] 
  open(filename,"r") do f
    while !eof(f)
      line = readline(f)
      if line == "" || line[1] == '#'
        continue
      end
      lhs = parse(line).args[1]
      if typeof(lhs) == Symbol
        push!(names_to_export,lhs)
      end
      raw_input *= line*"; "
    end
  end
  for name in names_to_export
    eval(:(export $name))
  end
  eval(parse(raw_input))
end

function load_defaults(filename::String=default_file)
  raw_input =""
  names_to_export = Symbol[] 
  open(filename,"r") do f
    while !eof(f)
      raw_input *= readline(f)*"; "
    end
  end
  raw_input *= "nothing ;"
  eval(parse(raw_input))
end

"Macro which allows you to define any default variable that will get loaded every time you use this package."
macro set_default(expr)
  expr2file(default_file,expr)
  load_defaults()
end

function define_def(default, expr1, expr2)
  if !isdefined(default)
    expr2file(default_file,expr1)
    init_defaults(default_file)
  else
    expr2file(default_file,expr2)
    load_defaults(default_file)
  end
end

"""
    add_default_pseudo_dir(pseudo_symbol::Symbol, dir::String)

Adds an entry inside the `default_pseudo_dirs` dictionary with flag `pseudo_symbol`.
"""
function add_default_pseudo_dir(pseudo_symbol::Symbol, dir::String)
  expr_ndef = :(default_pseudo_dirs = Dict{Symbol,String}($(Expr(:quote,pseudo_symbol)) => $dir)) 
  expr_def  = :(default_pseudo_dirs[$(QuoteNode(pseudo_symbol))] = $dir)
  define_def(:default_pseudo_dirs,expr_ndef,expr_def)
end

"""
    remove_default_pseudo_dir(pseudo_symbol::Symbol)

Removes entry with flag `pseudo_symbol` from the `default_pseudo_dirs` dictionary. 
"""
function remove_default_pseudo_dir(pseudo_symbol::Symbol)
  if isdefined(:default_pseudo_dirs) && haskey(default_pseudo_dirs,pseudo_symbol)
    pop!(default_pseudo_dirs,pseudo_symbol)
    if isempty(default_pseudo_dirs)
      rm_expr_lhs(default_file,:default_pseudo_dirs)
    else
      rm_expr_lhs(default_file,:(default_pseudo_dirs[$(QuoteNode(pseudo_symbol))]))
    end
  end
  load_defaults(default_file)
end

"""
    set_default_server(server::String)

Sets the default server.
"""
function set_default_server(server::String)
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
function configure_default_pseudos(server = get_default_server(), pseudo_dirs = get_default_pseudo_dirs())
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
 
  if !isdefined(:default_pseudos)
    expr2file(default_file,:(default_pseudos = Dict{Symbol,Dict{Symbol,Array{String,1}}}()))
    init_defaults(default_file)
  end
  
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

"""
     set_default_job_header(lines)
    
Sets the header that will get added to each job.tt file. 
"""
function set_default_job_header(lines)
  expr = :(default_job_header = $lines)
  expr2file(default_file,expr)
  if !isdefined(:default_job_header)
    init_defaults(default_file)
  else
    load_defaults(default_file)
  end
end

"""
     set_default_input(input::DFInput, calculation::Symbol)

Adds the input to the default inputs, writes it to a file in user_defaults folder to be read every time on load.
"""
function set_default_input(input::DFInput, calculation::Symbol)
  if !isdefined(:default_inputs)
    expr = :(default_inputs = Dict{Symbol,DFInput}())
    expr2file(default_file,expr)
    init_defaults(default_file)
  end
  filename = dirname(default_file)*"/"*String(calculation)
  if typeof(input) == WannierInput
    write_input(input,filename * ".win")
    expr2file(default_file,:(default_inputs[$(QuoteNode(calculation))] = read_wannier_input($filename * ".win",run_command=$(input.run_command))))
  elseif typeof(input) == QEInput
    write_input(input,filename * ".in")
    expr2file(default_file,:(default_inputs[$(QuoteNode(calculation))] = read_qe_input($filename * ".in",run_command = $(input.run_command))))
  end
  load_defaults(default_file)
end

"""
    remove_default_input(input::Symbol)
    
Remove the default input specified by the Symbol. Also removes the stored input file.
"""
function remove_default_input(input::Symbol)
  if haskey(default_inputs,input)
    input = pop!(default_inputs,input)
    if isempty(default_inputs)
      rm_expr_lhs(default_file,:default_inputs)
      default_inputs = nothing
    else
      rm_expr_lhs(default_file,:(default_inputs[$(QuoteNode(input))]))
    end
    rm(joinpath(@__DIR__,"../user_defaults/$(input.filename)"))
  else
    error("Default_calculations does not have an input with symbol $symbol.\n  Possible symbols are: $(keys(default_inputs))")
  end
end
  

