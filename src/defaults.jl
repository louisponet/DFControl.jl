"File with all the user defaults inside it"
const default_file = joinpath(@__DIR__, "../user_defaults/user_defaults.jl")

function init_defaults(filename::String)
    raw_input       = ""
    names_to_export = Symbol[]
    open(filename, "r") do f
        while !eof(f)
            line = readline(f)
            if line == "" || line[1] == '#'
                continue
            end
            lhs = parse(line).args[1]
            if typeof(lhs) == Symbol
                push!(names_to_export, lhs)
            end
            raw_input *= line * "; "
        end
    end
    for name in names_to_export
        eval(:(export $name))
    end
    eval(parse(raw_input))
end

function load_defaults(filename::String=default_file)
    raw_input       = ""
    names_to_export = Symbol[]
    open(filename, "r") do f
        while !eof(f)
            raw_input *= readline(f) * "; "
        end
    end
    raw_input *= "nothing ;"
    eval(parse(raw_input))
end

"Macro which allows you to define any default variable that will get loaded every time you use this package."
macro add_default(expr)
    expr2file(default_file, expr)
    eval(Main, expr)
    load_defaults()
end

function define_def(default, expr1, expr2)
    if !isdefined(default)
        expr2file(default_file, expr1)
        init_defaults(default_file)
    else
        expr2file(default_file, expr2)
        load_defaults(default_file)
    end
end

"""
    setdefault_pseudodir(pseudo_symbol::Symbol, dir::String)

Adds an entry inside the `default_pseudodirs` with flag `pseudo_symbol`, and adds it to the `user_defaults.jl` file.
"""
function setdefault_pseudodir(pseudo_symbol::Symbol, dir::String)
    expr_ndef = :(default_pseudo_dirs = Dict{Symbol,String}($(Expr(:quote, pseudo_symbol)) => $dir))
    expr_def  = :(default_pseudo_dirs[$(QuoteNode(pseudo_symbol))] = $dir)
    define_def(:default_pseudo_dirs, expr_ndef, expr_def)
end

"""
    removedefault_pseudodir(pseudo_symbol::Symbol)

Removes entry with flag `pseudo_symbol` from the `default_pseudodirs` and `user_defaults.jl` file.
"""
function removedefault_pseudodir(pseudo_symbol::Symbol)
    if isdefined(:default_pseudo_dirs) && haskey(DFControl.default_pseudo_dirs, pseudo_symbol)
        pop!(DFControl.default_pseudo_dirs, pseudo_symbol)
        rm_expr_lhs(default_file, :(default_pseudo_dirs[$(QuoteNode(pseudo_symbol))]))
        if isempty(DFControl.default_pseudo_dirs)
            rm_expr_lhs(default_file, :default_pseudo_dirs)
            default_pseudo_dirs = nothing
        end
    end
    load_defaults(default_file)
end

"""
    setdefault_server(server::String)

Sets the default server variable, and also adds it to the `user_defaults.jl` file.
"""
function setdefault_server(server::String)
    expr_ndef = :(default_server = $server)
    expr_def  = expr_ndef
    define_def(:default_server, expr_ndef, expr_def)
end

"""
    getdefault_server()

Returns the default server if it's defined. If it is not defined return "".
"""
function getdefault_server()
    if isdefined(:default_server)
        return default_server
    else
        return ""
    end
end

"""
    getdefault_pseudo_dirs()

Returns the default pseudo dirs if it's defined. If it is not defined return nothing.
"""
function getdefault_pseudodirs()
    if isdefined(:default_pseudo_dirs)
        return default_pseudo_dirs
    else
        return nothing
    end
end

function getdefault_pseudodir(pseudo_set)
    return getdefault_pseudodirs()[pseudo_set]
end

"""
    configure_defaultpseudos(server = getdefault_server(), pseudo_dirs=getdefault_pseudodirs())

Reads the specified `default_pseudo_dirs` on the `default_server` and sets up the `default_pseudos` variable, and also adds all the entries to the `user_defaults.jl` file.
"""
function configure_defaultpseudos(server = getdefault_server(), pseudo_dirs=getdefault_pseudodirs())
    if server == ""
        error("Either supply a valid server string or setup a default server through 'setdefault_server!()'.")
    end

    if pseudo_dirs == nothing
        error("Either supply valid pseudo directories or setup a default pseudo dir through 'setdefault_pseudodir()'.")
    end

    outputs = Dict{Symbol, String}()
    if typeof(pseudo_dirs) == String
        outputs[:default] = readstring(`ssh -t $server ls $pseudo_dirs`)
    elseif typeof(pseudo_dirs) <: Dict
        for (name, directory) in pseudo_dirs
            outputs[name] = readstring(`ssh -t $server ls $directory`)
        end
    end

    if !isdefined(:default_pseudos)
        expr2file(default_file, :(default_pseudos = Dict{Symbol, Dict{Symbol, Vector{String}}}()))
        init_defaults(default_file)
    end

    elsyms = Symbol[]
    for el in ELEMENTS
        expr2file(default_file,:(default_pseudos[$(QuoteNode(el.symbol))] = Dict{Symbol, Vector{String}}()))
        push!(elsyms, Symbol(titlecase(string(el.symbol))))
    end

    for (name, pseudo_string) in outputs
        pseudos = filter(x -> x != "", split(pseudo_string, "\n"))
        i = 1
        while i <= length(pseudos)
            pseudo  = pseudos[i]
            element = Symbol(titlecase(String(split(split(pseudo, ".")[1], "_")[1])))
            if element in elsyms
                t_expr  = :(String[$pseudo])
                j = 1
                while j + i <= length(pseudos) && Symbol(split(pseudos[i + j],".")[1]) == element
                    push!(t_expr.args,pseudos[i + j])
                    j += 1
                end
                i += j
                expr2file(default_file, :(default_pseudos[$(QuoteNode(element))][$(QuoteNode(name))] = $t_expr))
            else
                i+=1
            end
        end
    end
    load_defaults(default_file)
end

"""
    getdefault_pseudo(atom::Symbol, pseudo_setname=:default; pseudo_specifier=nothing)

Returns the pseudo potential string linked to the atom.
"""
function getdefault_pseudo(atom::Symbol, pseudo_setname=:default; pseudo_specifier="")
    if !isnull(tryparse(Int, String(atom)[end:end]))
        pp_atom = Symbol(String(atom)[1:end-1])
    else
        pp_atom = atom
    end
    if isdefined(:default_pseudos)
        if pseudo_specifier != ""
            return getfirst(x -> contains(x, pseudo_specifier), default_pseudos[pp_atom][pseudo_setname])
        else
            return default_pseudos[pp_atom][pseudo_setname][1]
        end
    end
end

"""
    setdefault_jobheader(lines)

Sets the header that will get added to each job.tt file, if no other header was specified.
"""
function setdefault_jobheader(lines)
    expr = :(default_job_header = $lines)
    expr2file(default_file,expr)
    if !isdefined(:default_job_header)
        init_defaults(default_file)
    else
        load_defaults(default_file)
    end
end

"""
    setdefault_input(input::dfinput, structure, calculation::Symbol)

Adds the input to the `default_inputs` variable, and writes it to a file in user_defaults folder to be read every time on load.
"""
function setdefault_input(input::DFInput{T}, structure, calculation::Symbol) where T
    if !isdefined(:default_inputs)
        expr = :(default_inputs = Dict{Symbol, Tuple{DFInput, Union{AbstractStructure, Void}}}())
        expr2file(default_file,expr)
        init_defaults(default_file)
    end
    filename = dirname(default_file) * "/" * String(calculation)
    runcommand = input.runcommand
    exec        = input.exec
    if T == Wannier90
        save(input, structure, filename * ".win")
        expr2file(default_file, :(default_inputs[$(QuoteNode(calculation))] = read_wannier_input($filename * ".win", runcommand = Exec($(runcommand.exec), $(runcommand.dir), Dict($(runcommand.flags...))), exec = Exec($(exec.exec), $(exec.dir), Dict($(exec.flags...))))))
    elseif T == QE
        save(input, structure, filename * ".in")
        expr2file(default_file, :(default_inputs[$(QuoteNode(calculation))] = read_qe_input($filename * ".in", runcommand = Exec($(runcommand.exec), $(runcommand.dir), Dict($(runcommand.flags...))), exec = Exec($(exec.exec), $(exec.dir), Dict($(exec.flags...))))))
    end
    load_defaults(default_file)
end

"""
    removedefault_input(input::Symbol)

Remove input from the `default_inputs` variable. Also removes the stored input file.
"""
function removedefault_input(input::Symbol)
    if haskey(DFControl.default_inputs, input)
        input_t = pop!(DFControl.default_inputs, input)[1]
        rm_expr_lhs(default_file, :(default_inputs[$(QuoteNode(input))]))
        if isempty(DFControl.default_inputs)
            rm_expr_lhs(default_file, :default_inputs)
            default_inputs = nothing
        end
        rm(joinpath(@__DIR__, "../user_defaults/$(input_t.filename)"))
    else
        error("Default_calculations does not have an input with symbol $symbol.\n  Possible symbols are: $(keys(default_inputs))")
    end
end

function getdefault_jobheader()
    if isdefined(:default_job_header)
        return default_job_header
    else
        return ""
    end
end

function findspecifier(str, strs::Vector{<:AbstractString})
    tmp = Char[]
    i = 1
    for s in strs
        if s == str
            continue
        end
        for (ch1, ch2) in zip(str, s)
            if ch1 != ch2
                push!(tmp, ch1)
            elseif !isempty(tmp)
                break
            end
        end
        !isempty(tmp) && break
        i += 1
    end
    testout = join(tmp)
    for s in strs[i+1:end]
        if contains(s, testout)
            return findspecifier(str, strs[i+1:end])
        end
    end
    return testout
end

function getpseudoset(elsym::Symbol, str::String)
    if !isdefined(:default_pseudos)
        return :none, ""
    else
        for (key, val) in default_pseudos[elsym]
            if length(val) == 1
                str == val[1] && return key, ""
            else
                for v in val
                    if v == str
                        return key, findspecifier(str, val)
                    end
                end
            end
        end
    end
    return :none, ""
end
