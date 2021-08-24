const EXEC_FILE = DFC.config_path("registered_execs.json")

function load_execs()
    return ispath(EXEC_FILE) ? JSON3.read(read(EXEC_FILE, String), Vector{Exec}) : Exec[]
end

function write_execs(execs)
    return JSON3.write(EXEC_FILE, execs)
end

function isrunnable(e::Exec)
    # To find the path to then ldd on
    fullpath = !isempty(e.dir) ? joinpath(e.dir, e.exec) : Sys.which(e.exec)
    if fullpath !== nothing && ispath(fullpath)
        out = Pipe()
        err = Pipe()
        #by definition by submitting something with modules this means the server has them
        if !isempty(e.modules)
            cmd = `echo "source /etc/profile && module load $(join(e.modules, " ")) && ldd $fullpath"`
        else
            cmd = `echo "source /etc/profile && ldd $fullpath"`
        end
        run(pipeline(pipeline(cmd, stdout=ignorestatus(`bash`)),stdout = out, stderr=err))
        close(out.in)
        close(err.in)

        stderr = String(read(err))
        stdout = String(read(out))
        # This basically means that the executable would run
        return !occursin("not found", stdout) && !occursin("not found", stderr)
    else
        return false
    end
end

"Verifies the validity of an executable and also registers it if it wasn't known before."
function verify_exec(e::Exec)
    valid = isrunnable(e)
    if valid
        maybe_register(e)
    end
    return valid
end

function maybe_register(e::Exec)
    known_execs = load_execs()
    preexisting = getfirst(x-> x.dir == e.dir && x.exec == e.exec,  known_execs)
    if preexisting === nothing
        empty!(e.flags)
        push!(known_execs, e)
    elseif length(preexisting.modules) > length(e.modules)
        preexisting.modules = e.modules
    end
    write_execs(known_execs)
end

"Finds all executables that are known with the same exec name."
function known_execs(exec::AbstractString)
    out = Exec[] 
    for e in load_execs()
        if e.exec == exec
            push!(out, e)
        end
    end
    return out
end

#these are all the control data, they hold the flags that guide the calculation
function pdos(calculation::Calculation, args...)
    @error "pdos reading not implemented for package $(eltype(calculation))."
end

function Emin_from_projwfc(calculation::Calculation, args...)
    @error "Emin_from_projwfc is not implemented for package $(eltype(calculation))."
end

function readoutput(c::Calculation; kwargs...)
    @error "Output parsing for package $(eltype(c)) not implemented."
end

"""
    outputdata(calculation::Calculation; extra_parse_funcs=[], print=true, overwrite=true)

If an output file exists for `calculation` this will parse it and return a `Dict` with the parsed data.
If `overwrite=false` and `calculation.outputdata` is not empty, this will be returned instead of reparsing the
output file.

`extra_parse_funcs` should be `Vector{Pair{String,Function}}`, where the string will be used to match `occursin(str, line)`
for each line of the output file. If this returns `true`, the function will be called as `func(results_dict, line, file)`.
The purpose is to allow for additional parsing that is not implemented, or for temporary non-standard values that are printed
while working on the DFT code, e.g. for debugging.

!!! note
    This only works for files where not all information is pulled out for yet,
    e.g. projwfc.x outputs are fully parsed already.

Example (from src/qe/fileio.jl):
```julia

function qe_parse_nat(results, line, f)
    results[:nat] = parse(Int, split(line)[end])
end

outputdata(job["scf"], extra_parse_funcs = ["number of atoms/cell" => qe_parse_nat])
```
"""
function outputdata(calculation::Calculation;
                    extra_parse_funcs::Vector{<:Pair{String}} = Pair{String}[],
                    print = true, overwrite = true)
    if ispath(Calculations.outpath(calculation)) 
        t = FileIO.readoutput(calculation; parse_funcs = extra_parse_funcs)
        return t === nothing ?
                              FileIO.parse_file(Calculations.outpath(calculation),
                                                extra_parse_funcs) : t
    end
    print &&
        (@warn "No output data or output file found for calculation: $(calculation.name).")
    return Dict{Symbol,Any}()
end

