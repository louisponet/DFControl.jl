const RUN_EXECS = ["mpirun", "mpiexec", "srun"]

isparseable(exec::Exec) = exec.exec ∈ parseable_execs()

hasflag(exec::Exec, s::Symbol) = findfirst(x -> x.symbol == s, exec.flags) != nothing

Base.hash(e::ExecFlag, h::UInt) = hash(e.symbol, hash(e.value, h))
Base.:(==)(e1::ExecFlag, e2::ExecFlag) = e1.symbol == e2.symbol && e1.value == e2.value
Base.eltype(::ExecFlag{T}) where {T} = T

function Base.:(==)(e1::Exec, e2::Exec)
    if e1.exec != e2.exec || e1.dir != e2.dir
        return false
    else
        for f in e1.flags
            f2 = getfirst(isequal(f), e2.flags)
            f2 === nothing && return false
            f == f2 && return true
        end
        return true
    end
end

allexecs() = vcat(RUN_EXECS, QE_EXECS, WAN_EXECS, ELK_EXECS)
parseable_execs() = vcat(QE_EXECS, WAN_EXECS, ELK_EXECS)
has_parseable_exec(l::String) = occursin(">", l) && any(occursin.(parseable_execs(), (l,)))

function calculationparser(exec::Exec)
    if is_qe_exec(exec)
        qe_read_calculation
    elseif is_wannier_exec(exec)
        wan_read_calculation
    elseif is_elk_exec(exec)
        elk_read_calculation
    end
end

function set_flags!(exec::Exec, flags...)
    for (f, val) in flags
        flag = isa(f, String) ? getfirst(x -> x.name == f, exec.flags) :
               getfirst(x -> x.symbol == f, exec.flags)
        if flag != nothing
            flag.value = convert(eltype(flag), val)
        else
            found = false
            #TODO generalize this
            for (f1, f2) in zip((is_qe_exec, is_wannier_exec, is_mpi_exec),
                    (qe_execflag, wan_execflag, mpi_flag))
                if f1(exec)
                    def_flag = f2(f)
                    if def_flag !== nothing
                        push!(exec.flags, ExecFlag(def_flag, val))
                        found = true
                        break
                    end
                end
            end
            if !found
                error("Flag $f was not found in allowed executable flags.")
            end
        end
    end
    return exec.flags
end

"""
    set_execflags!(calculation::DFCalculation, exec::String, flags::Pair...)
    set_execflags!(job::DFJob, exec::String, flags::Pair...)

Goes through the execs in `calculation` or `job` and sets the `flags`
if the executable name contains `exec`.
"""
function set_execflags!(calculation::DFCalculation, exec::String, flags::Pair...)
    for e in execs(calculation, exec)
        set_flags!(e, flags...)
    end
end

function set_execflags!(job::DFJob, exec::String, flags::Pair...)
    for i in job.calculations
        set_execflags!(i, exec, flags...)
    end
end

function rm_flags!(exec::Exec, flags...)
    for f in flags
        if isa(f, String)
            filter!(x -> x.name != f, exec.flags)
        elseif isa(f, Symbol)
            filter!(x -> x.symbol != f, exec.flags)
        end
    end
    return exec.flags
end

"""
    rm_execflags!(calculation::DFCalculation, exec::String, flags...)
    rm_execflags!(job::DFJob, exec::String, flags...)

Goes through the execs in `calculation` or `job` and removes the `flags`
if the executable name contains `exec`.
"""
function rm_execflags!(calculation::DFCalculation, exec::String, flags...)
    return rm_flags!.(execs(calculation, exec), flags...)
end

function rm_execflags!(job::DFJob, exec, flags...)
    return rm_execflags!.(job.calculations, (exec, flags)...)
end

set_dir!(exec::Exec, dir) = exec.dir = dir

"""
    set_execdir!(calculation::DFCalculation, exec::String, dir::String)
    set_execdir!(job::DFJob, exec::String, dir::String)
    
Goes through the execs in `calculation` or `job` and sets the `dir`
if the executable name contains `exec`.

Example:
```julia
set_execdir!(calculation, "pw.x", "/path/to/QE/bin")
set_execdir!(job, "pw.x", "/path/to/QE/bin")
```
"""
function set_execdir!(calculation::DFCalculation, exec::String, dir::String)
    return set_dir!.(execs(calculation, exec), dir)
end

function set_execdir!(job::DFJob, exec::String, dir::String)
    return set_execdir!.(job.calculations, exec, dir)
end

include("qe/execs.jl")
include("wannier90/execs.jl")
include("elk/execs.jl")

#### MPI Exec Functionality ###
include(joinpath(depsdir, "mpirunflags.jl"))
const MPI_FLAGS = _MPI_FLAGS()

mpi_flag(flag::AbstractString) = getfirst(x -> x.name == flag, MPI_FLAGS)
mpi_flag(flag::Symbol) = getfirst(x -> x.symbol == flag, MPI_FLAGS)

function mpi_flag_val(::Type{String}, line, i)
    v = line[i+1]
    return v, i + 2
end

function mpi_flag_val(::Type{Vector{String}}, line, i)
    tval = String[]
    while i + 1 <= length(line) && !occursin('-', line[i+1])
        push!(tval, line[i+1])
        i += 1
    end
    return tval, i + 1
end

function mpi_flag_val(::Type{Int}, line, i)
    if line[i+1][1] == '$'
        tval = line[i+1]
    else
        tval = parse(Int, line[i+1])
    end
    return tval, i + 2
end

function mpi_flag_val(::Type{Vector{Int}}, line, i)
    tval = Union{Int,String}[]
    while i + 1 <= length(line) && !occursin('-', line[i+1])
        if line[i+1][1] == '$'
            push!(tval, line[i+1])
        else
            push!(tval, parse(Int, line[i+1]))
        end
        i += 1
    end
    return tval, i + 1
end

function mpi_flag_val(::Type{Pair{Int,String}}, line, i)
    tval = Union{Int,String}[]
    while i + 1 <= length(line) && !occursin('-', line[i+1])
        if line[i+1][1] == '$'
            push!(tval, line[i+1])
        else
            tparse = tryparse(Int, line[i+1])
            if tparse != nothing
                push!(tval, tparse)
            else
                push!(tval, string(line[i+1]))
            end
        end
        i += 1
    end
    return tval, i + 1
end
function parse_mpi_flags(line::Vector{<:SubString})
    eflags = ExecFlag[]
    i = 1
    while i <= length(line)
        s = line[i]
        if s[1] != '-'
            break
        end
        if s[2] == '-'
            mflag = mpi_flag(strip(s, '-'))
        else
            mflag = mpi_flag(Symbol(strip(s, '-')))
        end
        @assert mflag !== nothing "$(strip(s, '-')) is not a recognized mpi_flag"
        val, i = mpi_flag_val(eltype(mflag), line, i)
        push!(eflags, ExecFlag(mflag, val))
    end
    return eflags
end

is_mpi_exec(exec::Exec) = occursin("mpi", exec.exec)

"""
    execs(calculation::DFCalculation, exec::String="")
    execs(job::DFJob, exec::String="")

Convenience function to filter all executables in `calculation` or `job` for which `exec` occurs in the executable name.
"""
function execs(calculation::DFCalculation, exec::String="")
    return filter(x -> occursin(exec, x.exec), calculation.execs)
end
function execs(job::DFJob, exec::String="")
    return filter(x -> occursin(exec, x.exec), vcat(execs.(calculations(job))...))
end

"""
    exec(calculation::DFCalculation, exec::String)
    exec(job::DFJob, exec::String)

Convenience function that returns the first executable in `calculation` or `job`
for which `exec` occursin the name.
"""
exec(calculation::DFCalculation, exec::String) = (es = execs(calculation, exec); isempty(es) ? nothing : es[1])
exec(job::DFJob, exec::String) = (es = execs(job, exec); isempty(es) ? nothing : es[1])

function Base.string(e::Exec)
    direxec = joinpath(e.dir, e.exec)
    str = "$direxec"
    for flag in e.flags
        str *= " $(join(fill('-', flag.minus_count)))$(flag.symbol)"
        if !isa(flag.value, AbstractString)
            for v in flag.value
                str *= " $v"
            end
        else
            str *= " $(flag.value)"
        end
    end
    return str
end

path(exec::Exec) = joinpath(exec.dir, exec.exec)

function verify_execs(job::DFJob, server::Server)
    proc = establish_connection(server)
    for e in unique(execs(job))
        if !isempty(e.dir)
            p = path(e)
            if !Distributed.remotecall_fetch(ispath, proc, p)
                error("$e is not a valid executable on server $(server.name)")
            end
        else
            if Distributed.remotecall_fetch(Sys.which, proc, e.exec) === nothing
                error("$e is not a valid executable on server $(server.name)")
            end
        end
    end
    Distributed.rmprocs(proc)
end
