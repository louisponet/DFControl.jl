include("qe/fileio.jl")
include("abinit/fileio.jl")
include("wannier90/fileio.jl")
include("elk/fileio.jl")
include("vasp/fileio.jl")

#--------------------Used by other file processing------------------#
#---------------------------BEGINNING GENERAL SECTION-------------------#
#Incomplete: only works with SBATCH right now
#---------------------------END GENERAL SECTION-------------------#

function expr2file(filename::String, expression::Expr)
    eq        = Symbol("=")
    lines     = readlines(filename)
    new_lines = String[]
    found     = false

    if expression.head != eq
        error("For now only writing of assignment expressions is possible.")
    end

    lhs = expression.args[1]
    rhs = expression.args[2]

    for line in lines
        if line == ""
            continue
        end

        expr = Meta.parse(line)
        if typeof(expr) == Nothing
            continue
        end
        if expr.head != eq
            continue
        end

        lhs_t = expr.args[1]
        rhs_t = expr.args[2]

        if lhs_t == lhs
            found = true
            push!(new_lines, "$(:($lhs = $rhs))")
        else
            push!(new_lines, "$expr")
        end
    end

    open(filename, "w") do f
        for line in new_lines
            write(f, line * "\n")
        end
        if !found
            write(f, "$expression\n")
        end
    end
end

function rm_expr_lhs(filename, lhs)
    lines       = readlines(filename)
    write_lines = String[]
    ind_2_rm    = 0

    for line in lines
        lhs_t = Meta.parse(line).args[1]
        if lhs_t == lhs
            continue
        else
            push!(write_lines, line)
        end
    end

    open(filename, "w") do f
        for line in write_lines
            write(f, line * "\n")
        end
    end
end

"LOL this absolutely is impossible to do for QE"
function writeabortfile(job::Job, calculation::Calculation{QE})
    abortpath = joinpath(job.dir, TEMP_CALC_DIR, "$(job.name).EXIT")
    open(abortpath, "w") do f
        return write(f, " \n")
    end
    return abortpath
end

function read_cutoffs_from_pseudofile(file::AbstractString)
    ecutwfc = 0.0
    ecutrho = 0.0
    open(file, "r") do f
        line = readline(f)
        i = 1
        while i < 100 #semi arbitrary cutoff to amount of lines read
            line = readline(f)
            if occursin("Suggested minimum cutoff for wavefunctions:", line)
                ecutwfc = parse(Float64, split(line)[end-1])
                ecutrho = parse(Float64, split(readline(f))[end-1])
                break
            end
            i += 1
        end
    end
    return ecutwfc, ecutrho
end

function write_xsf(filename::AbstractString, structure::DFC.Structure)
    open(filename, "w") do f
        write(f, "CRYSTAL\n")
        c = ustrip.(structure.cell')
        write(f, "PRIMVEC\n")
        write(f, "$(c[1,1]) $(c[1,2]) $(c[1,3])\n")
        write(f, "$(c[2,1]) $(c[2,2]) $(c[2,3])\n")
        write(f, "$(c[3,1]) $(c[3,2]) $(c[3,3])\n")
        write(f, "CONVVEC\n")
        write(f, "$(c[1,1]) $(c[1,2]) $(c[1,3])\n")
        write(f, "$(c[2,1]) $(c[2,2]) $(c[2,3])\n")
        write(f, "$(c[3,1]) $(c[3,2]) $(c[3,3])\n")
        write(f, "PRIMCOORD\n")
        write(f, "$(length(structure.atoms)) 1\n")
        for at in structure.atoms
            n = at.element.symbol
            p = ustrip.(at.position_cart)
            write(f, "$n $(p[1]) $(p[2]) $(p[3])\n")
        end
    end
end

function writelines(file, lines)
    open(file, "w") do f
        for l in lines
            write(f, l)
            write(f, "\n")
        end
    end
end
