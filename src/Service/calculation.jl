#these are all the control data, they hold the flags that guide the calculation
outfiles(c::Calculation) = filter(ispath, [Calculations.outpath(c)])

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
    if ispath(Calculations.outpath(calculation)) || !isempty(calculation.outdata)
        if !overwrite && !isempty(calculation.outdata)
            return calculation.outdata
        else
            t = FileIO.readoutput(calculation; parse_funcs = extra_parse_funcs)
            calculation.outdata = t === nothing ?
                                  FileIO.parse_file(Calculation.outpath(calculation), extra_parse_funcs) : t
            return calculation.outdata
        end
    end
    print &&
        (@warn "No output data or output file found for calculation: $(calculation.name).")
    return Dict{Symbol,Any}()
end

function pdos(c::Calculation{QE}, atsym::Symbol, magnetic::Bool, soc::Bool,
              filter_word = "")
    @assert isprojwfc(c) "Please specify a valid projwfc calculation."
    kresolved = haskey(c, :kresolveddos) && calculation[:kresolveddos]
    files = filter(x -> occursin("($atsym)", x) &&
                            occursin("#", x) &&
                            occursin(filter_word, x), searchdir(c.dir, "pdos"))
    @assert !isempty(files) "No pdos files found in calculation directory $(c.dir)"
    files = joinpath.((c,), files)
    energies, = kresolved ? qe_read_kpdos(files[1]) : qe_read_pdos(files[1])
    atdos = magnetic && !soc ? zeros(size(energies, 1), 2) : zeros(size(energies, 1))
    if kresolved
        for f in files
            if magnetic && !occursin(".5", f)
                tu = qe_read_kpdos(f, 2)[2]
                td = qe_read_kpdos(f, 3)[2]
                atdos[:, 1] .+= reduce(+, tu; dims = 2) ./ size(tu, 2)
                atdos[:, 2] .+= reduce(+, td; dims = 2) ./ size(tu, 2)
                # elseif occursin(".5", f)
            else
                t = qe_read_kpdos(f, 1)[2]
                atdos .+= (reshape(reduce(+, t; dims = 2), size(atdos, 1)) ./ size(t, 2))
            end
        end
    else
        for f in files
            if magnetic && !occursin(".5", f)
                atdos .+= qe_read_pdos(f)[2][:, 1:2]
                # elseif occursin(".5", f)
            else
                atdos .+= qe_read_pdos(f)[2][:, 1]
            end
        end
    end
    return (energies = energies, pdos = atdos)
end

issoccalc(calculation::Calculation{Wannier90}) = flag(calculation, :spinors) == true

for f in (:cp, :mv)
    @eval function Base.$f(i::Calculation{Wannier90}, dest::String; kwargs...)
        for glob in ("$(i.name)","UNK") # seedname should also cover generated pw2wannier90 files
            for f in searchdir(i, glob)
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
    end
end

