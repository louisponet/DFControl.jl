# Calculation utils
using DFControl: inpath, outpath
#

function readoutput(c::Calculation{QE}; kwargs...)
    return qe_read_output(c; kwargs...)
end

#asserts
function pdos(c::Calculation{QE}, atsym::Symbol, magnetic::Bool, soc::Bool,
              filter_word = "")
    @assert isprojwfc(c) "Please specify a valid projwfc calculation."
    kresolved = hasflag(c, :kresolveddos) && calculation[:kresolveddos]
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

