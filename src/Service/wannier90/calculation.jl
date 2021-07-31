issoccalc(calculation::DFCalculation{Wannier90}) = flag(calculation, :spinors) == true

function readoutput(calculation::DFCalculation{Wannier90}; kwargs...)
    return wan_read_output(outpath(calculation); kwargs...)
end

for f in (:cp, :mv)
    @eval function Base.$f(i::DFCalculation{Wannier90}, dest::String; kwargs...)
        for glob in ("$(name(i))", "UNK") # seedname should also cover generated pw2wannier90 files
            for f in searchdir(i, glob)
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
    end
end

