issoccalc(calculation::Calculation{Wannier90}) = flag(calculation, :spinors) == true

function readoutput(calculation::Calculation{Wannier90}; kwargs...)
    return wan_read_output(outpath(calculation); kwargs...)
end

for f in (:cp, :mv)
    @eval function Base.$f(i::Calculation{Wannier90}, dest::String; kwargs...)
        for glob in ("$(i.name)","UNK") # seedname should also cover generated pw2wannier90 files
            for f in searchdir(i, glob)
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
    end
end

