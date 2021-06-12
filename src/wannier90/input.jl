issoccalc(input::DFInput{Wannier90}) = flag(input, :spinors) == true

readoutput(input::DFInput{Wannier90}) = wan_read_output(outpath(input))

for f in (:cp, :mv)
    @eval function Base.$f(i::DFInput{Wannier90}, dest::String; kwargs...)
        for glob in ("$(name(i))", "UNK") # seedname should also cover generated pw2wannier90 files
            for f in searchdir(i, glob)
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
    end
end

