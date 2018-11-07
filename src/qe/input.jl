import ..DFControl: DFJob
import ..DFControl: infile, outfile, setkpoints!, isbandscalc, isnscfcalc, isscfcalc, isspincalc, readoutput, generate_jobinputs, sanitizeflags!, cleanflags!, fortstring

outfile(input::DFInput{QE})        = namewext(input, ".out")
infile(input::DFInput{QE})         = namewext(input, ".in")

kgrid(na, nb, nc) = reshape([[a, b, c, 1 / (na * nb * nc)] for a in collect(range(0, stop=1, length=na + 1))[1:end - 1], b in collect(range(0, stop=1, length=nb + 1))[1:end - 1], c in collect(range(0, stop=1, length=nc + 1))[1:end - 1]], (na * nb * nc))

function setkpoints!(input::DFInput{QE}, k_grid::NTuple{3, Int}; print=true) #nscf

    calc = flag(input, :calculation)
    print && calc != "'nscf'" && warn("Expected calculation to be 'nscf'.\nGot $calc.")
    setdata!(input, :k_points, kgrid(k_grid...), option = :crystal, print=print)
    prod(k_grid) > 100 && setflags!(input, :verbosity => "'high'", print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::NTuple{6, Int}; print=true) #scf
    calc = flag(input, :calculation)
    print && (calc != "'scf'" || !occursin("relax", calc)) && warn("Expected calculation to be 'scf', 'vc-relax', 'relax'.\nGot $calc.")
    setdata!(input, :k_points, [k_grid...], option = :automatic, print=print)
    prod(k_grid[1:3]) > 100 && setflags!(input, :verbosity => "'high'", print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::Vector{NTuple{4, T}}; print=true, k_option=:crystal_b) where T<:AbstractFloat
    calc = flag(input, :calculation)
    print && calc != "'bands'" && warn("Expected calculation to be 'bands', got $calc.")
    @assert in(k_option, [:tpiba_b, :crystal_b, :tpiba_c, :crystal_c]) error("Only $([:tpiba_b, :crystal_b, :tpiba_c, :crystal_c]...) are allowed as a k_option, got $k_option.")
    if k_option in [:tpiba_c, :crystal_c]
        @assert length(k_grid) == 3 error("If $([:tpiba_c, :crystal_c]...) is selected the length of the k_points needs to be 3, got length: $(length(k_grid)).")
    end
    num_k = 0.0
    for k in k_grid
        num_k += k[4]
    end
    if num_k > 100.
        setflags!(input, :verbosity => "'high'", print=print)
        if print
            @info "Verbosity is set to high because num_kpoints > 100,\n
                       otherwise bands won't get printed."
        end
    end
    setdata!(input, :k_points, k_grid, option=k_option, print=print)
    return input
end

function sanitizeflags!(input::DFInput{QE}, job::DFJob)
    cleanflags!(input)

    setflags!(input, :outdir => fortstring(job.server_dir), print=false)
    setflags!(input, :prefix => job.name, print=false)
    flag(input, :ecutwfc) == nothing && setflags!(input,  :ecutwfc => 25.0) #arbitrary default
    #TODO add all the required flags
end

isbandscalc(input::DFInput{QE}) = flag(input, :calculation) == "'bands'"
isnscfcalc(input::DFInput{QE}) = flag(input, :calculation) == "'nscf'"
isscfcalc(input::DFInput{QE}) = flag(input, :calculation) == "'scf'"
isspincalc(input::DFInput{QE}) = all(flag(input, :nspin) .!= [nothing, 1])

readoutput(input::DFInput{QE}) = read_output(outpath(input))

function generate_jobinputs(::Type{QE}, local_dir, structure, calculations, common_flags...)
    job_calcs = DFInput{QE}[]
    if typeof(common_flags) != Dict
        common_flags = Dict(common_flags)
    end
    for (calc, (excs, data)) in calculations
        calc_ = typeof(calc) == String ? Symbol(calc) : calc
        if in(calc_, [Symbol("vc-relax"), :relax, :scf])
            k_points = get(data, :k_points, [1, 1, 1, 0, 0, 0])
            k_option = :automatic
        elseif calc_ == :nscf
            k_points = kgrid(get(data, :k_points, [1, 1, 1])...)
            k_option = :crystal
        elseif calc_ == :bands
            k_points = get(data, :k_points, [[0., 0., 0., 1.]])
            num_k = 0.0
            for point in k_points
                num_k += point[4]
            end
            if num_k > 100.
                push!(data[:flags], :verbosity => "'high'")
            end
            k_option = :crystal_b
        end
        flags  = get(data, :flags, Dict{Symbol, Any}())
        if excs[2].exec == "pw.x"
            push!(flags, :calculation => "'$(string(calc_))'")
            datablocks = [InputData(:k_points, k_option, k_points)]
        else
            datablocks =  InputData[]
        end
        input_ = DFInput{QE}(string(calc_), local_dir,
                         Dict{Symbol, Any}(),
                         datablocks, excs, true)
        setflags!(input_, flags..., print=false)
        push!(job_calcs, input_)
    end
    return job_calcs
end
