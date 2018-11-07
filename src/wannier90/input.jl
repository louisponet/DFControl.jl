import ..DFControl: readoutput, setkpoints!, infile, outfile, generate_waninputs
import ..DFControl.QuantumEspresso: QE


readoutput(input::DFInput{Wan90}) = SymAnyDict()
infile(input::DFInput{Wan90})  = namewext(input, ".win")
outfile(input::DFInput{Wan90}) = namewext(input, ".wout")

"""
    setkpoints!(input::DFInput, k_grid)

Sets the kpoints of the input. Will automatically generate the kgrid values if necessary.
"""
function setkpoints!(input::DFInput{Wan90}, k_grid::NTuple{3, Int}; print=true)
    setflags!(input, :mp_grid => [k_grid...], print=print)
    setdata!(input, :kpoints, kgrid(k_grid..., :wan), print=print)
    return input
end

kgrid(na, nb, nc) = reshape([[a, b, c] for a in collect(range(0, stop=1, length=na + 1))[1:end - 1], b in collect(range(0, stop=1, length=nb + 1))[1:end - 1], c in collect(range(0, stop=1, length=nc + 1))[1:end - 1]],(na * nb * nc))

"""
    generate_waninputs(nscf::DFInput{QE}, wanflags...;
                     pw2wanexec=Exec("pw2wannier90.x", nscf.execs[2].dir),
                     wanexec=Exec("wannier90.x", nscf.execs[2].dir),
                     bands=readbands(nscf),
                     print=true)

Generates the input files for Wannier90 starting from an already performed nscf calculation.
"""
function generate_waninputs(nscf::DFInput{QE}, wanflags...;
                     pw2wanexec=Exec("pw2wannier90.x", nscf.execs[2].dir),
                     wanexec=Exec("wannier90.x", nscf.execs[2].dir),
                     print=true)

    spin = isspincalc(nscf)
    if spin
        pw2wannames = ["pw2wan_up", "pw2wan_dn"]
        wannames = ["wanup", "wandn"]
        print && info("Spin polarized calculation found (inferred from nscf input).")
    else
        pw2wannames = ["pw2wan"]
        wannames = ["wan"]
    end

    if flag(nscf, :nosym) != true
        print && info("'nosym' flag was not set in the nscf calculation.\nIf this was not intended please set it and rerun the nscf calculation.\nThis generally gives errors because of omitted kpoints, needed for pw2wannier90.x")
    end

    wanflags = SymAnyDict(wanflags)
    wanflags[:mp_grid] = kakbkc(data(nscf, :k_points).data)
    print && info("mp_grid=$(join(wanflags[:mp_grid]," ")) (inferred from nscf input).")

    pw2wanflags = SymAnyDict(:prefix => flag(nscf, :prefix), :outdir => flag(nscf, :outdir) == nothing ? "'./'" : flag(nscf, :outdir))
    if haskey(wanflags, :write_hr)
        pw2wanflags[:write_amn] = true
        pw2wanflags[:write_mmn] = true
    end
    if haskey(wanflags, :wannier_plot)
        pw2wanflags[:write_unk] = true
    end

    kdata = InputData(:kpoints, :none, kgrid(wanflags[:mp_grid]...))
    inputs = DFInput[]
    for (pw2wanfil, wanfil) in zip(pw2wannames, wannames)
        push!(inputs, DFInput{Wan90}(wanfil, job.local_dir, copy(wanflags), [kdata], [Exec(), wanexec], true))
        push!(inputs, DFInput{QE}(pw2wanfil, job.local_dir, copy(pw2wanflags), InputData[], [nscf.execs[1], pw2wanexec], true))
    end

    setflags!(n, flags...) = setflags!(findfirst(x-> name(x) == n, inputs), flags...)
    if spin
        setflags!("pw2wan_up",  :spin_component => "'up'")
        setflags!("pw2wan_dn",  :spin_component => "'down'")
        setflags!("wanup", :spin => "'up'")
        setflags!("wandn", :spin => "'down'")
    end
    return inputs
end
