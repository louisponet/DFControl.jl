import Base: parse

#this is all pretty hacky with regards to the new structure and atom api. can for sure be a lot better!
"Quantum espresso card option parser"
function cardoption(line)
    sline = split(line)
    if length(sline) < 2 && lowercase(sline[1]) == "k_points"
        return :tpiba
    else
        return Symbol(match(r"((?:[a-z][a-z0-9_]*))", sline[2]).match)
    end
end

function qe_parse_time(str::AbstractString)
    s = findfirst(isequal('s'), str)
    ms = findfirst(isequal('.'), str)
    m = findfirst(isequal('m'), str)
    m = m === nothing ? 0 : m
    h = findfirst(isequal('h'), str)
    h = h === nothing ? 0 : h
    t = Millisecond(0)
    if s !== nothing
        t += Second(parse(Int, str[m+1:ms-1])) + Millisecond(parse(Int, str[ms+1:s-1]))
    end
    if m != 0
        t += Minute(parse(Int, str[h+1:m-1]))
    end
    if h != 0
        t += Hour(parse(Int, str[1:h-1]))
    end
    return t
end

function qe_read_output(calculation::DFCalculation{QE}, args...; kwargs...)
    if isprojwfc(calculation)
        return qe_read_projwfc_output(calculation, args...; kwargs...)
    elseif ishp(calculation)
        return qe_read_hp_output(calculation, args...; kwargs...)
    elseif ispw(calculation)
        return qe_read_pw_output(outpath(calculation), args...; kwargs...)
    end
end

function parse_Hubbard_block(f)
    # Each of these will have n Hubbard typ elements at the end
    ids = Int[]
    traces = NamedTuple{(:up, :down, :total),NTuple{3,Float64}}[]
    eigvals = (up = Vector{Float64}[], down = Vector{Float64}[])
    eigvec = (up = Matrix{Float64}[], down = Matrix{Float64}[])
    occupations = (up = Matrix{Float64}[], down = Matrix{Float64}[])
    magmoms = Float64[]
    line = readline(f)
    cur_spin = :up
    while strip(line) != "--- exit write_ns ---"
        line = readline(f)
        if line[1:4] == "atom"
            sline = split(line)
            push!(ids, parse(Int, sline[2]))
            push!(traces,
                  NamedTuple{(:up, :down, :total)}(parse.(Float64,
                                                          (sline[end-2], sline[end-1],
                                                           sline[end]))))
            for spin in (:up, :down)
                readline(f) #should be spin1
                readline(f)# should be eigvals
                push!(eigvals[spin], parse.(Float64, split(readline(f))))
                dim = length(eigvals[spin][1])
                readline(f) #eigvectors
                tmat = zeros(dim, dim)
                for i in 1:dim
                    tmat[i, :] = parse.(Float64, split(readline(f)))
                end
                push!(eigvec[spin], tmat)
                readline(f) #occupations
                for i in 1:dim
                    tmat[i, :] = parse.(Float64, split(readline(f)))
                end
                push!(occupations[spin], tmat)
            end
            push!(magmoms, parse(Float64, split(readline(f))[end]))
        end
    end
    return [(id = i, trace = t, eigvals = (up = val_up, down = val_down),
             eigvecs = (up = vec_up, down = vec_down),
             occupations = (up = occ_up, down = occ_down), magmom = m)
            for (i, t, val_up, val_down, vec_up, vec_down, occ_up, occ_down, m) in
                zip(ids, traces, eigvals.up, eigvals.down, eigvec.up, eigvec.down,
                    occupations.up, occupations.down, magmoms)]
end

function qe_parse_polarization(out, line, f)
    s_line = split(line)
    P      = parse(Float64, s_line[3])
    mod    = parse(Float64, s_line[5][1:end-1])
    readline(f)
    s_line             = parse.(Float64, split(readline(f))[6:2:10])
    out[:polarization] = Point3{Float64}(P * s_line[1], P * s_line[2], P * s_line[3])
    return out[:pol_mod] = mod
end

function qe_parse_lattice_parameter(out, line, f)
    out[:in_alat] = ustrip(uconvert(Ang, parse(Float64, split(line)[5]) * 1a₀))
    return out[:alat] = :crystal
end

function qe_parse_n_KS(out, line, f)
    return out[:n_KS_states] = parse(Int, split(line)[5])
end

function qe_parse_crystal_axes(out, line, f)
    m = Mat3(reshape([parse.(Float64, split(readline(f))[4:6]);
                      parse.(Float64, split(readline(f))[4:6]);
                      parse.(Float64, split(readline(f))[4:6])], (3, 3))')
    out[:cell_parameters] = copy(m)
    return out[:in_cell] = m
end

function qe_parse_reciprocal_axes(out, line, f)
    cell_1 = parse.(Float64, split(readline(f))[4:6]) .* 2π / out[:in_alat]
    cell_2 = parse.(Float64, split(readline(f))[4:6]) .* 2π / out[:in_alat]
    cell_3 = parse.(Float64, split(readline(f))[4:6]) .* 2π / out[:in_alat]
    return out[:in_recip_cell] = Mat3([cell_1 cell_2 cell_3])
end
function qe_parse_atomic_species(out, line, f)
    if !haskey(out, :atsyms)
        line = readline(f)
        out[:atsyms] = Symbol[]
        while !isempty(line)
            push!(out[:atsyms], Symbol(strip_split(line)[1]))
            line = readline(f)
        end
    end
end

function qe_parse_nat(out, line, f)
    return out[:nat] = parse(Int, split(line)[end])
end

function qe_parse_crystal_positions(out, line, f)
    readline(f)
    readline(f)
    out[:in_cryst_positions] = Tuple{Symbol,Point3{Float64}}[] # in crystal coord
    for i in 1:out[:nat]
        sline = split(readline(f))
        push!(out[:in_cryst_positions],
              (Symbol(sline[2]),
               Point3(parse(Float64, sline[7]), parse(Float64, sline[8]),
                      parse(Float64, sline[9]))))
    end
end

function qe_parse_cart_positions(out, line, f)
    readline(f)
    readline(f)
    out[:in_cart_positions] = Tuple{Symbol,Point3{Float64}}[] # in crystal coord
    for i in 1:out[:nat]
        sline = split(readline(f))
        push!(out[:in_cart_positions],
              (Symbol(sline[2]),
               Point3(parse(Float64, sline[7]), parse(Float64, sline[8]),
                      parse(Float64, sline[9]))))
    end
end

function qe_parse_pseudo(out, line, f)
    !haskey(out, :pseudos) && (out[:pseudos] = Dict{Symbol,Pseudo}())
    pseudopath = readline(f) |> strip |> splitdir
    return out[:pseudos][Symbol(split(line)[5])] = Pseudo(pseudopath[2], pseudopath[1])
end

function qe_parse_fermi(out, line, f)
    sline = split(line)
    if occursin("energy is", line)
        out[:fermi] = parse(Float64, sline[5])
    elseif occursin("up/dw", line)
        sline            = split(line)
        out[:fermi_up]   = parse(Float64, sline[7])
        out[:fermi_down] = parse(Float64, sline[8])
        out[:fermi]      = min(out[:fermi_down], out[:fermi_up])
    end
end

function qe_parse_highest_lowest(out, line, f)
    sline = split(line)
    if occursin(line, "lowest")
        high = parse(Float64, sline[7])
        low = parse(Float64, sline[8])
        out[:fermi] = high
        out[:highest_occupied] = high
        out[:lowest_unoccupied] = low
    else
        out[:fermi] = parse(Float64, sline[5])
        out[:highest_occupied] = out[:fermi]
    end
end

function qe_parse_total_energy(out, line, f)
    if haskey(out, :total_energy)
        push!(out[:total_energy], parse(Float64, split(line)[end-1]))
    else
        out[:total_energy] = [parse(Float64, split(line)[end-1])]
    end
end

function qe_parse_k_cryst(out, line, f)
    if length(split(line)) == 2
        out[:k_cryst] = (v = Vec3{Float64}[], w = Float64[])
        line = readline(f)
        while line != "" && !occursin("--------", line)
            parsed = parse_k_line(line, Float64)
            push!(out[:k_cryst].v, parsed.v)
            push!(out[:k_cryst].w, parsed.w)
            line = readline(f)
        end
    end
end

function qe_parse_k_cart(out, line, f)
    if length(split(line)) == 5
        line = readline(f)
        alat = out[:in_alat]
        out[:k_cart] = (v = Vec3{typeof(2π / alat)}[], w = Float64[])
        while line != "" && !occursin("--------", line)
            tparse = parse_k_line(line, Float64)
            push!(out[:k_cart].v, tparse.v .* 2π / alat)
            push!(out[:k_cart].w, tparse.w)
            line = readline(f)
        end
    end
end

function qe_parse_k_eigvals(out, line, f)
    tmp = Float64[]
    readline(f)
    line = readline(f)
    while line != "" && !occursin("--------", line)
        append!(tmp, parse_line(Float64, line))
        line = readline(f)
    end
    if haskey(out, :k_eigvals)
        push!(out[:k_eigvals], tmp)
    else
        out[:k_eigvals] = [tmp]
    end
end

#! format: off
function qe_parse_cell_parameters(out, line, f)
    out[:alat]            = occursin("angstrom", line) ? :angstrom : parse(Float64, split(line)[end][1:end-1])
    out[:cell_parameters] = Mat3(reshape([parse.(Float64, split(readline(f)));
                                          parse.(Float64, split(readline(f)));
                                          parse.(Float64, split(readline(f)))], (3, 3))')
end
#! format: on
function qe_parse_atomic_positions(out, line, f)
    out[:pos_option] = cardoption(line)
    line = readline(f)
    atoms = Tuple{Symbol,Point3{Float64}}[]
    while length(atoms) < out[:nat]
        s_line = split(line)
        key    = Symbol(s_line[1])
        push!(atoms, (key, Point3(parse.(Float64, s_line[2:end])...)))
        line = readline(f)
    end
    return out[:atomic_positions] = atoms
end

function qe_parse_total_force(out, line, f)
    sline = split(line)
    force = parse(Float64, sline[4])
    scf_contrib = parse(Float64, sline[end])
    if !haskey(out, :total_force)
        out[:total_force] = [force]
    else
        push!(out[:total_force], force)
    end
    if !haskey(out, :scf_correction)
        out[:scf_correction] = [scf_contrib]
    else
        push!(out[:scf_correction], scf_contrib)
    end
end

function qe_parse_scf_iteration(out, line, f)
    sline = split(line)
    it = length(sline[2]) == 1 ? parse(Int, sline[3]) :
         sline[2][2:end] == "***" ? out[:scf_iteration][end] + 1 :
         parse(Int, sline[2][2:end])
    if !haskey(out, :scf_iteration)
        out[:scf_iteration] = [it]
    else
        push!(out[:scf_iteration], it)
    end
    if it == 1
        out[:scf_converged] = false
        haskey(out, :scf_steps) ? out[:scf_steps] += 1 : out[:scf_steps] = 1
    end
end

function qe_parse_colin_magmoms(out, line, f)
    key = :colin_mag_moments
    out[key] = Float64[]
    line = readline(f)
    while !isempty(line)
        push!(out[key], parse.(Float64, split(line)[end]))
        line = readline(f)
    end
end

function qe_parse_scf_accuracy(out, line, f)
    key = :accuracy
    acc = parse(Float64, split(line)[5])
    if haskey(out, key)
        push!(out[key], acc)
    else
        out[key] = [acc]
    end
end

function qe_parse_total_magnetization(out, line, f)
    key = :total_magnetization
    mag = parse(Float64, split(line)[end-2])
    if haskey(out, key)
        push!(out[key], mag)
    else
        out[key] = [mag]
    end
end

function qe_parse_magnetization(out, line, f)
    if !haskey(out, :magnetization)
        out[:magnetization] = Vec3{Float64}[]
    end
    atom_number = parse(Int, split(line)[3])
    readline(f)
    if length(out[:magnetization]) < atom_number
        push!(out[:magnetization], parse(Vec3{Float64}, split(readline(f))[3:5]))
    else
        out[:magnetization][atom_number] = parse(Vec3{Float64}, split(readline(f))[3:5])
    end
end

function qe_parse_Hubbard(out, line, f)
    if !haskey(out, :Hubbard)
        out[:Hubbard] = [parse_Hubbard_block(f)]
    else
        push!(out[:Hubbard], parse_Hubbard_block(f))
    end
end

function qe_parse_timing(out, line, f)
    out[:timing] = TimingData[]
    curparent = ""
    while !occursin("PWSCF", line)
        isempty(line) && (line = readline(f); continue)
        sline = split(line)
        if line[end] == ':' # descent into call case
            curparent = String(sline[end][1:end-1])
        elseif length(sline) == 9 # normal call
            td = TimingData(String(sline[1]), qe_parse_time(sline[3]),
                            qe_parse_time(sline[5]), parse(Int, sline[8]), TimingData[])
            push!(out[:timing], td)
            if !isempty(curparent) # Child case
                if curparent[1] == '*'
                    if td.name[1] == 'c' || td.name[1] == 'r'
                        curparent = replace(curparent, '*' => td.name[1])
                        parent = getfirst(x -> x.name == curparent, out[:timing])
                        curparent = replace(curparent, td.name[1] => '*')
                    else
                        parent = getfirst(x -> occursin(curparent[2:end], x.name),
                                          out[:timing])
                    end
                else
                    parent = getfirst(x -> x.name == curparent, out[:timing])
                end
                push!(parent.children, td)
            end
        elseif sline[1] == "PWSCF" # Final PWSCF report
            push!(out[:timing],
                  TimingData("PWSCF", qe_parse_time(sline[3]), qe_parse_time(sline[5]), 1,
                             TimingData[]))
        end
        line = strip(readline(f))
    end
    # cleanup
    for td in out[:timing]
        id = findfirst(x -> x == ':', td.name)
        td.name = id !== nothing ? td.name[id+1:end] : td.name
    end
end

function qe_parse_starting_magnetization(out, line, f)
    readline(f)
    out[:starting_magnetization] = Dict{Symbol, Vec3}()
    line = readline(f)
    while !isempty(line)
        sline = split(line)
        atsym = Symbol(sline[1])
        mag = parse.(Float64, sline[2:end])
        out[:starting_magnetization][atsym] = length(mag) == 1 ?  Vec3(0.0, 0.0, mag[1]) : Vec3(mag...)
        line = readline(f)
    end
end

function qe_parse_starting_simplified_dftu(out, line, f)
    readline(f)
    out[:starting_simplified_dftu] = Dict{Symbol, DFTU}()
    line = readline(f)
    while !isempty(line)
        sline = split(line)
        atsym = Symbol(sline[1])
        L = parse(Int, sline[2])
        vals = parse.(Float64, sline[3:end])
        out[:starting_simplified_dftu][atsym] = DFTU{Float64}(l = L, U = vals[1], α=vals[2], J0=vals[3], β=vals[4])
        line = readline(f)
    end
end


const QE_PW_PARSE_FUNCTIONS = ["C/m^2" => qe_parse_polarization,
                               "lattice parameter" => qe_parse_lattice_parameter,
                               "number of Kohn-Sham states" => qe_parse_n_KS,
                               "crystal axes" => qe_parse_crystal_axes,
                               "reciprocal axes" => qe_parse_reciprocal_axes,
                               "atomic species   valence    mass" => qe_parse_atomic_species,
                               "number of atoms/cell" => qe_parse_nat,
                               "Crystallographic axes" => qe_parse_crystal_positions,
                               "PseudoPot" => qe_parse_pseudo,
                               "the Fermi energy is" => qe_parse_fermi,
                               "highest occupied" => qe_parse_highest_lowest,
                               "total energy  " => qe_parse_total_energy,
                               "SPIN UP" => (x, y, z) -> x[:colincalc] = true,
                               "cryst." => qe_parse_k_cryst, "cart." => qe_parse_k_cart,
                               "bands (ev)" => qe_parse_k_eigvals,
                               "End of self-consistent" => (x, y, z) -> haskey(x,
                                                                               :k_eigvals) &&
                                   empty!(x[:k_eigvals]),
                               "End of band structure" => (x, y, z) -> haskey(x,
                                                                              :k_eigvals) &&
                                   empty!(x[:k_eigvals]),
                               "CELL_PARAMETERS (" => qe_parse_cell_parameters,
                               "ATOMIC_POSITIONS (" => qe_parse_atomic_positions,
                               "Total force" => qe_parse_total_force,
                               "iteration #" => qe_parse_scf_iteration,
                               "Magnetic moment per site" => qe_parse_colin_magmoms,
                               "estimated scf accuracy" => qe_parse_scf_accuracy,
                               "total magnetization" => qe_parse_total_magnetization,
                               "convergence has been" => (x, y, z) -> x[:scf_converged] = true,
                               "Begin final coordinates" => (x, y, z) -> x[:converged] = true,
                               "atom number" => qe_parse_magnetization,
                               "--- enter write_ns ---" => qe_parse_Hubbard,
                               "init_run" => qe_parse_timing,
                               "Starting magnetic structure" => qe_parse_starting_magnetization,
                               "Simplified LDA+U calculation" => qe_parse_starting_simplified_dftu,
                               ]

"""
    qe_read_pw_output(filename::String; parse_funcs::Vector{Pair{String}}=Pair{String,<:Function}[])

Reads a pw quantum espresso calculation, returns a dictionary with all found data in the file.
The additional `parse_funcs` should be of the form:
`func(out_dict, line, f)` with `f` the file. 
"""
function qe_read_pw_output(filename::String;
                           parse_funcs::Vector{<:Pair{String}} = Pair{String}[])
    out = parse_file(filename, QE_PW_PARSE_FUNCTIONS; extra_parse_funcs = parse_funcs)
    if haskey(out, :in_alat) &&
       haskey(out, :in_cell) &&
       (haskey(out, :in_cart_positions) || haskey(out, :in_cryst_positions))
        cell_data = InputData(:cell_parameters, :alat, pop!(out, :in_cell))
        if haskey(out, :in_cryst_positions)
            atoms_data = InputData(:atomic_positions, :crystal,
                                   pop!(out, :in_cryst_positions))
        else
            atoms_data = InputData(:atomic_positions, :alat, pop!(out, :in_cart_positions))
        end
        pseudo_data = InputData(:atomic_species, :none, out[:pseudos])
        tmp_flags = Dict{Symbol,Any}(:ibrav => 0)
        tmp_flags[:A] = out[:in_alat]
        out[:initial_structure] = extract_structure!("initial", tmp_flags, cell_data,
                                                     out[:atsyms], atoms_data, pseudo_data)
        # Add starting mag and DFTU
        if haskey(out, :starting_magnetization)
            set_magnetization!(out[:initial_structure], pairs(out[:starting_magnetization])...; print=false)
        end
        if haskey(out, :starting_simplified_dftu)
            dftus = out[:starting_simplified_dftu]
            for (atsym, dftu) in dftus
                for a in out[:initial_structure][atsym]
                    a.dftu = dftu
                end
            end
        end
            
    end

    # Process final Structure
    if haskey(out, :pos_option) && haskey(out, :alat) && haskey(out, :cell_parameters)
        pseudo_data = InputData(:atomic_species, :none, out[:pseudos])
        tmp_flags = Dict{Symbol,Any}(:ibrav => 0)
        if haskey(out, :alat)
            tmp_flags[:A] = out[:alat] == :angstrom ? 1.0 :
                            (out[:alat] == :crystal ? out[:in_alat] :
                             conversions[:bohr2ang] * out[:alat])
        else
            tmp_flags[:A] = 1.0
        end
        cell_data = InputData(:cell_parameters, :alat, out[:cell_parameters])
        atoms_data = InputData(:atomic_positions, out[:pos_option], out[:atomic_positions])
        out[:final_structure] = extract_structure!("final", tmp_flags, cell_data,
                                                   out[:atsyms], atoms_data, pseudo_data)
        # Add starting mag and DFTU
        if haskey(out, :starting_magnetization)
            set_magnetization!(out[:initial_structure], pairs(out[:starting_magnetization])...; print=false)
        end
        if haskey(out, :starting_simplified_dftu)
            dftus = out[:starting_simplified_dftu]
            for (atsym, dftu) in dftus
                for a in out[:initial_structure][atsym]
                    a.dftu = dftu
                end
            end
        end
    end

    #process bands
    if haskey(out, :k_eigvals) &&
       !isempty(out[:k_eigvals]) &&
       haskey(out, :k_cart) &&
       haskey(out, :in_recip_cell)
        if !haskey(out, :k_cryst) && haskey(out, :in_recip_cell) && haskey(out, :k_cart)
            out[:k_cryst] = (v = (out[:in_recip_cell]^-1,) .* out[:k_cart].v,
                             w = out[:k_cart].w)
        end
        if get(out, :colincalc, false)
            out[:bands_up]   = [DFBand(out[:k_cart].v, out[:k_cryst].v, zeros(length(out[:k_cart].v))) for i in 1:length(out[:k_eigvals][1])]
            out[:bands_down] = [DFBand(out[:k_cart].v, out[:k_cryst].v, zeros(length(out[:k_cart].v))) for i in 1:length(out[:k_eigvals][1])]
        else
            out[:bands] = [DFBand(out[:k_cart].v, out[:k_cryst].v,
                                  zeros(length(out[:k_cart].v)))
                           for i in 1:length(out[:k_eigvals][1])]
        end
        for i in 1:length(out[:k_eigvals])
            for i1 in 1:length(out[:k_eigvals][i])
                if get(out, :colincalc, false)
                    if i <= length(out[:k_cart].v)
                        out[:bands_up][i1].eigvals[i] = out[:k_eigvals][i][i1]
                    else
                        out[:bands_down][i1].eigvals[i-length(out[:k_cart].v)] = out[:k_eigvals][i][i1]
                    end
                else
                    out[:bands][i1].eigvals[i] = out[:k_eigvals][i][i1]
                end
            end
        end
    end
    out[:converged] = get(out, :converged, false) ? true :
                      get(out, :scf_converged, false) && !haskey(out, :total_force)
    if haskey(out, :scf_iteration)
        out[:n_scf] = length(findall(i -> out[:scf_iteration][i+1] < out[:scf_iteration][i],
                                     1:length(out[:scf_iteration])-1))
    end
    for f in
        (:in_cart_positions, :in_alat, :in_cryst_positions, :alat, :pos_option, :pseudos,
         :cell_parameters, :in_recip_cell, :scf_converged, :atsyms, :nat, :k_eigvals,
         :k_cryst, :k_cart, :starting_simplified_dftu, :starting_magnetization)
        pop!(out, f, nothing)
    end
    return out
end

"""
    qe_read_kpdos(filename::String,column=1;fermi=0)

Reads the k_resolved partial density of states from a Quantum Espresso projwfc output file.
Only use this if the flag kresolveddos=true in the projwfc calculation file!! The returned matrix can be readily plotted using heatmap() from Plots.jl!
Optional calculation: column = 1 (column of the output, 1 = first column after ik and E)
fermi  = 0 (possible fermi offset of the read energy values)
Return:         Array{Float64,2}(length(k_points),length(energies)) ,
(ytickvals,yticks)
"""
function qe_read_kpdos(filename::String, column = 1; fermi = 0)
    read_tmp = readdlm(filename, Float64; comments = true)
    zmat     = zeros(typeof(read_tmp[1]), Int64(read_tmp[end, 1]), div(size(read_tmp)[1], Int64(read_tmp[end, 1])))
    for i1 in 1:size(zmat)[1]
        for i2 in 1:size(zmat)[2]
            zmat[i1, i2] = read_tmp[size(zmat)[2]*(i1-1)+i2, 2+column]
        end
    end

    yticks    = collect(Int(div(read_tmp[1, 2] - fermi, 1)):1:Int(div(read_tmp[end, 2] - fermi, 1)))
    ytickvals = [findfirst(x -> norm(yticks[1] + fermi - x) <= 0.1, read_tmp[:, 2])]
    for (i, tick) in enumerate(yticks[2:end])
        push!(ytickvals,
              findnext(x -> norm(tick + fermi - x) <= 0.1, read_tmp[:, 2], ytickvals[i]))
    end

    return unique(read_tmp[:, 2]), zmat', ytickvals, yticks
end

"""
    qe_read_pdos(filename::String)

Reads partial dos file.
"""
function qe_read_pdos(filename::String)
    read_tmp = readdlm(filename; skipstart = 1)
    energies = read_tmp[:, 1]
    values   = read_tmp[:, 2:end]

    return energies, values
end

function qe_read_projwfc_output(c::DFCalculation{QE}, args...; kwargs...)
    out = Dict{Symbol,Any}()
    pdos_files = searchdir(c, ".pdos_")
    if flag(c, :kresolveddos) == true
        out[:heatmaps]  = Vector{Matrix{Float64}}()
        out[:ytickvals] = Vector{Vector{Float64}}()
        out[:yticks]    = Vector{Vector{Float64}}()
        for f in pdos_files
            th, vals, ticks = qe_read_kpdos(f, args...)
            push!(out[:heatmaps], th)
            push!(out[:ytickvals], vals)
            push!(out[:yticks], ticks)
        end
    else
        out[:pdos] = NamedTuple{(:energies, :values),
                                Tuple{Vector{Float64},Matrix{Float64}}}[]
        for f in pdos_files
            energs, vals = qe_read_pdos(f, args...)
            push!(out[:pdos], (energies = energs, values = vals))
        end
    end
    out[:states], out[:bands] = qe_read_projwfc(outpath(c))
    return out
end

"""
    qe_read_projwfc(filename::String)

Reads the output file of a projwfc.x calculation.
Each kpoint will have as many energy dos values as there are bands in the scf/nscf calculation that
generated the density upon which the projwfc.x was called.
Returns:
    states: [(:atom_id, :wfc_id, :j, :l, :m),...] where each j==0 for a non spin polarized calculation.
    kpdos : kpoint => [(:e, :ψ, :ψ²), ...] where ψ is the coefficient vector in terms of the states.
"""
function qe_read_projwfc(filename::String)
    lines = readlines(filename) .|> strip

    i_prob_sizes = findfirst(x -> !isempty(x) && x[1:4] == "Prob", lines)
    istart = findfirst(x -> x == "Atomic states used for projection", lines) + 2

    natomwfc = 0
    nx       = 0
    nbnd     = 0
    nkstot   = 0
    npwx     = 0
    nkb      = 0
    for i in i_prob_sizes+1:istart-3
        l = lines[i]
        if isempty(l)
            break
        end
        sline = split(l)
        v = parse(Int, sline[3])
        if sline[1] == "natomwfc"
            natomwfc = v
        elseif sline[1] == "nx"
            nx = v
        elseif sline[1] == "nbnd"
            nbnd = v
        elseif sline[1] == "nkstot"
            nkstot = v
        elseif sline[1] == "npwx"
            npwx = v
        elseif sline[1] == "nkb"
            nkb = v
        end
    end

    state_tuple = NamedTuple{(:atom_id, :wfc_id, :l, :j, :m),
                             Tuple{Int,Int,Float64,Float64,Float64}}
    states = state_tuple[]
    for i in 1:natomwfc
        l = replace_multiple(lines[i+istart], "(" => " ", ")" => " ", "," => "", "=" => " ",
                             ":" => "", "#" => " ") |> split
        if length(l) == 11 #spinpolarized
            push!(states,
                  state_tuple((parse.(Int, (l[4], l[7]))..., parse(Float64, l[9]), 0.0,
                               parse(Float64, l[11]))))
        else #not spin polarized
            push!(states,
                  state_tuple((parse.(Int, (l[4], l[7]))...,
                               parse.(Float64, (l[9], l[11], l[13]))...)))
        end
    end
    ETuple = NamedTuple{(:e, :ψ, :ψ²),Tuple{Float64,Vector{Float64},Float64}}
    kdos = Pair{Vec3{Float64},Vector{ETuple}}[]
    while length(kdos) < nkstot
        istart   = findnext(x -> occursin("k = ", x), lines, istart + 1)
        k        = Vec3(parse.(Float64, split(lines[istart])[3:end]))
        etuples  = ETuple[]
        istop_ψ  = istart - 1
        istart_ψ = istart
        while length(etuples) < nbnd
            eline    = replace_multiple(lines[istop_ψ+2], "=" => "", "(" => " ", ")" => " ")
            e        = parse(Float64, split(eline)[end-1])
            coeffs   = zeros(length(states))
            istart_ψ = findnext(x -> !isempty(x) && x[1:3] == "===", lines, istop_ψ + 1) + 1
            istop_ψ  = findnext(x -> !isempty(x) && x[2:4] == "psi", lines, istart_ψ) - 1
            for i in istart_ψ:istop_ψ
                l = replace_multiple(lines[i], "psi =" => " ", "*[#" => " ", "]+" => " ",
                                     "]" => " ") |>
                    strip |>
                    split
                for k in 1:2:length(l)
                    coeffs[parse(Int, l[k+1])] = parse(Float64, l[k])
                end
            end
            ψ² = parse(Float64, split(lines[istop_ψ+1])[end])
            push!(etuples, (e = e, ψ = coeffs, ψ² = ψ²))
        end
        push!(kdos, k => etuples)
    end
    nkstot = length(kdos)
    nbnd   = length(last(kdos[1]))
    bands  = [DFBand(nkstot) for i in 1:nbnd]
    for b in bands
        b.extra[:ψ]  = Vector{Vector{Float64}}(undef, nkstot)
        b.extra[:ψ²] = Vector{Float64}(undef, nkstot)
    end

    for (i, (k, energies)) in enumerate(kdos)
        for (ie, etuple) in enumerate(energies)
            bands[ie].k_points_cryst[i] = k
            bands[ie].k_points_cart[i]  = zero(Vec3{Float64})
            bands[ie].eigvals[i]        = etuple.e
            bands[ie].extra[:ψ][i]      = etuple.ψ
            bands[ie].extra[:ψ²][i]     = etuple.ψ²
        end
    end
    return states, bands
end

function qe_parse_pert_at(out, line, f)
    sline = split(line)
    nat = parse(Int, sline[3])
    out[:pert_at] = []
    readline(f)
    for i in 1:nat
        sline = split(readline(f))
        push!(out[:pert_at],
              (name = Symbol(sline[2]),
               position = Point3(parse.(Float64, sline[end-3:end-1])...)))
    end
end

function qe_parse_Hubbard_U(out, line, f)
    out[:Hubbard_U] = []
    readline(f)
    readline(f)
    for i in 1:length(out[:pert_at])
        sline = split(readline(f))
        push!(out[:Hubbard_U],
              (orig_name = Symbol(sline[3]), new_name = Symbol(sline[6]),
               U = parse(Float64, sline[7])))
    end
end

const QE_HP_PARSE_FUNCS = ["will be perturbed" => qe_parse_pert_at,
                           "Hubbard U parameters:" => qe_parse_Hubbard_U]

function qe_read_hp_output(c::DFCalculation{QE}; parse_funcs = Pair{String,<:Function}[])
    out = parse_file(outpath(c), QE_HP_PARSE_FUNCS; extra_parse_funcs = parse_funcs)

    hubbard_file = joinpath(dir(c), """$(get(c, :prefix, "pwscf")).Hubbard_parameters.dat""")
    if ispath(hubbard_file)
        merge(out,
              parse_file(hubbard_file, QE_HP_PARSE_FUNCS; extra_parse_funcs = parse_funcs))
    end
    return out
end

function alat(flags, pop = false)
    if haskey(flags, :A)
        alat = pop ? pop!(flags, :A) : flags[:A]
        alat *= 1Ang
    elseif haskey(flags, :celldm_1)
        alat = pop ? pop!(flags, :celldm_1) : flags[:celldm_1]
        alat *= 1a₀
    elseif haskey(flags, :celldm)
        alat = pop ? pop!(flags, :celldm)[1] : flags[:celldm][1]
        alat *= 1a₀
    else
        error("Cell option 'alat' was found, but no matching flag was set. \n
               The 'alat' has to  be specified through 'A' or 'celldm(1)'.")
    end
    return alat
end

#TODO handle more fancy cells
function extract_cell!(flags, cell_block)
    if cell_block != nothing
        _alat = 1.0Ang
        if cell_block.option == :alat
            @assert pop!(flags, :ibrav) == 0 "Only ibrav = 0 allowed for now."
            _alat = alat(flags)

        elseif cell_block.option == :bohr
            _alat = 1u"a₀"
        end

        return (_alat .* cell_block.data)'
    end
end

function qe_DFTU(speciesid::Int, parsed_flags::SymAnyDict)
    U  = 0.0
    J0 = 0.0
    J  = [0.0]
    α  = 0.0
    β  = 0.0
    if haskey(parsed_flags, :Hubbard_U) && length(parsed_flags[:Hubbard_U]) >= speciesid
        U = parsed_flags[:Hubbard_U][speciesid]
    end
    if haskey(parsed_flags, :Hubbard_J0) && length(parsed_flags[:Hubbard_J0]) >= speciesid
        J0 = parsed_flags[:Hubbard_J0][speciesid]
    end
    if haskey(parsed_flags, :Hubbard_J) && length(parsed_flags[:Hubbard_J]) >= speciesid
        J = Float64.(parsed_flags[:Hubbard_J][:, speciesid])
    end
    if haskey(parsed_flags, :Hubbard_alpha) &&
       length(parsed_flags[:Hubbard_alpha]) >= speciesid
        α = parsed_flags[:Hubbard_alpha][speciesid]
    end
    if haskey(parsed_flags, :Hubbard_beta) &&
       length(parsed_flags[:Hubbard_beta]) >= speciesid
        β = parsed_flags[:Hubbard_beta][speciesid]
    end
    return DFTU{Float64}(; U = U, J0 = J0, α = α, β = β, J = sum(J) == 0 ? [0.0] : J)
end

degree2π(ang) = ang / 180 * π

function qe_magnetization(atid::Int, parsed_flags::SymAnyDict)
    θ = haskey(parsed_flags, :angle1) && length(parsed_flags[:angle1]) >= atid ?
        parsed_flags[:angle1][atid] : 0.0
    θ = degree2π(θ)
    ϕ = haskey(parsed_flags, :angle2) && length(parsed_flags[:angle2]) >= atid ?
        parsed_flags[:angle2][atid] : 0.0
    ϕ = degree2π(ϕ)

    start = haskey(parsed_flags, :starting_magnetization) &&
            length(parsed_flags[:starting_magnetization]) >= atid ?
            parsed_flags[:starting_magnetization][atid] : 0.0
    if start isa AbstractVector
        return Vec3{Float64}(start...)
    else
        return start * Vec3{Float64}(sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))
    end
end

function extract_atoms!(parsed_flags, atsyms, atom_block, pseudo_block,
                        cell::Mat3{LT}) where {LT<:Length}
    atoms = Atom{Float64,LT}[]

    option = atom_block.option
    if option == :crystal || option == :crystal_sg
        primv = cell
        cell  = Mat3(Matrix(1.0I, 3, 3))
    elseif option == :alat
        primv = alat(parsed_flags, true) * Mat3(Matrix(1.0I, 3, 3))
    elseif option == :bohr
        primv = 1a₀ .* Mat3(Matrix(1.0I, 3, 3))
    else
        primv = 1Ang .* Mat3(Matrix(1.0I, 3, 3))
    end
    for (at_sym, pos) in atom_block.data
        if haskey(pseudo_block.data, at_sym)
            pseudo = pseudo_block.data[at_sym]
        else
            elkey = getfirst(x -> x != at_sym && element(x) == element(at_sym),
                             keys(pseudo_block.data))
            pseudo = elkey !== nothing ? pseudo_block.data[elkey] :
                     (@warn "No Pseudo found for atom '$at_sym'.\nUsing Pseudo()."; Pseudo())
        end
        speciesid = findfirst(isequal(at_sym), atsyms)
        push!(atoms,
              Atom(; name = at_sym, element = element(at_sym), position_cart = primv * pos,
                   position_cryst = ustrip.(inv(cell) * pos), pseudo = pseudo,
                   magnetization = qe_magnetization(speciesid, parsed_flags),
                   dftu = qe_DFTU(speciesid, parsed_flags)))
    end

    return atoms
end

function extract_structure!(name, parsed_flags, cell_block, atsyms, atom_block,
                            pseudo_block)
    if atom_block == nothing
        return nothing
    end
    cell = extract_cell!(parsed_flags, cell_block)
    atoms = extract_atoms!(parsed_flags, atsyms, atom_block, pseudo_block, cell)
    return Structure(name, cell, atoms)
end

function separate(f, A::AbstractVector{T}) where {T}
    true_part = T[]
    false_part = T[]
    while length(A) > 0
        t = pop!(A)
        if f(t)
            push!(true_part, t)
        else
            push!(false_part, t)
        end
    end
    return reverse(true_part), reverse(false_part)
end

"""
    qe_read_calculation(filename, T=Float64; execs=[Exec(exec="pw.x")], run=true, structure_name="noname")

Reads a Quantum Espresso calculation file. The `QE_EXEC` inside execs gets used to find which flags are allowed in this calculation file, and convert the read values to the correct Types.
Returns a `DFCalculation{QE}` and the `Structure` that is found in the calculation.
"""
function qe_read_calculation(filename; execs = [Exec(; exec = "pw.x")], run = true,
                             structure_name = "noname")
    @assert ispath(filename) "$filename is not a valid path."
    t_lines = read(filename) |>
              String |>
              x -> split(x, "\n") .|> x -> cut_after(x, '!') .|> x -> cut_after(x, '#')
    lines = join(t_lines, "\n") |>
            # x -> replace(x, "," => "\n")  |>
            # x -> replace(x, "," => " ")   |>
            x -> replace(x, "'" => " ") |>
                 x -> split(x, "\n") .|>
                      strip |>
                      x -> filter(y -> !occursin("&", y), x) |>
                           x -> filter(y -> !(occursin("/", y) && length(y) == 1), x) |>
                                x -> filter(!isempty, x)

    exec = getfirst(x -> x.exec ∈ QE_EXECS, execs)

    flaglines, lines = separate(x -> occursin("=", x), lines)
    flaglines = strip_split.(flaglines, "=")
    easy_flaglines, difficult_flaglines = separate(x -> !occursin("(", x[1]), flaglines)
    parsed_flags = SymAnyDict()
    #easy flags
    for (f, v) in easy_flaglines
        sym = Symbol(f)
        typ = flagtype(QE, exec, sym)
        if eltype(typ) <: Bool
            v = replace(lowercase(v), "." => "")
        elseif eltype(typ) <: Number
            v = replace(v, "d" => "e")
        elseif typ === Nothing
            @warn "Flag $sym not found in allowed flags and will be ignored."
            continue
        end
        if typ <: AbstractArray
            typ = eltype(typ)
        end
        tval = typ != String ? parse.((typ,), split(v)) : v
        parsed_flags[sym] = length(tval) == 1 ? tval[1] : tval
    end

    used_lineids = Int[]
    findcard(s) = findfirst(l -> occursin(s, lowercase(l)), lines)
    i_species = findcard("atomic_species")
    i_cell = findcard("cell_parameters")
    i_positions = findcard("atomic_positions")
    if i_species !== nothing && i_cell !== nothing && i_positions !== nothing
        push!(used_lineids, i_species)
        nat  = parsed_flags[:nat]
        ntyp = parsed_flags[:ntyp]

        pseudos = InputData(:atomic_species, :none, Dict{Symbol,Pseudo}())
        pseudo_dir = string(pop!(parsed_flags, :pseudo_dir, "./"))
        atsyms = Symbol[]
        for k in 1:ntyp
            push!(used_lineids, i_species + k)
            sline = strip_split(lines[i_species+k])
            atsym = Symbol(sline[1])
            pseudos.data[atsym] = Pseudo(sline[end], pseudo_dir)
            push!(atsyms, atsym)
        end

        append!(used_lineids, [i_cell, i_cell + 1, i_cell + 2, i_cell + 3])
        cell_block = InputData(:cell_parameters, cardoption(lines[i_cell]),
                               Mat3([parse(Float64, split(lines[i_cell+k])[j])
                                     for k in 1:3, j in 1:3]))

        push!(used_lineids, i_positions)
        atom_block = InputData(:atomic_positions, cardoption(lines[i_positions]),
                               Tuple{Symbol,Point3{Float64}}[])
        for k in 1:nat
            push!(used_lineids, i_positions + k)
            sline = split(lines[i_positions+k])
            atsym = Symbol(sline[1])
            point = Point3(parse.(Float64, sline[2:4]))
            push!(atom_block.data, (atsym, point))
        end

        #the difficult flags, can only be present if atomic stuff is found
        for (f, v) in difficult_flaglines
            try
                _s = split(replace(replace(replace(f, "(" => " "), ")" => " "), "," => " "))

                sym = Symbol(_s[1])
                ids = parse.(Int, _s[2:end])
                typ = flagtype(QE, exec, sym)
                v = replace(v, "d" => "e")
                if typ === Nothing
                    @warn "Flag $f in file $filename not found in allowed flags for $(exec.exec)"
                    continue
                end
                parsedval = parse.((eltype(typ),), split(v))
                if !haskey(parsed_flags, sym)
                    if typ <: AbstractMatrix
                        parsed_flags[sym] = length(parsedval) == 1 ?
                                            zeros(eltype(typ), ntyp, 10) :
                                            fill(zeros(eltype(typ), length(parsedval)),
                                                 ntyp, 10) #arbitrary limit
                    elseif typ <: AbstractVector
                        parsed_flags[sym] = length(parsedval) == 1 ?
                                            zeros(eltype(typ), ntyp) :
                                            fill(zeros(eltype(typ), length(parsedval)),
                                                 ntyp)
                    elseif sym == :Hubbard_J
                        parsed_flags[sym] = zeros(eltype(typ), 3, nat)
                    elseif sym == :starting_ns_eigenvalue
                        #7 and 4 are the largest possible, and in the end if they are not filled
                        #they won't be written in the new input anyway
                        parsed_flags[sym] = zeros(eltype(typ), 7, 4, nat) 
                    end
                end
                parsed_flags[sym][ids...] = length(parsedval) == 1 ? parsedval[1] :
                                            parsedval
            catch e
                @warn "Parsing error of flag $f in file $filename." exception = e
            end
        end
        structure = extract_structure!(structure_name, parsed_flags, cell_block, atsyms,
                                       atom_block, pseudos)
        delete!.((parsed_flags,), [:ibrav, :nat, :ntyp, :A, :celldm_1, :celldm])
        delete!.((parsed_flags,),
                 [:Hubbard_U, :Hubbard_J0, :Hubbard_alpha, :Hubbard_beta, :Hubbard_J])
        delete!.((parsed_flags,), [:starting_magnetization, :angle1, :angle2, :nspin]) #hubbard and magnetization flags

    else
        structure = nothing
    end

    datablocks = InputData[]
    i = findcard("k_points")
    if i !== nothing
        append!(used_lineids, [i, i + 1])
        k_option = cardoption(lines[i])
        if k_option == :automatic
            s_line = split(lines[i+1])
            k_data = parse.(Int, s_line)
        else
            nks    = parse(Int, lines[i+1])
            k_data = Vector{NTuple{4,Float64}}(undef, nks)
            for k in 1:nks
                push!(used_lineids, i + 1 + k)
                k_data[k] = (parse.(Float64, split(lines[i+1+k]))...,)
            end
        end
        push!(datablocks, InputData(:k_points, k_option, k_data))
    end
    #Now we need to deal with lines that were not flaglines and not inside the usual cards.
    #Ultimately the goal would be to actually take the info read from the documentation to
    #understand the syntax of the cards, and read them accordingly. In the meantime we do
    #it like this.
    remlines = String[]
    for i in 1:length(lines)
        if i ∈ used_lineids
            continue
        end
        push!(remlines, lines[i])
    end
    !isempty(remlines) && push!(datablocks, InputData(:noname, :nooption, remlines))

    pop!.((parsed_flags,), [:prefix, :outdir], nothing)
    dir, file = splitdir(filename)
    return DFCalculation{QE}(name = splitext(file)[1], dir = dir, flags = parsed_flags, data = datablocks, execs = execs, run = run),
           structure
end

function qe_writeflag(f, flag, value)
    if isa(value, Vector)
        for i in 1:length(value)
            if !iszero(value[i])
                if length(value[i]) == 1
                    write(f, "  $flag($i) = $(value[i])\n")
                else
                    write(f, "  $flag($i) =")
                    for v in value[i]
                        write(f, " $v")
                    end
                    write(f, "\n")
                end
            end
        end
    elseif isa(value, Array)
        cids = CartesianIndices(value)
        for i in eachindex(value)
            if !iszero(value[i])
                write(f, "  $(flag)$(Tuple(cids[i])) = $(value[i])\n")
            end
        end
    elseif isa(value, AbstractString)
        write(f, "  $flag = '$value'\n")
    else
        write(f, "  $flag = $value\n")
    end
end

"""
    save(calculation::DFCalculation{QE}, structure, filename::String=inpath(calculation))

Writes a Quantum Espresso calculation file.
"""
function save(calculation::DFCalculation{QE}, structure,
              filename::String = inpath(calculation); relative_positions = true)
    if haskey(flags(calculation), :calculation)
        set_flags!(calculation,
                   :calculation => replace(calculation[:calculation], "_" => "-");
                   print = false)
    end
    open(filename, "w") do f
        if exec(calculation, "ph.x") !== nothing
            write(f, "--\n")
        end
        writeflag(flag_data) = qe_writeflag(f, flag_data[1], flag_data[2])
        write_dat(data) = write_data(f, data)

        controls = Dict{Symbol,Dict{Symbol,Any}}()

        for (flag, val) in calculation.flags
            block, variable = qe_block_variable(calculation, flag)
            if !haskey(controls, block.name)
                controls[block.name] = Dict{Symbol,Any}()
            end
            controls[block.name][flag] = val
        end

        #Here we try to figure out the correct order of the control blocks
        # first we find the order of the pw.x calculations, the rest should follow.
        blocks2file = []
        for name in [:control, :system, :electrons, :ions, :cell]
            push!(blocks2file, name => pop!(controls, name, nothing))
        end
        for name in keys(controls)
            push!(blocks2file, name => pop!(controls, name, nothing))
        end
        filter!(x -> x[2] !== nothing, blocks2file)
        for (name, flags) in blocks2file
            write(f, "&$name\n")
            if name == :system
                nat  = length(atoms(structure))
                ntyp = length(unique(atoms(structure)))
                # A     = 1.0
                ibrav = 0
                write(f, "  ibrav = $ibrav\n")
                # write(f,"  A = $A\n")
                write(f, "  nat = $nat\n")
                write(f, "  ntyp = $ntyp\n")
            end
            map(writeflag, [(flag, data) for (flag, data) in flags])
            write(f, "/\n\n")
        end
        if exec(calculation, "pw.x") !== nothing
            write_structure(f, calculation, structure;
                            relative_positions = relative_positions)
        end
        for dat in calculation.data
            if dat.name != :noname
                if dat.option != :none
                    write(f, "$(uppercase(String(dat.name))) ($(dat.option))\n")
                else
                    write(f, "$(uppercase(String(dat.name)))\n")
                end
            end
            if dat.name == :k_points && dat.option != :automatic
                write(f, "$(length(dat.data))\n")
                write_dat(dat.data)
            else
                write_dat(dat.data)
            end
            write(f, "\n")
        end
    end
    #TODO handle writing hubbard and magnetization better
    delete!.((calculation.flags,),
                    (:Hubbard_U, :Hubbard_J0, :Hubbard_J, :Hubbard_alpha, :Hubbard_beta,
                     :starting_magnetization, :angle1, :angle2, :pseudo_dir))
    return
end

function write_structure(f, calculation::DFCalculation{QE}, structure;
                         relative_positions = true)
    unique_at    = unique(atoms(structure))
    pseudo_lines = String[]
    atom_lines   = String[]
    for at in unique_at
        push!(pseudo_lines,
              "$(name(at)) $(element(at).atomic_weight)   $(pseudo(at).name)\n")
    end

    for at in atoms(structure)
        push!(atom_lines, position_string(QE, at; relative = relative_positions))
    end

    write(f, "ATOMIC_SPECIES\n")
    write.((f,), pseudo_lines)

    write(f, "\n")
    write(f, "CELL_PARAMETERS (angstrom)\n")
    write_cell(f, (ustrip.(uconvert.(Ang, cell(structure))))')
    write(f, "\n")

    if relative_positions
        write(f, "ATOMIC_POSITIONS (crystal) \n")
    else
        write(f, "ATOMIC_POSITIONS (angstrom) \n")
    end
    write.((f,), atom_lines)
    write(f, "\n")
    return 
end

function qe_generate_pw2wancalculation(calculation::DFCalculation{Wannier90},
                                       nscf_calculation::DFCalculation{QE}, runexecs)
    flags = Dict()
    flags[:prefix] = nscf_calculation[:prefix]
    flags[:seedname] = "$(name(calculation))"
    flags[:outdir] = nscf_calculation[:outdir]
    flags[:wan_mode] = "standalone"
    flags[:write_mmn] = true
    flags[:write_amn] = true
    if flag(calculation, :spin) !== nothing
        flags[:spin_component] = flag(calculation, :spin)
    end
    if flag(calculation, :spinors) !== nothing
        flags[:write_spn] = flag(calculation, :spinors)
    end
    if flag(calculation, :wannier_plot) !== nothing
        flags[:write_unk] = flag(calculation, :wannier_plot)
    end
    if any(flag(calculation, :berry_task) .== ("morb"))
        flags[:write_uHu] = true
    end
    pw2wanexec = Exec("pw2wannier90.x", runexecs[2].dir)
    run = get(calculation.flags, :preprocess, false) && calculation.run
    return DFCalculation{QE}(name = "pw2wan_$(flags[:seedname])", dir = dir(calculation), flags = flags,
                             data = InputData[], execs = [runexecs[1], pw2wanexec], run = run)
end
