"""
    read_abi_input(filename::String, T=Float64)

Returns an ABINIT input. We assume that jdtset is on a seperate line. 
If # DATASET is supplied as first line, it will make sure that amount of datasets are read no matter what ndtset and jdtset are.
ndtset and jdtset will be taken into account to decide which calculations will be marked as 'should run'.
"""
function read_abi_input(filename::String, T=Float64; run_command= "", pseudos=[""])
    datasets = read_abi_datasets(filename, T)
    structures = extract_structure(datasets...)
    # structures = [structures[1] ; filter(x -> x[1] != structures[1][1] && x[2] != structures[1][2], structures)]
    inputs = Array{AbinitInput,1}()
    jdtset = Int[]
    ndtset = 0
    println(datasets)
    for (i, data) in enumerate(datasets)
        file, ext= splitext(filename)
        if length(datasets) > 1
            newfile = file * "$i" * ext
        else
            newfile = file * ext
        end
        if haskey(data, :jdtset)
            jdtset = pop!(data,:jdtset)
        end
        run = i in jdtset

        push!(inputs, AbinitInput(newfile, datasets, [structures[i][1], structures[i][2],AbinitDataBlock(:pseudos, :pseudos, pseudos)], run_command, run))
    end
    return inputs
end

@inline function filter_comment(line)
    for c in ['#', '!']
        t = findfirst(line,c)
        if t == 0
            line = line
        elseif t == 1
            return ""
        else
            line = line[1:t - 1]
        end
    end
    return line
end

function expand_star_syntax(s::String)
    s = strip(s)
    if !contains(s,"*")
        return s
    end
    # Handle e.g `pawecutdg*`
    if typeof(parse(string(s[1]))) == Symbol && s[end - 1] == '*' 
        return s
    end
    flag, rest = split(s)
    spl = strip_split(rest, "*")

    # Handle "*2" case i.e. return "*2"
    if length(spl) == 1
        @assert s[1] != '*' "error in expanding start syntax"
        return string(s[1])
    end

    val = spl[2]
    out = "$flag $val"
    len = parse(string(spl[1]))
    for i = 1:len-1
        out *= " $val"
    end
    return out
end

@inline function eval_abinit_operator(s)
    spl = split(s)
    out = "$(spl[1])"
    for p in spl[2:end]
        out *= " $(eval(parse(string(p))))"
    end
    return out
end

"series like 10 ecut+ 5 are not supported"
function read_abi_datasets(filename::String, T=Float64)
    lines = readlines(filename) .|>
        strip                   .|>
        lowercase               .|>
        filter_comment           
    lines = filter(x -> length(x) > 1, lines) .|>
        expand_star_syntax                    .|>
        eval_abinit_operator     

    lines = split(join(lines," "))

    datasets = Array{OrderedDict{Symbol, Any},1}()
    dataset = OrderedDict{Symbol,Any}()
    flag = Symbol("start")
    for (l, line) in enumerate(lines)
        if contains(line, "?") || contains(line, ":")
            continue
        end
        #convert to the correct abinit units
        if line in ABI_UNIT_NAMES 
            dataset[flag] = convert_2abi(dataset[flag], line)
            continue
        end
        if typeof(parse(line)) == Symbol
            j = tryparse(Int, string(line[end]))
            if !isnull(j)
                dtset = get(j)+1
                while dtset > length(datasets)
                    push!(datasets, OrderedDict{Symbol, Any}())
                end

                dataset = datasets[dtset]
                flag = Symbol(line[1:end-1])
            else
                dtset = 1
                while dtset > length(datasets)
                    push!(datasets, OrderedDict{Symbol, Any}())
                end
                dataset = datasets[dtset]
                flag = Symbol(line)
            end
        else
            if flag == :ndtset
                continue
            end
            if haskey(dataset,flag)
                push!(dataset[flag], parse(line))
            else
                dataset[flag] = [parse(line)]
            end
        end
    end
    #make flags that are length one just be their value and that they are in the abivars
    for (d,data) in enumerate(datasets)
        for (flag, val) in data
            flag_type = get_abi_flag_type(flag)
            if flag_type == Void
                error("Flag $flag in dataset $d not found in abinit variables.")
            elseif flag_type != typeof(val)
                try
                    val = convert.(flag_type, val)
                catch
                    error("Flag type of flag $flag not correct, found: $(typeof(val)) expected:$flag_type.")
                end
            elseif flag == :spgroup || flag == :nobj
                error("Abinit spgroup builder is not supported. Structure must be given explicitly!")
            end
            if length(val) == 1
                data[flag] = val[1]
            end
        end
    end

    return datasets
end

function extract_structure(abi_datasets...)
    #possible different structures for different inputs
    out_structures = Array{Array{AbinitDataBlock, 1}, 1}()
    natom = 1
    znucl = 1
    ntypat =1
    typat = [1]
    for data in abi_datasets
        acell = pop!(data, :acell, [1.0,1.0,1.0])
        rprimd = zeros(3,3)
        if haskey(data, :angdeg) 
            if haskey(data, :rprim)
                error("rprim and angdeg cannot be used at the same time.")
            end
            angdeg = pop!(data, :angdeg)
            if !all(x -> x > 0, angdeg)
                error("Angles must be > 0 but got $angdeg")
            elseif sum(angdeg)
                error("The sum of angdeg must be lower that 360, angdeg $angdeg")
            end
            tol12 = 1e-12
            rprim = zeros(3,3)
            if abs(angdeg[1] -angdeg[2]) < tol12 && abs(angdeg[2] - angdeg[3]) < tol12 &&
            abs(angdeg[1]-90.) + abs(angdeg[2]-90.) + abs(angdeg[3] -90.) > tol12

            # Treat the case of equal angles (except all right angles):
            # generates trigonal symmetry wrt third axis
                cosang = cos(pi * angdeg[1]/180.0)
                a2 = 2.0 / 3.0 * (1.0 - cosang)
                aa = sqrt(a2)
                cc = sqrt(1.0-a2)
                rprim[1,1] = aa     
                rprim[1,2] = 0.0              
                rprim[1,3] = cc
                rprim[2,1] = -0.5 * aa
                rprim[2,2] = sqrt(3.0) * 0.5 * aa
                rprim[2,3] = cc
                rprim[3,1] = -0.5 * aa
                rprim[3,2] = -sqrt(3.0) * 0.5 * aa
                rprim[3,3] = cc
            else:
                # Treat all the other cases
                rprim[1,1] = 1.0
                rprim[2,1] = cos(pi * angdeg[3] / 180.)
                rprim[2,2] = sin(pi * angdeg[3] / 180.)
                rprim[3,1] = cos(pi * angdeg[2] / 180.)
                rprim[3,2] = (cos(pi * angdeg[1] / 180.0) - rprim[2,1] * rprim[3,1]) / rprim[2,2]
                rprim[3,3] = sqrt(1.0 - rprim[3,1]^2 - rprim[3,2]^2)
            end

            rprimd = [acell[i] * rprim[i,:] for i=1:3]
        else
            rprimd = [acell[i] * reshape(pop!(data,:rprim, eye(3)),3,3)[i,:] for i=1:3]
        end
        
        ntypat = pop!(data, :ntypat, ntypat) 
        natom = pop!(data, :natom, natom)
        znucl = pop!(data, :znucl, znucl)
        #some kind of npsp not supported
        typat = pop!(data, :typat, typat)
        typat = typat[1:ntypat]
        
        atom_positions = OrderedDict{Symbol,Array{Point3D,1}}()
        if haskey(data, :xred)
            atoms = reshape(data[:xred],(div(length(data[:xred]),3), 3))
            pop!(data, :xred)
            for i=1:length(typat)
                typ = typat[i]
                atom = Array(atoms[i,:])
                nz = length(filter( x-> x==znucl[typ], znucl))
                z_sym = nz == 1 ? get_element_sym(znucl[typ]) : Symbol(String(get_element_sym(znucl[typ])) * "$typ")
                println(rprimd) 
                if haskey(atom_positions, z_sym)
                    push!(atom_positions[z_sym], Point3D(atom' * rprimd))
                else
                    atom_positions[z_sym] = [Point3D(atom' * rprimd)]
                end
            end


        elseif haskey(data, :xcart)
            atoms = reshape(data[:xcart],(div(length(data[:xcart],3),3), 3))
            pop!(data, :xcart)
            for (typ, atom) in zip(typat, rows(atoms))
                z = znucl[typ]
                nz = length(filter( x-> x==z, znucl))
                z_sym = nz == 1 ? get_element_sym(z) : Symbol(String(get_element_sym(z)) * "$typ")
                
                if haskey(atom_positions, z_sym)
                    push!(atom_positions[z_sym], Point3D(atom))
                else
                    atom_positions[z_sym] = [Point3D(atom)]
                end
            end

        elseif haskey(data, :xangst)
            atoms = reshape(data[:xangst],(div(length(data[:xangst],3),3), 3))
            pop!(data, :xangst)
            for (typ, atom) in zip(typat, rows(atoms))
                z = znucl[typ]
                nz = length(filter( x-> x==z, znucl))
                z_sym = nz == 1 ? get_element_sym(z) : Symbol(String(get_element_sym(z)) * "$typ")
                
                if haskey(atom_positions, z_sym)
                    push!(atom_positions[z_sym], abi_conversions[:ang] * Point3D(atom))
                else
                    atom_positions[z_sym] = [abi_conversions[:ang] * Point3D(atom)]
                end
            end
  
        else
            try 
                rprimd = out_structures[1][1].data
                atoms = out_structures[1][2].data
            catch 
                error("xred|xcart|xangs must be given in input.")
            end
        end
        push!(out_structures, [AbinitDataBlock(:cell_parameters, :bohr, rprimd), AbinitDataBlock(:atomic_positions, :cart, atoms)])
    end
    return out_structures
end


function write_abi_input(input::AbinitInput, filename::String=input.filename)
    flags = input.flags
    open(filename, "w") do f
        write_flag(flag_data) = write_flag_line(f, flag_data[1], flag_data[2], "")
        write(f, input.structure[:abi_string])
        write(f, "\n")
        write_flag.(collect(flags))
    end
end

#question: Does it matter that we might be writing redundant data such as spinat etc?
"""
    write_abi_datasets(inputs::Array{AbinitInput,1}, directory)

Takes all the inputs, sees which ones have the same structure and constructs input files for each seperate sturcture.
The filename of the first input of a certain structure is used as file for all the datasets.
"""
function write_abi_datasets(inputs::Array{AbinitInput,1}, directory)
    input_groups = Array{Array{AbinitInput,1},1}([[pop!(inputs)]])
    while length(inputs) != 0
        input = pop!(inputs)
        for group in input_groups
            if group[1].structure == input.structure
                push!(group, input)
                break
            end
        end
    end
    
    filenames    = String[]
    run_commands = String[]
    pseudos      = Array{Array{String,1},1}()

    for group in input_groups
        run_indices = Int[]
        for (i, _input) in enumerate(reverse(group))
            if _input.run 
                push!(run_indices, i)
            end
        end
        if length(run_indices) == length(group)
            jdtset = ""
            ndtset = length(group)
        else
            jdtset = join(["$c" for c in run_indices], " ")
            ndtset = length(run_indices)
        end
        file, ext = splitext(group[end].filename)  
        push!(filenames, file[1:end-1] * ext)
        push!(run_commands, group[end].run_command)
        push!(pseudos, get_data(group[end], :pseudos))
        open(directory * file[1:end-1] * ext, "w") do f
            write(f, "# DATASETS $(length(group))\n")
            write(f, "ndtset $ndtset\n")
            if jdtset != "" 
                write(f, "jdtset $jdtset\n")
            end
            write(f, group[1].structure[:abi_string])
            write(f, "\n")
            for (t, dt) in enumerate(reverse(group))
                write_flag(flag_data) = write_flag_line(f, flag_data[1], flag_data[2], "", t)
                write(f, "#=============== BEGIN DATASET $t ===============#\n")
                write_flag.(collect(dt.flags))
                write(f, "#================ END DATASET $t ================#\n")
            end
        end
    end
    return zip(filenames, pseudos, run_commands)
end

#very stupid
function read_abi_output(filename::String, T=Float64)
    if contains(filename, "FATBANDS")
        return read_abi_fatbands(filename, T)
    elseif contains(filename, "EBANDS.agr")
        return read_abi_ebands(filename, T)
    elseif contains(filename, "_EIG")
        return read_abi_eig(filename, T)
    else
        error("Please supply a file with FATBANDS, EBANDS.agr or _EIG in the filename.")
    end
end

"Reads an abinit FATBANDS output file and returns the found DFBands, with pdos values in the data field. K-points are not given (they aren't present in the output file)."
function read_abi_fatbands(filename::String, T=Float64)
    bands = DFBand[]
    open(filename, "r") do f
        while !eof(f)
            line = readline(f)
            if contains(line, "BAND number")
                extra   = OrderedDict{Symbol,Any}(:pdos => T[])
                eigvals = Array{T,1}()
                line    = readline(f)
                while line != "" && line != "&"
                    eigval, pdos = parse.(T, split(line)[2:end])
                    push!(eigvals, eigval)
                    push!(extra[:pdos], pdos)
                    line = readline(f)
                end
                push!(bands, DFBand{T}([T[0.0, 0.0, 0.0] for i=1:length(eigvals)], [T[0.0, 0.0, 0.0] for i=1:length(eigvals)], eigvals, extra))
            end
        end
    end
    return bands
end

"Reads an abinit EBANDS.agr output file and returns the found DFBands. K-points only given in crystalline coordinates."
function read_abi_ebands(filename::String, T=Float64)
    bands = DFBand[]
    open(filename, "r") do f
        k_points_cryst = Array{Array{T,1},1}()
        while !eof(f)
            line = readline(f)
            if contains(line,"List of k-points")
                line = readline(f)
                while line[1] != '@'
                    k_point = parse.(T, replace.(strip.(strip.(strip.(split(line)[4:end], '['), ']'), ','), 'E', 'e'))
                    push!(k_points_cryst, k_point)
                    line = readline(f)
                end
                
            elseif contains(line,"target")
                readline(f)
                eigvals = T[]
                line = readline(f) 
                while line[1] != '&'
                    push!(eigvals, parse(T, split(line)[2]))
                    line = readline(f)
                end
                push!(bands, DFBand([T[0.0, 0.0, 0.0] for i=1:length(eigvals)], k_points_cryst, eigvals))
            end
        end
    end
    return bands
end

"Reads and abinit _EIG output file and returns the found DFBands. K-points are only given in crystalline coordinates."
function read_abi_eig(filename::String, T=Float64)
    bands = DFBand[]
    k_points_array = Array{Array{T,1},1}()
    eigvals_array  = Array{Array{T,1},1}()
    open(filename, "r") do f
        while !eof(f)
            line = readline(f)
            if contains(line, "kpt#")
                push!(k_points_array, parse.(T, split(line)[7:9]))
                read_abi_eig_block!(f, k_points_array, eigvals_array, T)
            end
        end
    end
    for eigvals in eigvals_array
        push!(bands, DFBand([T[0.0, 0.0, 0.0] for i=1:length(eigvals)], k_points_array, eigvals))
    end
    return bands
end 

function read_abi_eig_block!(f, k_points_array, eigvals_array, T=Float64)
    first_run = isempty(eigvals_array)
    if eof(f)
        return
    end
    i = 1
    while !eof(f)
        line = readline(f)
        if contains(line,"Eigenvalues")
            continue
        elseif contains(line, "kpt#")
            push!(k_points_array, parse.(T, split(line)[7:9]))
            read_abi_eig_block!(f, k_points_array, eigvals_array, T)
        else
            for val in parse.(T, split(line))
                if first_run 
                    push!(eigvals_array, [val])
                else
                    push!(eigvals_array[i], val)
                end
                i += 1
            end
        end
    end
end


