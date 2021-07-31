#
# function expand_star_syntax(s::String)
#     s = strip(s)
#     if !occursin("*", s)
#         return s
#     end
#     # Handle e.g `pawecutdg*`
#     if typeof(Meta.parse(string(s[1]))) == Symbol && s[end - 1] == '*'
#         return s
#     end
#     flag, rest = split(s)
#     spl = strip_split(rest, "*")
#
#     # Handle "*2" case i.e. return "*2"
#     if length(spl) == 1
#         @assert s[1] != '*' "error in expanding start syntax"
#         return string(s[1])
#     end
#
#     val = spl[2]
#     out = "$flag $val"
#     len = Meta.parse(string(spl[1]))
#     for i = 1:len-1
#         out *= " $val"
#     end
#     return out
# end
#
# @inline function eval_abinit_operator(s)
#     spl = split(s)
#     out = "$(spl[1])"
#     for p in spl[2:end]
#         out *= " $(eval(Meta.parse(string(p))))"
#     end
#     return out
# end
#
# "series like 10 ecut+ 5 are not supported"
# function read_abi_datasets(filename::String, T=Float64)
#     lines = readlines(filename) .|>
#         strip                   .|>
#         lowercase               .|>
#         filter_comment
#     lines = filter(x -> length(x) > 1, lines) .|>
#         expand_star_syntax                    .|>
#         eval_abinit_operator
#
#     lines = split(join(lines," "))
#
#     datasets = Vector{Dict{Symbol, Any}}()
#     dataset = Dict{Symbol,Any}()
#     flag = Symbol("start")
#     for (l, line) in enumerate(lines)
#         if occursin("?", line) || occursin(":", line)
#             continue
#         end
#         #convert to the correct abinit units
#         if line in ABI_UNIT_NAMES
#             dataset[flag] = convert_2abi(dataset[flag], line)
#             continue
#         end
#         if typeof(Meta.parse(line)) == Symbol
#             j = tryparse(Int, string(line[end]))
#             if !isnull(j)
#                 dtset = get(j)+1
#                 while dtset > length(datasets)
#                     push!(datasets, Dict{Symbol, Any}())
#                 end
#
#                 dataset = datasets[dtset]
#                 flag = Symbol(line[1:end-1])
#             else
#                 dtset = 1
#                 while dtset > length(datasets)
#                     push!(datasets, Dict{Symbol, Any}())
#                 end
#                 dataset = datasets[dtset]
#                 flag = Symbol(line)
#             end
#         else
#             if flag == :ndtset
#                 continue
#             end
#             if haskey(dataset,flag)
#                 push!(dataset[flag], Meta.parse(line))
#             else
#                 dataset[flag] = [Meta.parse(line)]
#             end
#         end
#     end
#     #make flags that are length one just be their value and that they are in the abivars
#     for (d,data) in enumerate(datasets)
#         for (flag, val) in data
#             flag_type = abi_flag_type(flag)
#             if flag_type == Nothing
#                 error("Flag $flag in dataset $d not found in abinit variables.")
#             elseif flag_type != typeof(val)
#                 try
#                     val = convert.(flag_type, val)
#                 catch
#                     error("Flag type of flag $flag not correct, found: $(typeof(val)) expected:$flag_type.")
#                 end
#             elseif flag == :spgroup || flag == :nobj
#                 error("Abinit spgroup builder is not supported. Structure must be given explicitly!")
#             end
#             if length(val) == 1
#                 data[flag] = val[1]
#             end
#         end
#     end
#
#     return datasets
# end
#
# #Everything is in angstrom, bohr lengths begone foul beasts!
# function extract_structures!(abi_datasets...; structure_name = "NoName")
#     #possible different structures for different calculations
#     out_structures = Vector{Vector{AbinitDataBlock}}()
#     out_structures = Vector{Structure}()
#     natom = 1
#     znucl = 1
#     ntypat =1
#     typat = [1]
#     for data in abi_datasets
#
#         acell = pop!(data, :acell, [1.0, 1.0, 1.0])
#         rprimd = zeros(3, 3)
#         if haskey(data, :angdeg)
#             if haskey(data, :rprim)
#                 error("rprim and angdeg cannot be used at the same time.")
#             end
#             angdeg = pop!(data, :angdeg)
#             if !all(x -> x > 0, angdeg)
#                 error("Angles must be > 0 but got $angdeg")
#             elseif sum(angdeg)
#                 error("The sum of angdeg must be lower that 360, angdeg $angdeg")
#             end
#             tol12 = 1e-12
#             rprim = zeros(3,3)
#             if abs(angdeg[1] -angdeg[2]) < tol12 && abs(angdeg[2] - angdeg[3]) < tol12 &&
#             abs(angdeg[1]-90.) + abs(angdeg[2]-90.) + abs(angdeg[3] -90.) > tol12
#
#             # Treat the case of equal angles (except all right angles):
#             # generates trigonal symmetry wrt third axis
#                 cosang = cos(pi * angdeg[1]/180.0)
#                 a2 = 2.0 / 3.0 * (1.0 - cosang)
#                 aa = sqrt(a2)
#                 cc = sqrt(1.0-a2)
#                 rprim[1,1] = aa
#                 rprim[1,2] = 0.0
#                 rprim[1,3] = cc
#                 rprim[2,1] = -0.5 * aa
#                 rprim[2,2] = sqrt(3.0) * 0.5 * aa
#                 rprim[2,3] = cc
#                 rprim[3,1] = -0.5 * aa
#                 rprim[3,2] = -sqrt(3.0) * 0.5 * aa
#                 rprim[3,3] = cc
#             else
#                 # Treat all the other cases
#                 rprim[1,1] = 1.0
#                 rprim[2,1] = cos(pi * angdeg[3] / 180.)
#                 rprim[2,2] = sin(pi * angdeg[3] / 180.)
#                 rprim[3,1] = cos(pi * angdeg[2] / 180.)
#                 rprim[3,2] = (cos(pi * angdeg[1] / 180.0) - rprim[2,1] * rprim[3,1]) / rprim[2,2]
#                 rprim[3,3] = sqrt(1.0 - rprim[3,1]^2 - rprim[3,2]^2)
#             end
#             for i = 1:3
#                 rprimd[i,:] *= acell[i]
#             end
#         else
#             rprimd = reshape(pop!(data, :rprim, eye(3)), (3,3))
#             for i=1:3
#                 rprimd[i,:] *= acell[i]
#             end
#         end
#         rprimd *= conversions[:bohr2ang]
#         rprimd = Mat3(rprimd)
#         ntypat = pop!(data, :ntypat, ntypat)
#         natom = pop!(data, :natom, natom)
#         znucl = pop!(data, :znucl, znucl)
#         #some kind of npsp not supported
#         typat = pop!(data, :typat, typat)
#         typat = typat[1:natom]
#
#         s_atoms = Atom{eltype(rprimd)}[]
#         if haskey(data, :xred)
#             xred = pop!(data, :xred)
#             atoms = reshape(xred, (div(length(xred), 3), 3))
#
#             for i=1:length(typat)
#                 typ = typat[i]
#                 atom = Array(atoms[i,:])
#                 z = znucl[typ]
#                 nz = length(filter( x-> x==z, znucl))
#                 _element = element(z)
#                 z_sym = nz == 1 ? _element.symbol : Symbol(String(_element.symbol) * "$typ")
#                 push!(s_atoms, Atom(z_sym, _element, rprimd' * Point3(atom)))
#             end
#
#
#         elseif haskey(data, :xcart)
#             xcart = pop!(data, :xcart)
#             atoms = reshape(xcart, (div(length(xcart), 3), 3))
#             for (typ, atom) in zip(typat, rows(atoms))
#                 z = znucl[typ]
#                 nz = length(filter( x-> x==z, znucl))
#                 _element = element(z)
#                 z_sym = nz == 1 ? _element.symbol : Symbol(String(_element.symbol) * "$typ")
#                 push!(s_atoms, Atom(z_sym, _element, conversions[:bohr2ang] * Point3(atom)))
#             end
#
#         elseif haskey(data, :xangst)
#             xangst = pop!(data, :xangst)
#             atoms = reshape(xangst, (div(length(xangst), 3), 3))
#             for (typ, atom) in zip(typat, rows(atoms))
#                 z = znucl[typ]
#                 nz = length(filter( x-> x==z, znucl))
#                 element = element(z)
#                 z_sym = nz == 1 ? _element.symbol : Symbol(String(_element.symbol) * "$typ")
#
#                 push!(s_atoms, Atom(z_sym, _element, Point3(atom)))
#             end
#
#         else
#             try
#                 rprimd = out_structures[1].cell
#                 atoms = out_structures[1].atoms
#             catch
#                 error("xred|xcart|xangs must be given in calculation.")
#             end
#         end
#         # up till here everything was in bohr, not anymore.
#         if !isempty(s_atoms)
#             push!(out_structures, Structure(structure_name, rprimd, s_atoms))
#         end
#     end
#     return out_structures
# end
#
# """
#     read_abi_calculation(filename::String, T=Float64)
#
# Returns an ABINIT calculation. We assume that jdtset is on a seperate line.
# If # DATASET is supplied as first line, it will make sure that amount of datasets are read no matter what ndtset and jdtset are.
# `ndtset` and `jdtset` will be taken into account to decide which calculations will be marked as 'should run'.
# Either all structures of the calculation are the same or all are different, otherwise there will be an error.
# """
# function read_abi_calculation(filename::String, T=Float64; runcommand= "", pseudos=[""], structure_name = "NoName")
#     datasets = read_abi_datasets(filename, T)
#     structures = extract_structures!(datasets..., structure_name = structure_name)
#     calculations = Array{AbinitInput,1}()
#     jdtset = Int[]
#     ndtset = 0
#     @assert length(datasets) == length(structures) || length(structures) == 1 "All structures of the calculation have to be the same or all different."
#     for (i, data) in enumerate(datasets)
#         structure = length(structures) == 1 ? structures[1] : structures[i]
#         file, ext= splitext(splitdir(filename)[2])
#         if length(datasets) > 1
#             newfile = file * "$i" * ext
#         else
#             newfile = file * ext
#         end
#         if haskey(data, :jdtset)
#             jdtset = pop!(data, :jdtset)
#         end
#         run = i in jdtset
#         if isempty(data)
#             continue
#         end
#         push!(calculations, AbinitInput(newfile, structure, data, [AbinitDataBlock(:pseudos, :pseudos, pseudos)], runcommand, run))
#     end
#     return calculations
# end
#
# function write_abi_structure(f, structure)
#     write(f, "acell 1.0 1.0 1.0\n")
#     write(f, "xangst\n")
#     for at in structure.atoms
#         write(f, position_cart(at))
#         write(f, "\n")
#     end
#     write(f, "rprim\n")
#     write_cell(f, structure.cell)
#
#     unique = unique(structure.atoms)
#     write(f, "nat $(length(structure.atoms))\n")
#     write(f, "ntypat $(length(unique))\n")
#     write(f, "typat\n")
#     i = 1
#     for at in structure.atoms
#         write(f, "$(findfirst(x->name(x) == name(at), unique)) ")
#         if i == 3
#             i=0
#             write(f,"\n")
#         end
#         i+=1
#     end
#     write(f, "\n")
#     write(f, "znucl ")
#     for at in unique
#         write(f, "$(element(at).Z) ")
#     end
#     write(f, "\n")
# end
#
# #question: Does it matter that we might be writing redundant data such as spinat etc?
# """
#     write_abi_datasets(calculations::Array{AbinitInput,1}, directory)
#
# Takes all the calculations, sees which ones have the same structure and constructs calculation files for each seperate sturcture.
# The filename of the first calculation of a certain structure is used as file for all the datasets.
# """
# function write_abi_datasets(calculations::Vector{DFCalculation{Abinit}}, directory)
#     calculation_groups = Vector{Vector{DFCalculation{Abinit}}}([[calculations[end]]])
#     for calculation in reverse(calculations)[2:end]
#         for group in calculation_groups
#             if group[1].structure == calculation.structure
#                 push!(group, calculation)
#                 break
#             end
#         end
#     end
#
#     filenames    = String[]
#     run_commands = String[]
#     pseudos      = Vector{Vector{String}}()
#
#     for group in calculation_groups
#         run_indices = Int[]
#         for (i, _calculation) in enumerate(reverse(group))
#             if _calculation.run
#                 push!(run_indices, i)
#             end
#         end
#         if length(run_indices) == length(group)
#             jdtset = ""
#             ndtset = length(group)
#         else
#             jdtset = join(["$c" for c in run_indices], " ")
#             ndtset = length(run_indices)
#         end
#         file, ext = splitext(group[end].filename)
#         push!(filenames, file[1:end-1] * ext)
#         push!(run_commands, group[end].runcommand)
#         push!(pseudos, data(group[end], :pseudos))
#         open(directory * file[1:end-1] * ext, "w") do f
#             write(f, "# DATASETS $(length(group))\n")
#             write(f, "ndtset $ndtset\n")
#             if jdtset != ""
#                 write(f, "jdtset $jdtset\n")
#             end
#             write_abi_structure(f, group[1].structure)
#             write(f, "\n")
#             for (t, dt) in enumerate(reverse(group))
#                 write_flag(flag_data) = write_flag_line(f, flag_data[1], flag_data[2], "", t)
#                 write(f, "#=============== BEGIN DATASET $t ===============#\n")
#                 write_flag.(collect(dt.flags))
#                 write(f, "#================ END DATASET $t ================#\n")
#             end
#         end
#     end
#     return zip(filenames, pseudos, run_commands)
# end
#
# #very stupid
function read_abi_output(filename::String, T = Float64)
    if occursin("FATBANDS", filename)
        return read_abi_fatbands(filename, T)
    elseif occursin("EBANDS.agr", filename)
        return read_abi_ebands(filename, T)
    elseif occursin("_EIG", filename)
        return read_abi_eig(filename, T)
    else
        error("Please supply a file with FATBANDS, EBANDS.agr or _EIG in the filename.")
    end
end
#
# "Reads an abinit FATBANDS output file and returns the found DFBands, with pdos values in the data field. K-points are not given (they aren't present in the output file)."
# function read_abi_fatbands(filename::String, T=Float64)
#     bands = DFBand[]
#     open(filename, "r") do f
#         while !eof(f)
#             line = readline(f)
#             if occursin("BAND number", line)
#                 extra   = Dict{Symbol,Any}(:pdos => T[])
#                 eigvals = Vector{T}()
#                 line    = readline(f)
#                 while line != "" && line != "&"
#                     eigval, pdos = Meta.parse.(T, split(line)[2:end])
#                     push!(eigvals, eigval)
#                     push!(extra[:pdos], pdos)
#                     line = readline(f)
#                 end
#                 push!(bands, DFBand{T}([T[0.0, 0.0, 0.0] for i=1:length(eigvals)], [T[0.0, 0.0, 0.0] for i=1:length(eigvals)], eigvals, extra))
#             end
#         end
#     end
#     return bands
# end
#
# "Reads an abinit EBANDS.agr output file and returns the found DFBands. K-points only given in crystalline coordinates."
function read_abi_ebands(filename::String, T = Float64)
    bands = DFBand[]
    open(filename, "r") do f
        k_points_cryst = Vector{Vector{T}}()
        while !eof(f)
            line = readline(f)
            if occursin("List of k-points", line)
                line = readline(f)
                while line[1] != '@'
                    k_point = parse.(T,
                                     replace.(strip.(strip.(strip.(split(line)[4:end], '['),
                                                            ']'), ','), 'E', 'e'))#ohno the replace here!
                    push!(k_points_cryst, k_point)
                    line = readline(f)
                end

            elseif occursin("target", line)
                readline(f)
                eigvals = T[]
                line = readline(f)
                while line[1] != '&'
                    push!(eigvals, parse(T, split(line)[2]))
                    line = readline(f)
                end
                push!(bands,
                      DFBand([T[0.0, 0.0, 0.0] for i in 1:length(eigvals)], k_points_cryst,
                             eigvals))
            end
        end
    end
    return bands
end
#
# "Reads and abinit _EIG output file and returns the found DFBands. K-points are only given in crystalline coordinates."
# function read_abi_eig(filename::String, T=Float64)
#     bands = DFBand[]
#     k_points_array = Vector{Vector{T}}()
#     eigvals_array  = Vector{Vector{T}}()
#     open(filename, "r") do f
#         while !eof(f)
#             line = readline(f)
#             if occursin("kpt#", line)
#                 push!(k_points_array, Meta.parse.(T, split(line)[7:9]))
#                 read_abi_eig_block!(f, k_points_array, eigvals_array, T)
#             end
#         end
#     end
#     for eigvals in eigvals_array
#         push!(bands, DFBand([T[0.0, 0.0, 0.0] for i=1:length(eigvals)], k_points_array, eigvals))
#     end
#     return bands
# end
#
# function read_abi_eig_block!(f, k_points_array, eigvals_array, T=Float64)
#     first_run = isempty(eigvals_array)
#     if eof(f)
#         return
#     end
#     i = 1
#     while !eof(f)
#         line = readline(f)
#         if occursin("Eigenvalues", line)
#             continue
#         elseif occursin("kpt#", line)
#             push!(k_points_array, Meta.parse.(T, split(line)[7:9]))
#             read_abi_eig_block!(f, k_points_array, eigvals_array, T)
#         else
#             for val in Meta.parse.(T, split(line))
#                 if first_run
#                     push!(eigvals_array, [val])
#                 else
#                     push!(eigvals_array[i], val)
#                 end
#                 i += 1
#             end
#         end
#     end
# end
