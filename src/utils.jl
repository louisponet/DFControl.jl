"""
Searches a directory for all files containing the key.

Input: path::String,
       key::String
"""
search_dir(path::String,key) = filter(x->contains(x,key), readdir(path))

"""
Parse an array of strings into an array of a type.

Input:  T::Type,
        array::Array{String,1}
Return: Array{T,1}
"""
parse_string_array(T::Type,array) = map(x->(v = tryparse(T,x); isnull(v) ? 0.0 : get(v)),array)

"""
Parse a line for occurrences of type T.

Input:  T::Type,
        line::String
Return: Array{T,1}
"""
parse_line(T::Type,line::String) = parse_string_array(T,split(line))

"""
Mutatatively applies the fermi level to all eigvals in the band. If fermi is a quantum espresso scf output file it will try to find it in there.

Input:  band::Band,
        fermi::Union{String,AbstractFloat}
"""
function apply_fermi_level!(band::Band,fermi::Union{String,AbstractFloat})
  if typeof(fermi)==String
    fermi = read_fermi_from_qe_file(fermi)
  end
  for i=1:size(band.eigvals)[1]
    band.eigvals[i] -= fermi
  end
end

#@Cleanup Test if this is still necessary with the .f() syntax!!!
"""
Same as above but for an array of bands. Is this even necessary?
"""
function apply_fermi_level!{T<:Band}(bands::Array{T},fermi)
  for band in bands
    apply_fermi_level!(band,fermi)
  end
end

"""
Same as above but not mutatatively.
"""
function apply_fermi_level(band::Band,fermi)
  T = typeof(band.eigvals[1])
  if typeof(fermi)==String
    fermi = read_fermi_from_qe_file(fermi)
  end
  out=deepcopy(band)
  for i1=1:size(band.eigvals)[1]
    out.eigvals[i1] = band.eigvals[i1]-T(fermi)
  end
  return out
end
function apply_fermi_level{T<:Band}(bands::Array{T},fermi)
  out = similar(bands)
  for (i,band) in enumerate(bands)
    out[i] = apply_fermi_level(band,fermi)
  end
  return out
end

"""
Makes sure that a directory string ends with "/".

Input:  directory::String
Return: String
"""
function form_directory(directory::String)
  if directory[end]!='/'
    return directory*"/"
  else
    return directory
  end
end

"""
    gen_k_grid(na, nb, nc, input, T=Float32)

Returns an array of k-grid points that are equally spaced, input can be either `:wan` or `:nscf`, the returned grids are appropriate as inputs for wannier90 or an nscf calculation respectively.
"""
function gen_k_grid(na, nb, nc, input, T=Float32)
  if input == :wan || typeof(input) == WannierInput
    return [T[a,b,c] for a in collect(linspace(0,1,na+1))[1:end-1],b in collect(linspace(0,1,nb+1))[1:end-1],c in collect(linspace(0,1,nc+1))[1:end-1]]
  elseif input == :nscf || typeof(input) == QEInput
    return [T[a,b,c,1/(na*nb*nc)] for a in collect(linspace(0,1,na+1))[1:end-1],b in collect(linspace(0,1,nb+1))[1:end-1],c in collect(linspace(0,1,nc+1))[1:end-1]] 
  end
end

strip_split(line,args...) = strip.(split(line,args...))

#Incomplete for now only QE flags are returned
"""
    print_qe_flags(namelist_symbol::Symbol)

Prints the possible Quantum Espresso input flags and their type for a given input namelist.
"""
function print_qe_flags(namelist_symbol::Symbol)
  for block in QEControlFlags
    if block.name == namelist_symbol
      display(block)
    end
  end
end

"""
    print_qe_namelists()

Prints all the possible Quantum Espresso input namelists.
"""
function print_qe_namelists()
  for block in QEControlFlags
    println(block.name)
  end
end
