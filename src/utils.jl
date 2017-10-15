"""
Searches a directory for all files containing the key.

Input: path::String,
       key::String
"""
search_dir(path::String,key::String) = filter(x->contains(x,key), readdir(path))

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
  else
    for i=1:size(band.eigvals)[1]
      band.eigvals[i] -= fermi
    end
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
  if typeof(fermi)==String
    fermi = read_fermi_from_SCF_file(fermi)
  end
  out=deepcopy(band)
  for i1=1:size(band.eigvals)[1]
    out.k_points[i1,:] = band.k_points[i1,:]
    out.eigvals[i1] = band.eigvals[i1]-fermi
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