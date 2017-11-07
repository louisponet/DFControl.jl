# All the methods to change the inpût control flags, if you want to implement another kind of calculation add a similar one here!

"""
    change_flags!(input::QEInput, new_flag_data::Dict{Symbol,<:Any})

Changes the flags inside the input to the new ones if they are already defined and if the new ones have the same type.
"""
function change_flags!(input::QEInput, new_flag_data::Dict{Symbol,<:Any})
  found_keys = Symbol[]
  for block in input.control_blocks
    for (flag,value) in new_flag_data
      if haskey(block.flags,flag)
        old_data = block.flags[flag]
        if !(flag in found_keys) push!(found_keys,flag) end
        if typeof(block.flags[flag]) == typeof(new_flag_data[flag])
          block.flags[flag] = new_flag_data[flag]
          println("$(input.filename):\n -> $(block.name):\n  -> $flag:\n      $old_data changed to: $(new_flag_data[flag])")
          println("")
        else
          println("$(input.filename):\n -> $(block.name):\n  -> $flag:\n     type mismatch old:$old_data ($(typeof(old_data))), new: $(new_flag_data[flag]) ($(typeof(new_flag_data[flag])))\n    Change not applied.")
          println("")
        end
      end
    end
  end
  return found_keys
end

"""
    change_flags!(input::WannierInput, new_flag_data::Dict{Symbol,<:Any})

Changes the flags inside the input to the new ones if they are already defined and if the new ones have the same type.
"""
function change_flags!(input::WannierInput, new_flag_data::Dict{Symbol,<:Any})
  found_keys = Symbol[]
  for (flag,value) in new_flag_data
    if haskey(input.flags,flag)
      old_data = input.flags[flag]
      if !(flag in found_keys) push!(found_keys,flag) end
      if typeof(input.flags[flag]) == typeof(new_flag_data[flag])
        input.flags[flag] = new_flag_data[flag]
        println("$(input.filename):\n -> $flag:\n      $old_data changed to: $(new_flag_data[flag])")
        println("")
      else
        println("$(input.filename):\n -> $flag:\n     type mismatch old:$old_data ($(typeof(old_data))), new: $(new_flag_data[flag]) ($(typeof(new_flag_data[flag])))\n    Change not applied.")
        println("")
      end
    end
  end
  return found_keys
end

"""
    change_data!(input::DFInput, block_name::Symbol, new_block_data)

Changes the data of the specified 'DataBlock' to the new data if it has the correct type.
"""
function change_data!(input::DFInput, block_name::Symbol, new_block_data)
  for data_block in input.data_blocks
    if data_block.name == block_name
      if typeof(data_block.data) != typeof(new_block_data) 
        warn("Overwritten data of type '$(typeof(data_block.data))' with type '$(typeof(new_block_data))'.")
      end
      old_data = data_block.data
      data_block.data = new_block_data
      println("Block data '$(data_block.name)' in input  '$(input.filename)' is now:")
      display(data_block.data)
      println("")
    end
  end
end

"""
    get_flag(input::QEInput, flag::Symbol)

Returns the value of the flag.
"""
function get_flag(input::QEInput, flag::Symbol)
  for block in input.control_blocks
    if haskey(block.flags,flag)
      return block.flags[flag]
    end
  end
end

"""
    get_flag(input::WannierInput, flag::Symbol)

Returns the value of the flag.
"""
function get_flag(input::WannierInput, flag::Symbol)
  if haskey(input.flags,flag)
    return input.flags[flag]
  end
end

"""
    get_data(input::DFInput, block_symbol::Symbol)

Returns the specified 'DataBlock'.
"""
function get_data(input::DFInput, block_symbol::Symbol)
  block = get_block(input,block_symbol)
  if block != nothing && typeof(block) <: DataBlock
    return block.data
  else 
    error("No `DataBlock` with name '$block_symbol' found. ")
  end
end

"""
    get_block(input::QEInput, block_symbol::Symbol)

Returns the block with name `block_symbol`.
"""
function get_block(input::QEInput, block_symbol::Symbol)
  for block in [input.control_blocks; input.data_blocks]
    if block.name == block_symbol
      return block
    end
  end
  return nothing
end

"""
    get_block(input::WannierInput, block_symbol::Symbol)

Returns the block with name `block_symbol`.
"""
function get_block(input::WannierInput, block_symbol::Symbol)
  for block in input.data_blocks
    if block.name == block_symbol
      return block
    end
  end
  return nothing
end

#here comes the code for all the setting of flags of different inputs
"""
    add_flags!(input::QEInput, control_block_name::Symbol, flag_dict)

Adds the flags inside the dictionary to the 'ControlBlock'.
"""
function add_flags!(input::QEInput, control_block_name::Symbol, flag_dict)
  for block in input.control_blocks
    if block.name == control_block_name
      block.flags = merge((x,y) -> typeof(x) == typeof(y) ? y : x,block.flags,flag_dict)
      println("New input of block '$(block.name)' of calculation '$(input.filename)' is now:")
      display(block.flags)
      println("\n")
    end
  end
end

"""
    add_flags!(input::WannierInput, flag_dict)

Adds the flags inside the dictionary to the 'ControlBlock'.
"""
function add_flags!(input::WannierInput, flag_dict)
  input.flags = merge((x,y) -> typeof(x) == typeof(y) ? y : x,input.flags,flag_dict)
  println("New input of calculation '$(input.filename)' is now:")
  display(input.flags)
  println("\n")
end

#removes an input control flag, if you want to implement another input add a similar function here!
"""
    remove_flags!(input::QEInput, flags)

Remove the specified flags.
"""
function remove_flags!(input::QEInput, flags)
  for block in input.control_blocks
    if typeof(flags)<:Array{Symbol,1}
      for flag in flags
        if haskey(block.flags,flag)
          pop!(block.flags,flag)
          println("Removed flag '$flag' from block '$(block.name)' in input '$(input.filename)'")
        end
      end
    else
      if haskey(block.flags,flags)
        pop!(block.flags,flags)
        println("Removed flag '$flags' from block '$(block.name)' in input '$(input.filename)'")
      end
    end
  end
end

"""
    remove_flags!(input::WannierInput, flags)

Remove the specified flags.
"""
function remove_flags!(input::WannierInput, flags)
  if typeof(flags) <: Array{Symbol,1}
    for flag in flags
      if haskey(input.flags,flag)
        pop!(input.flags,flag,false)
        println("Removed flag '$flag' from input '$(input.filename)'")
      end
    end
  else
    if haskey(input.flags,flags)
      pop!(input.flags,flags,false)
      println("Removed flag '$flags' from input '$(input.filename)'")
    end
  end
end

"""
    get_blocks(input::DFInput) 

Returns all the blocks inside a calculation.
"""
function get_blocks(input::DFInput) 
  out = Block[]
  for block in  getfield.(input,filter(x->contains(String(x),"block"),fieldnames(input)))
    push!(out,block...)
  end
  return out 
end

"""
    print_block(input::DFInput, block_name::Symbol)

Print the information of a 'Block' inside the input.
"""
function print_block(input::DFInput, block_name::Symbol)
  input_blocks = get_blocks(input)
  found = false
  for block in input_blocks
    if block.name == block_name
      print_filename(input)
      display(block) 
      found = true
    end
  end
  return found
end

"""
    print_blocks(input::DFInput)

Print the information on all Blocks inside the input.
"""
function print_blocks(input::DFInput)
  print_filename(input)
  get_blocks(input) |> display
  println("")
end

"""
    print_filename(input::DFInput)

Prints the filename associated with the input.
"""
function print_filename(input::DFInput)
  println("Input file: $(input.filename)")
end

"""
    print_info(input::DFInput)

Prints general info of the input.
"""
function print_info(input::DFInput)
  println("Filename: $(input.filename)")
  if (:control_blocks in fieldnames(input))
    println("  Control Blocks:")
    for (i,block) in enumerate(input.control_blocks)
      println("    $i: $(block.name)")
    end
  end
  println("  Data Blocks:")
  for (i,block) in enumerate(input.data_blocks)
    println("    $i: $(block.name)")
  end
  println("  Run command: $(input.run_command)")
  println("  Runs: $(input.run)")
end

"""
    print_flags(input::DFInput)

Prints all the flags of the input.
"""
function print_flags(input::DFInput)
  println("Filename: $(input.filename)")
  if (:control_blocks in fieldnames(input))
    for block in input.control_blocks
      println("  $(block.name):")
      for (flag,value) in block.flags
        println("    $flag => $value")
      end
    end
  end
  if (:flags in fieldnames(input))
    for (flag,value) in input.flags
      println("  $flag => $value")
    end
  end
end

"""
    print_flag(input::DFInput, flag)

Prints information of the flag inside the input.
"""
function print_flag(input::DFInput, flag)
  if (:control_blocks in fieldnames(input))
    for block in input.control_blocks
      if haskey(block.flags,flag)
        println("Filename: $(input.filename)")
        println("  Block Name: $(block.name)")
        println("    $flag => $(block.flags[flag])")
        println("")
      end
    end
  end
  if (:flags in fieldnames(input))
    if haskey(input.flags,flag)
      println("Filename: $(input.filename)")
      println("  $flag => $(input.flags[flag])")
      println("")
    end
  end
end

#finish this macro
# macro change_data_switch(switches::Array{Tuple{Type,Symbol},1})
#   return quote 
#     if typeof(esc(input)) == switches[1][1]
#       change_data!(esc(input),switches[1][2],)
#in all of these functions, we use default names, can again be shorter...
#we only implement fractional atoms stuff does anyone even use anything else?
#Add more calculations here!
"""
    change_atoms!(input::DFInput, atoms::Dict{Symbol,<:Array{<:Point3D,1}}, pseudo_set_name=:default; pseudo_fuzzy = nothing)

Changes the data inside the 'DataBlock' that holds the data of the atom positions.
If 'default_pseudos' is defined it will look for the pseudo set and also write the correct values inside the 'DataBlock' that defines which pseudopotentials to use.
"""
function change_atoms!(input::DFInput, atoms::Dict{Symbol,<:Array{<:Point3D,1}}, pseudo_set_name=nothing,pseudo_fuzzy = nothing)
  if typeof(input) == WannierInput
    change_data!(input,:atoms_frac,atoms)
  elseif typeof(input) == QEInput
    change_data!(input,:atomic_positions,atoms)
    if isdefined(:default_pseudos) && pseudo_set_name != nothing
      atomic_species_dict = Dict{Symbol,String}()
      for atom in keys(atoms)
        atomic_species_dict[atom] = get_default_pseudo(atom,pseudo_set_name,pseudo_fuzzy=pseudo_fuzzy)
      end
      change_data!(input,:atomic_species,atomic_species_dict)
    end
    if isdefined(:default_pseudo_dirs) && pseudo_set_name != nothing
      change_flags!(input,Dict(:pseudo_dir => "'$(default_pseudo_dirs[pseudo_set_name])'"))
    end
  end
end

"""
    change_cell_parameters!(input::DFInput, cell_param::Array{AbstractFloat,2})

Changes the cell parameters `DataBlock`.
"""
function change_cell_parameters!(input::DFInput, cell_param::Array{<:AbstractFloat,2})
  assert(size(cell_param)==(3,3))
  if typeof(input) == WannierInput
    change_data!(input,:unit_cell_cart,cell_param)
  elseif typeof(input) == QEInput
    change_data!(input,:cell_parameters,cell_param)
  end
end

"""
    change_k_points!(input::DFInput,calc_filename,k_points)

Changes the data in the k point `DataBlock` inside the specified calculation.
"""
function change_k_points!(input::DFInput,k_points)
  if typeof(input) == WannierInput
    change_data!(input,:kpoints,k_points)
  elseif typeof(input) == QEInput
    change_data!(input,:k_points,k_points)
  end
end