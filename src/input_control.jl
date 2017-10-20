# All the methods to change the inpût control flags, if you want to implement another kind of calculation add a similar one here!
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

function get_flag(input::QEInput,flag::Symbol)
  for block in input.control_blocks
    if haskey(block.flags,flag)
      return block.flags[flag]
    end
  end
end

function get_flag(input::WannierInput,flag::Symbol)
  if haskey(input.flags,flag)
    return input.flags[flag]
  end
end

function get_data(input::DFInput,block_symbol::Symbol)
  for block in input.data_blocks
    if block.name == block_symbol
      return block.data
    end
  end
end

#here comes the code for all the setting of flags of different inputs
function set_flags!(input::QEInput, control_block_name::Symbol, flag_dict)
  for block in input.control_blocks
    if block.name == control_block_name
      block.flags = merge((x,y) -> typeof(x) == typeof(y) ? y : x,block.flags,flag_dict)
      println("New input of block '$(block.name)' of calculation '$(input.filename)' is now:")
      display(block.flags)
      println("\n")
    end
  end
end

function set_flags!(input::WannierInput, flag_dict)
  input.flags = merge((x,y) -> typeof(x) == typeof(y) ? y : x,input.flags,flag_dict)
  println("New input of calculation '$(input.filename)' is now:")
  display(input.flags)
  println("\n")
end

#removes an input control flag, if you want to implement another input add a similar function here!
function remove_flags!(input::QEInput,flags)
  for block in input.control_blocks
    if typeof(flags)<:Array{Symbol,1}
      for flag in flags
        if haskey(block.flags,flag)
          pop!(block.flags,flag)
          println("Removed flag '$flag' from block '$(block.name)' in input '$(input.filename)'")
        end
      end
    else
      if haskey(block.flags,flag)
        pop!(block.flags,flag)
        println("Removed flag '$flag' from block '$(block.name)' in input '$(input.filename)'")
      end
    end
  end
end

function remove_flags!(input::WannierInput,flags)
  if typeof(flags) <: Array{Symbol,1}
    for flag in flags
      if haskey(input.flags,flag)
        pop!(input.flags,flag,false)
        println("Removed flag '$flag' from input '$(input.filename)'")
      end
    end
  else
    if haskey(input.flags,flag)
      pop!(input.flags,flag,false)
      println("Removed flag '$flag' from input '$(input.filename)'")
    end
  end
end

function print_block(input::DFInput, block_name::Symbol)
  input_blocks = getfield.(input,filter(x->contains(String(x),"block"),fieldnames(input)))
  found = false
  for blocks in input_blocks
    for block in blocks
      if block.name == block_name
        println("Input file: $(input.filename)")
        display(block) 
        found = true
      end
    end
  end
  return found
end

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

function print_flag(input::DFInput,flag)
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
      println("  $flag => $(block.flags[flag])")
      println("")
    end
  end
end