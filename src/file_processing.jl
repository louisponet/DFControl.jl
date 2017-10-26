function parse_k_line(line,T)
  splt = split(line)
  k1   = parse(T,splt[5])
  k2   = parse(T,splt[6])
  k3   = parse(T,splt[7][1:1:end-2])
  return [k1,k2,k3]
end
function write_flag_line(f,flag,data)
  write(f,"   $flag = ")
  if typeof(data) <: Array
    write(f,"$(data[1])")
    for x in data[2:end]
      write(f," $x")
    end
    write(f,"\n")
  else #this should work for anything singular valued data such as bools, ''s and other types
    write(f,"$data\n")
  end
end
"""
    read_qe_bands_file(filename::String, T=Float32)

Reads the output file of a 'bands' calculation in Quantum Espresso.
Returns an array of DFBands each with the same k_points and their respective energies.
"""
function read_qe_bands_file(filename::String, T=Float32)
  f = open(filename)
  k_points_cart  = Array{Array{T,1},1}()
  k_points_cryst = Array{Array{T,1},1}()
  k_eigvals      = Array{Array{T,1},1}()
  while !eof(f)
    line = readline(f)
    if contains(line,"celldm(1)")
      alat_bohr = parse(T,split(line)[2])
      prefac    = T(2pi/alat_bohr * 1.889725)
    end
    if contains(line,"cryst.") && length(split(line))==2
      line = readline(f)
      while length(line) != 0
        push!(k_points_cryst,parse_k_line(line,T))
        line = readline(f)
      end
    end
    if contains(line,"cart.") && length(split(line)) == 5
      line = readline(f)
      while line != ""
        push!(k_points_cart,prefac*parse_k_line(line,T))
        line = readline(f)
      end
    end
    if contains(line,"k") && contains(line,"PWs)")
      tmp = T[]
      readline(f)
      line = readline(f)
      while line != ""
        append!(tmp,parse_line(T,line))
        line = readline(f)
      end
      push!(k_eigvals,tmp)
    end
  end
  out = Array{DFBand{T},1}()
  for i1=1:length(k_eigvals[1])
    eig_band = T[]
    for i=1:size(k_points_cart)[1]
      push!(eig_band,k_eigvals[i][i1])
    end
    push!(out,DFBand(k_points_cart,k_points_cryst,eig_band))
  end
  close(f)
  return out
end

"""
    read_ks_from_qe_bands_file(filename::String,T=Float32)

Read k-points from a Quantum Espresso bands output file in cartesian (2pi/alat in Angstrom^-1!) and crystalline coordinates.
Returns (k_points_cart,k_points_cryst).
"""
function read_ks_from_qe_bands_file(filename::String,T=Float32)
  open(filename) do f
    out_cryst = Array{Array{T,1},1}()
    out_cart  = Array{Array{T,1},1}()
    alat_bohr = zero(T)
    prefac    = zero(T)
    while !eof(f)
      line = readline(f)
      if contains(line,"celldm(1)")
        alat_bohr = parse(T,split(line)[2])
        prefac    = T(2pi/alat_bohr * 1.889725)
      elseif contains(line,"cryst.") && length(split(line))==2
        line = readline(f)
        while line != ""
          push!(out_cryst,parse_k_line(line,T))
          line = readline(f)
        end
      elseif contains(line,"cart.") && length(split(line)) == 5
        line = readline(f)
        while line != ""
          push!(out_cart,prefac*parse_k_line(line,T))
          line = readline(f)
        end
      end
    end
    return out_cart,out_cryst
  end
end

"""
    read_fermi_from_qe_file(filename::String,T=Float32)

Reads the Fermi level from a Quantum Espresso scf calculation output file
(if there is one).
"""
function read_fermi_from_qe_file(filename::String,T=Float32)
  out = nothing
  open(filename) do f
    while !eof(f)
      line = split(readline(f))
      if "Fermi" in line
        out= parse(T,line[5])
      elseif "lowest" in line &&  "unoccupied" in line
        out = Dict(:lowest_unoccupied => parse(T,line[end]),:highest_occupied => parse(T,line[end-1]))
      end
    end
  end
  if out==nothing
    error("Couldn't find the Fermi level in file $filename. Is this an scf output file?")
  else
    return out
  end
end

"""
    read_qe_kpdos(filename::String,column=1;fermi=0)

Reads the k_resolved partial density of states from a Quantum Espresso projwfc output file.
Only use this if the flag kresolveddos=true in the projwfc input file!!

The returned matrix can be readily plotted using heatmap() from Plots.jl!

Optional input: column = 1 (column of the output, 1 = first column after ik and E)
                fermi  = 0 (possible fermi offset of the read energy values)
Return:         Array{Float64,2}(length(k_points),length(energies)) ,
                (ytickvals,yticks)
"""
function read_qe_kpdos(filename::String,column=1;fermi=0)
  read_tmp = readdlm(filename)
  zmat = zeros(typeof(read_tmp[1]),Int64(read_tmp[end,1]),size(read_tmp)[1]/Int64(read_tmp[end,1]))
  for i1=1:size(zmat)[1]
    for i2=1:size(zmat)[2]
      zmat[i1,i2] = read_tmp[size(zmat)[2]*(i1-1)+i2,2+column]
    end
  end

  yticks = collect(Int64(div(read_tmp[1,2]-fermi,1)):1:Int64(div(read_tmp[end,2]-fermi,1)))
  ytickvals=[findfirst(x->norm(yticks[1]+fermi-x)<=0.1,read_tmp[:,2])]
  for (i,tick) in enumerate(yticks[2:end])
    push!(ytickvals,findnext(x->norm(tick+fermi-x)<=0.1,read_tmp[:,2],ytickvals[i]))
  end
  return  zmat',(ytickvals,yticks)
end



"""
    read_qe_input(filename,T=Float32)

Reads a Quantum Espresso input file.
Returns a DFInput.
"""
function read_qe_input(filename,T=Float32;run_command ="",run=true)
  function get_card_option(line,length)
    if contains(line,"{")
      return Symbol(strip(split(line,"{")[end],'}'))
    elseif contains(line,"(")
      return Symbol(strip(split(line,"(")[end],')'))
    end
  end
  control_blocks = Array{QEControlBlock,1}()
  data_blocks    = Array{QEDataBlock,1}()

  open(filename) do f
    line = readline(f)
    while !eof(f)
      @label start_label
      if contains(line,"&")
        c_block_name = Symbol(lowercase(strip(strip(line),'&')))
        flag_dict = Dict{Symbol,Any}()
        line = readline(f)
        while strip(line)!="/"
          if contains(line,"!")
            line = readline(f)
            continue
          end
          split_line = filter(x->x!="",strip.(split(line,",")))
          for s in split_line
            key,val = String.(strip.(split(s,"=")))
            #parse for ints
            if !isnull(tryparse(Int,val))
              flag_dict[Symbol(key)] = get(tryparse(Int,val))
            #parse for T<:AbstractFloat
            elseif !isnull(tryparse(T,val))
              flag_dict[Symbol(key)] = get(tryparse(T,val))
            #parse for Bools
            elseif !isnull(tryparse(Bool,strip(val,'.')))
              flag_dict[Symbol(key)] = get(tryparse(Bool,strip(val,'.')))
            else
              flag_dict[Symbol(key)] = val
            end
          end
          line = readline(f)
        end
        push!(control_blocks,QEControlBlock(c_block_name,flag_dict))
        @goto start_label

      elseif contains(line,"CELL_PARAMETERS")||contains(line,"cell_parameters")
        cell_unit = get_card_option(line,15)
        cell = Matrix{T}(3,3)
        cell[1,1:3] = parse.(T,split(readline(f)))
        cell[2,1:3] = parse.(T,split(readline(f)))
        cell[3,1:3] = parse.(T,split(readline(f)))
        line = readline(f)
        push!(data_blocks,QEDataBlock(:cell_parameters,cell_unit,cell))
        @goto start_label

      elseif contains(line,"ATOMIC_SPECIES")||contains(line,"atomic_species")
        line   = readline(f)
        pseudos = Dict{Symbol,String}()
        while length(split(line))==3
          pseudos[Symbol(split(line)[1])] = split(line)[end]
          line = readline(f)
        end
        push!(data_blocks,QEDataBlock(:atomic_species,:none,pseudos))
        @goto start_label

      elseif contains(line,"ATOMIC_POSITIONS")||contains(line,"atomic_positions")
        option = get_card_option(line,16)
        atoms  = Dict{Symbol,Array{Point3D{T},1}}()
        line   = readline(f)
        while length(split(line))==4
          s_line   = split(line)
          atom     = Symbol(s_line[1])
          position = Point3D(parse(T,s_line[2]),parse(T,s_line[3]),parse(T,s_line[4]))
          if !haskey(atoms,atom)
            atoms[atom]=[position]
          else
            push!(atoms[atom],position)
          end
          line = readline(f)
        end
        push!(data_blocks,QEDataBlock(:atomic_positions,option,atoms))
        @goto start_label

      elseif contains(line,"K_POINTS")||contains(line,"k_points")
        k_option = get_card_option(line,8)
        line     = readline(f)
        if k_option == :automatic
          s_line = split(line)
          k_data = parse.(Int,s_line)
        else
          nks = parse(Int,line)
          k_data = Array{Array{T,1},1}(nks)
          for i=1:nks
            k_data[i] = parse.(T,split(readline(f)))
          end
        end
        push!(data_blocks,QEDataBlock(:k_points,k_option,k_data))
        @goto start_label
      end
      line = readline(f)
    end
  end
  return QEInput(splitdir(filename)[2],control_blocks,data_blocks,run_command,run)
end

#can I use @generated here?
function write_block_data(f,data)
  if typeof(data)<:Array{Vector{Float32},1} #k_points
    for x in data
      for y in x
        write(f," $y")
      end
      write(f,"\n")
    end
  elseif typeof(data) <:Array{Int,1}
    for x in data
      write(f," $x")
    end
    write(f,"\n")
  elseif typeof(data) <: Matrix
    im,jm = size(data)
    for i=1:im
      for j=1:jm
        write(f," $(data[i,j])")
      end
      write(f,"\n")
    end
  elseif typeof(data)<:Dict{Symbol,<:Any}
    for (key,value) in data
      if typeof(value) == String
        write(f,"$key $(Float32(ELEMENTS[key].atomic_weight))   $value\n")
      elseif typeof(value) <: Array{<:Point3D,1}
        for at in value
          write(f,"$key $(at.x) $(at.y) $(at.z)\n")
        end
      end
    end
  end
end

"""
    write_qe_input(filename::String,df_input::DFInput)

Writes a Quantum Espresso input file.
"""
function write_qe_input(input::QEInput,filename::String=input.filename)
  open(filename,"w") do f
    write_flag(flag_data) = write_flag_line(f,flag_data[1],flag_data[2])
    write_block(data)      = write_block_data(f,data)
    for block in input.control_blocks
      write(f,"&$(block.name)\n")
      # block.flags |> write_flag
      map(write_flag,[(flag,data) for (flag,data) in block.flags])
      write(f,"/\n\n")
    end

    for block in input.data_blocks
      if block.option != :none
        write(f,"$(uppercase(String(block.name))) ($(block.option))\n")
      else
        write(f,"$(uppercase(String(block.name)))\n")
      end
      if block.name == :k_points && block.option != :automatic
        write(f,"$(length(block.data))\n")
        write_block(block.data)
      else
        write_block(block.data)
      end
      write(f,"\n")
      #write the other cards depending on whether or not they are there
    end
  end
end
#---------------------------END QUANTUM ESPRESSO SECTION----------------#
#---------------------------START WANNIER SECTION ----------------------#

"""
   read_wannier_input(filename::String, T=Float32)

Reads a DFInput from a wannier90 input file.
"""
function read_wannier_input(filename::String, T=Float32; run_command="", run=true, preprocess=true)
  flags = Dict{Symbol,Any}()
  data_blocks = Array{WannierDataBlock,1}()
  open(filename,"r") do f
    line = readline(f)
    while !eof(f)
      @label start_label
      if contains(line,"!") || line=="" || contains(lowercase(line),"end")
        line = readline(f)
        continue
      end
      if contains(lowercase(line),"begin")
        block_name = Symbol(split(lowercase(line))[end])

        if block_name == :projections
          line = readline(f)
          if line == "random"
            push!(data_blocks,WannierDataBlock(:projections,:random,nothing))
            @goto start_label
          else
            proj_dict = Dict{Symbol,Array{Symbol,1}}()
            while !contains(lowercase(line),"end")
              if contains(line,"!") || line == ""
                line = readline(f)
                continue
              end
              split_line = strip_split(line,':')
              atom = Symbol(split_line[1])
              projections = [Symbol(proj) for proj in strip_split(split_line[2],';')]
              proj_dict[atom] = projections
              line = readline(f)
            end
            push!(data_blocks,WannierDataBlock(:projections,:none,proj_dict))
            @goto start_label
          end

        elseif block_name == :kpoint_path
          line = readline(f)
          k_path_array = Array{Tuple{Symbol,Array{T,1}},1}()
          while !contains(lowercase(line),"end")
            if contains(line,"!")
              line = readline(f)
              continue
            end
            split_line = split(line)
            push!(k_path_array,(Symbol(split_line[1]),parse_string_array(T,split_line[2:4])))
            push!(k_path_array,(Symbol(split_line[5]),parse_string_array(T,split_line[6:8])))
            line = readline(f)
          end
          push!(data_blocks, WannierDataBlock(:kpoint_path,:none,k_path_array))
          @goto start_label

        elseif block_name == :unit_cell_cart
          line = readline(f)
          if length(split(line)) == 1
            option = Symbol(lowercase(line))
          else
            option = :ang
          end
          cell_param = Matrix{T}(3,3)
          for i=1:3
            cell_param[i,:] = parse_line(T,readline(f))
          end
          push!(data_blocks, WannierDataBlock(:unit_cell_cart,option,cell_param))
          line = readline(f)
          @goto start_label

        elseif block_name ==:atoms_frac || block_name == :atoms_cart
          line = readline(f)
          atoms = Dict{Symbol,Array{Point3D{T},1}}()
          option = Symbol(split(String(block_name),"_")[end])
          while !contains(lowercase(line),"end")
            split_line = strip_split(line)
            atom = Symbol(split_line[1])
            position = Point3D(parse_string_array(T,split_line[2:4]))
            if !haskey(atoms,atom)
              atoms[atom] = [position]
            else
              push!(atoms[atom],position)
            end
            line = readline(f)
          end
          push!(data_blocks, WannierDataBlock(block_name,option,atoms))
          @goto start_label

        elseif block_name == :kpoints
          line = readline(f)
          k_points = Array{Array{T,1},1}()
          while !contains(lowercase(line),"end")
            if line==""
              line = readline(f)
              continue
            end
            push!(k_points,parse_line(T,line))
            line = readline(f)
          end
          push!(data_blocks, WannierDataBlock(:kpoints,:none,k_points))
          @goto start_label
        end

      else
        if contains(line,"mp_grid")
          flags[:mp_grid] = parse_string_array(Int,split(split(line,'=')[2]))
        else
          split_line = strip_split(line,'=')
          flag = Symbol(split_line[1])
          value = split_line[2]
          if  lowercase(value)=="t" || lowercase(value) == "true"
            flags[flag] = true
          elseif lowercase(value) == "f" || lowercase(value) == "false"
            flags[flag] = false
          elseif !isnull(tryparse(Int,value))
            flags[flag] = get(tryparse(Int,value))
          elseif !isnull(tryparse(T,value))
            flags[flag] = get(tryparse(T,value))
          else
            flags[flag] = value
          end
        end
      end
      line = readline(f)
    end
  end
  return WannierInput(splitdir(filename)[2],flags,data_blocks,run_command,run,preprocess)
end

"""
    write_wannier_input(filename::String,input::DFInput)

Writes the input of a wannier90 input file.
"""
function write_wannier_input(input::WannierInput,filename::String=input.filename)
  open(filename,"w") do f
    for (flag,value) in input.flags
      write_flag_line(f,flag,value)
    end
    write(f,"\n")
    for block in input.data_blocks
      write(f,"begin $(block.name)\n")
      if block.name == :kpoint_path
        for i = 1:2:length(block.data)
          letter1, k_points1 = block.data[i]
          letter2, k_points2 = block.data[i+1]
          write(f,"$letter1 $(k_points1[1]) $(k_points1[2]) $(k_points1[3]) $letter2 $(k_points2[1]) $(k_points2[2]) $(k_points2[3])\n")
        end
        # write(f,"\n")

      elseif block.name == :projections
        for (atom,symbols) in block.data
          write(f,"$atom: $(symbols[1])")
          for sym in symbols[2:end]
            write(f,";$sym")
          end
          write(f,"\n")
        end
        write(f,"\n")

      elseif block.name == :unit_cell_cart
        matrix = block.data
        write(f,"$(block.option)\n")
        write(f,"$(matrix[1,1]) $(matrix[1,2]) $(matrix[1,3])\n")
        write(f,"$(matrix[2,1]) $(matrix[2,2]) $(matrix[2,3])\n")
        write(f,"$(matrix[3,1]) $(matrix[3,2]) $(matrix[3,3])\n")

      elseif block.name == :atoms_frac || block.name == :atoms_cart
        for (atom,points) in block.data
          for point in points
            write(f,"$atom $(point.x) $(point.y) $(point.z)\n")
          end
        end

      elseif block.name == :kpoints
        for k in block.data
          write(f,"$(k[1]) $(k[2]) $(k[3])\n")
        end
      end
      write(f,"end $(block.name)\n\n")
    end
  end
end
#---------------------------END WANNIER SECTION ------------------------#
#---------------------------BEGINNING GENERAL SECTION-------------------#
"""
    write_input(filename::String,df_input::DFInput)

Writes the input file for a DFInput.
Backend of DFInput decides what writing function is called.
"""
function write_input(df_input::DFInput,filename::String=df_input.filename)
  if typeof(df_input) == QEInput
    write_qe_input(df_input,filename)
  elseif typeof(df_input) == WannierInput
    write_wannier_input(df_input,filename)
  end
end

"""
    write_job_files(df_job::DFJob)

Writes all the input files and job file that are linked to a DFJob.
"""
function write_job_files(df_job::DFJob)
  files_to_remove = search_dir(df_job.local_dir,".in")
  new_filenames   = String[]

  open(df_job.local_dir*"job.tt","w") do f
    write(f,"#!/bin/bash\n","#SBATCH -N 1\n","#SBATCH --ntasks-per-node=24 \n","#SBATCH --time=24:00:00 \n","#SBATCH -J $(df_job.name) \n",
          "#SBATCH -p defpart\n\n","module load open-mpi/gcc/1.10.2\n","module load mkl/2016.1.056\n","\n")

    for i=1:length(df_job.calculations)
      calculation = df_job.calculations[i]
      run_command = calculation.run_command
      filename = calculation.filename
      write_input(calculation,df_job.local_dir*filename)
      should_run  = calculation.run
      if typeof(calculation) == WannierInput
        if calculation.preprocess
          run_command *= " -pp"
        end
        if !should_run
          write(f,"#$run_command $(filename[1:end-4])\n")
        else
          write(f,"$run_command $(filename[1:end-4])\n")
        end
      else
        if !should_run
          write(f,"#$run_command <$filename> $(split(filename,".")[1]).out \n")
        else
          write(f,"$run_command <$filename> $(split(filename,".")[1]).out \n")
        end
      end
      push!(new_filenames,filename)
    end
  end

  for file in files_to_remove
    if !in(file,new_filenames)
      rm(df_job.local_dir*file)
    end
  end
end

"""
    read_job_file(job_file::String)

Reads and returns the name, input files, run_commands and whether or not they need to be commented out.
All files that are read contain "in".
This reads QE and wannier90 inputs for now.
"""
function read_job_file(job_file::String)
  input_files  = Array{String,1}()
  output_files = Array{String,1}()
  run_commands = Array{String,1}()
  should_run   = Array{Bool,1}()
  name     = ""
  open(job_file,"r") do f
    while !eof(f)
      line = readline(f)
      if line == ""
        continue
      end
      if contains(line,".x ")
        if contains(line,"#")
          push!(should_run,false)
          line = line[2:end]
        else
          push!(should_run,true)
        end

        s_line = split(line)
        i = 0
        for (j,s) in enumerate(s_line)
          if contains(s,".x")
            i = j
            break
          end
        end
        run_command = prod([s*" " for s in s_line[1:i]])
        if contains(s_line[i+1],"-")
          run_command *= s_line[i+1]
          i+=1
        end
        push!(run_commands,run_command)
        #handles QE and Wannier.
        in_out = strip_split(prod(s_line[i+1:end]),">")
        if length(in_out) == 2
          push!(input_files,strip(in_out[1],'<'))
          push!(output_files,in_out[2])
        else
          push!(input_files,strip(in_out[1],'<'))
        end 
      elseif contains(line,"#SBATCH") && contains(line,"-J")
        name = split(line)[end]
      end
    end
  end
  return name,input_files,output_files,run_commands,should_run
end

#---------------------------END GENERAL SECTION-------------------#

#Incomplete: we should probably handle writing an array of expressions as well!
function expr2file(filename::String, expression::Expr)
  eq        = Symbol("=")
  lines     = readlines(filename)
  new_lines = String[]
  found     = false

  if expression.head != eq
    error("For now only writing of assignment expressions is possible.")
  end

  lhs = expression.args[1]
  rhs = expression.args[2]

  for line in lines
    if line == ""
      continue
    end
    
    expr = parse(line)
    if expr.head != eq
      continue
    end
    
    lhs_t = expr.args[1]
    rhs_t = expr.args[2]
    if lhs == :(default_pseudos[:Ge])
    end

    if lhs_t == lhs
      found = true
      push!(new_lines,"$(:($lhs = $rhs))")
    else
      push!(new_lines,"$expr")
    end
  end

  open(filename,"w") do f
    for line in new_lines
      write(f,line*"\n")
    end
    if !found
      write(f,"$expression\n")
    end
  end
end

function init_defaults(filename::String)
  raw_input =""
  names_to_export = Symbol[] 
  open(filename,"r") do f
    while !eof(f)
      line = readline(f)
      lhs = parse(line).args[1]
      if typeof(lhs) == Symbol
        push!(names_to_export,lhs)
      end
      raw_input *= line*"; "
    end
  end
  for name in names_to_export
    eval(:(export $name))
  end
  eval(parse(raw_input))
end

function load_defaults(filename::String)
  raw_input =""
  names_to_export = Symbol[] 
  open(filename,"r") do f
    while !eof(f)
      raw_input *= readline(f)*"; "
    end
  end
  raw_input *= "nothing ;"
  eval(parse(raw_input))
end



