function parse_k_line(line,T)
  splt = split(line)
  k1   = parse(T,splt[5])
  k2   = parse(T,splt[6])
  k3   = parse(T,splt[7][1:1:end-2])
  return [k1,k2,k3]
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
  out = zero(T)
  open(filename) do f
    while !eof(f)
      line = split(readline(f))
      if "Fermi" in line
        out= parse(T,line[5])
      end
    end
  end
  if out==zero(T)
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
function read_qe_input(filename,T=Float32)
  function get_card_option(line,length)
    return Symbol(strip(strip(line[length+1:end]),['(',')','{','}']))
  end
  control_blocks = Dict{Symbol,Dict{Symbol,Any}}()
  cell_param     = Dict{Symbol,Any}()
  pseudos        = Dict{Symbol,String}()
  atoms          = Dict{Symbol,Array{<:Point3D,1}}()
  k_dict         = Dict{Symbol,Any}()

  open(filename) do f
    while !eof(f)
      line = readline(f)
      @label startlabel
      if contains(line,"&")
        block_symbol = Symbol(lowercase(strip(strip(line),'&')))
        tmp_dict = Dict{Symbol,Any}()
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
              tmp_dict[Symbol(key)] = get(tryparse(Int,val))
            #parse for T<:AbstractFloat
            elseif !isnull(tryparse(T,val))
              tmp_dict[Symbol(key)] = get(tryparse(T,val))
            #parse for Bools
            elseif !isnull(tryparse(Bool,strip(val,'.')))
              tmp_dict[Symbol(key)] = get(tryparse(Bool,strip(val,'.')))
            else
              tmp_dict[Symbol(key)] = val
            end
          end
          line = readline(f)
        end
        control_blocks[block_symbol]=tmp_dict
        @goto startlabel


      elseif contains(line,"CELL_PARAMETERS")||contains(line,"cell_parameters")
        cell_unit = get_card_option(line,15)
        cell = Matrix{T}(3,3)
        cell[1,1:3] = parse.(T,split(readline(f)))
        cell[2,1:3] = parse.(T,split(readline(f)))
        cell[3,1:3] = parse.(T,split(readline(f)))
        cell_param  = Dict(cell_unit => cell)
        line = readline(f)
        @goto startlabel

      elseif contains(line,"ATOMIC_SPECIES")||contains(line,"atomic_species")
        line = readline(f)
        while length(split(line))==3
          pseudos[Symbol(split(line)[1])] = split(line)[end]
          line = readline(f)
        end
        @goto startlabel

#TODO fix that atoms can only be fractional!
      elseif contains(line,"ATOMIC_POSITIONS")||contains(line,"atomic_positions")
        pos_unit = get_card_option(line,16)
        line = readline(f)
        while length(split(line))==4
          line_segs = split(line)
          atom = Symbol(line_segs[1])
          position=Point3D(parse(T,line_segs[2]),parse(T,line_segs[3]),parse(T,line_segs[4]))
          if !haskey(atoms,atom)
            atoms[atom]=[position]
          else
            push!(atoms[atom],position)
          end
          line = readline(f)
        end
        @goto startlabel

      elseif contains(line,"K_POINTS")||contains(line,"k_points")
        k_option = get_card_option(line,8)
        line     = readline(f)
        k_dict[:option] = k_option
        if k_option == :automatic
          line_segs = split(line)
          k_dict[:nk1],k_dict[:nk2],k_dict[:nk3],k_dict[:sk1],k_dict[:sk2],k_dict[:sk3] = parse.(Int,line_segs)
        else
          nks = parse(Int,line)
          k_point_array = Array{Array{T,1},1}(nks)
          for i=1:nks
            k_point_array[i] = parse.(T,split(readline(f)))
          end
          k_dict[:points] = k_point_array
        end
        @goto startlabel
      end

    end
  end
  return DFInput(:QE,control_blocks,pseudos,cell_param,atoms,k_dict)
end


"""
    write_qe_input(filename::String,df_input::DFInput)

Writes a Quantum Espresso input file.
"""
function write_qe_input(filename::String,df_input::DFInput)
  if df_input.backend != :QE
    error("Tried to write a QE input file with a DfInput of backend: $(df_input.backend).")
  else
    open(filename,"w") do f
      write(f,"!Generated from DFTools at $(now())\n")

      # write the control blocks first
      for (key,value) in df_input.control_blocks
        write(f,"&$(String(key))\n")
        for (input_var,var_val) in df_input.control_blocks[key]
          write(f,"  $(String(input_var)) = $(var_val)\n")
        end
        write(f,"/\n\n")
      end

      #write the other cards depending on whether or not they are there
      if !isempty(df_input.atoms)
        write(f,"ATOMIC_SPECIES\n")
        for atom in keys(df_input.atoms)
    #@Incomplete Once default pseudos are implemented put it in the default slot here!!!!!
          pseudo = get(df_input.pseudos,atom,"PUT PSEUDO HERE")
          write(f,"$(String(atom))   $(Float32(ELEMENTS[atom].atomic_weight))   $pseudo\n")
        end
        write(f,"\n")
      end

      if !isempty(df_input.cell_param)
        cell_format = collect(keys(df_input.cell_param))[1]
        cell = df_input.cell_param[cell_format]
        write(f,"CELL_PARAMETERS ($cell_format)\n")
        write(f,"$(cell[1,1])  $(cell[1,2]) $(cell[1,3])\n")
        write(f,"$(cell[2,1])  $(cell[2,2]) $(cell[2,3])\n")
        write(f,"$(cell[3,1])  $(cell[3,2]) $(cell[3,3])\n")
        write(f,"\n")
      end

      if !isempty(df_input.atoms)
        write(f,"ATOMIC_POSITIONS (crystal)\n")
        for (atom,points) in df_input.atoms
          for point in points
            write(f,"$(String(atom))    $(point.x)    $(point.y)    $(point.z) \n")
          end
        end
        write(f,"\n")
      end

      if !isempty(df_input.k_points)
        k_option = df_input.k_points[:option]
        write(f,"K_POINTS ($k_option)\n")

        if k_option == :automatic
          write(f,"$(df_input.k_points[:nk1])  $(df_input.k_points[:nk2])  $(df_input.k_points[:nk2])  $(df_input.k_points[:sk1])  $(df_input.k_points[:sk2])  $(df_input.k_points[:sk3])\n")
        elseif k_option == :crystal_b
          points = df_input.k_points[:points]
          write(f,"$(size(points)[1])\n")
          for point in points
            write(f,"$(point[1])  $(point[2])  $(point[3]) $(Int(point[4]))\n")
          end
        else
          points = df_input.k_points[:points]
          write(f,"$(size(points)[1])\n")
          for point in points
            write(f,"$(point[1])  $(point[2])  $(point[3]) $(point[4])\n")
          end
        end

      end
    end
  end
end

#@Incomplete: does not read wan input files now
"""
    read_qe_inputs_from_job_file(job_file)

Reads a job file and all the mentioned Quantum Espresso input files.
Returns an array of String pairs that are (calculation,input_file)
"""
function read_qe_inputs_from_job_file(job_file)
  inputs = Array{Tuple{String,String},1}()
  open(job_file,"r") do f
    while !eof(f)
      line = readline(f)
      if contains(line,"pw.x")||contains(line,"projwfc.x")
        split_line = split(line)
        input_file = strip.(strip.(filter(x->contains(x,".in"),split_line),'<'),'>')[1]
        run_command = filter(x->contains(x,".x"),split_line)[1]
        push!(inputs,(run_command,input_file))
      end
    end
  end
  return inputs
end
#---------------------------END QUANTUM ESPRESSO SECTION----------------#
#---------------------------START WANNIER SECTION ----------------------#

#We also use DFInput to store wannier input, the non-named block in wannier is stored in the :control
#dictionary.
function read_wannier_input(filename::String, T=Float32)
  strip_split(line,args...) = strip.(split(line,args...))
  control_blocks = Dict{Symbol,Any}()
  control_blocks[:control] = Dict{Symbol,Any}()
  atoms = Dict{Symbol,Union{Symbol,Array{Point3D{T},1}}}()
  k_points = Dict{Symbol,Array{T,1}}()
  cell_param_dict = Dict{Symbol,Matrix{T}}()
  open(filename,"r") do f
    while !eof(f)
      line = readline(f)
      if contains(line,"!")||line==""
        continue
      end
      if contains(line,"Begin") || contains(line,"begin")
        block_name = Symbol(lowercase(split(line)[end]))
        line = readline(f)

        if block_name == :projections
          if line == "random"
            control_blocks[block_name] = Dict(:random=>Symbol[])
          else
            proj_dict = Dict{Symbol,Array{Symbol,1}}()
            while !(contains(line,"End") || contains(line,"end"))
              if contains(line,"!")
                line = readline(f)
                continue
              end
              split_line = strip_split(line,':')
              atom = Symbol(split_line[1])
              projections = [Symbol(proj) for proj in strip_split(split_line[2],';')]
              proj_dict[atom] = projections
              line = readline(f)
            end
            control_blocks[block_name] = proj_dict
          end
        elseif block_name == :kpoint_path
          k_path_array = Array{Tuple{Symbol,Array{T,1}},1}()
          while !(contains(line,"End") || contains(line,"end"))
            if contains(line,"!")
              line = readline(f)
              continue
            end
            split_line = split(line)
            push!(k_path_array,(Symbol(split_line[1]),parse_string_array(T,split_line[2:4])))
            push!(k_path_array,(Symbol(split_line[5]),parse_string_array(T,split_line[6:8])))
            line = readline(f)
          end
          control_blocks[block_name] = k_path_array
        elseif block_name == :unit_cell_cart
          if length(split(line)) == 1
            option = Symbol(lowercase(line))
          else
            option = :ang
          end
          cell_param = Matrix{T}(3,3)
          for i=1:3
            cell_param[i,:] = parse_line(T,readline(f))
          end
          cell_param_dict[option] = cell_param
          line = readline(f)

        elseif block_name ==:atoms_frac || block_name == :atoms_cart
          atoms[:option] = block_name
          while !(contains(line,"end") || contains(line,"End"))
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

        elseif block_name == :kpoints
          i=1
          while !(contains(line,"end") || contains(line,"End"))
            k_points[Symbol(i)] = parse_line(T,line)
            line = readline(f)
            i+=1
          end
        end

      else
        if contains(line,"mp_grid")
          control_blocks[:control][:mp_grid] = parse_string_array(Int,split(split(line,'=')[2]))
        else
          split_line = strip_split(line,'=')
          flag = Symbol(split_line[1])
          value = split_line[2]
          if  lowercase(value)=="t" || lowercase(value) == "true"
            control_blocks[:control][flag] = true
          elseif lowercase(value) == "f" || lowercase(value) == "false"
            control_blocks[:control][flag] = false
          elseif !isnull(tryparse(Int,value))
            control_blocks[:control][flag] = get(tryparse(Int,value))
          elseif !isnull(tryparse(T,value))
            control_blocks[:control][flag] = get(tryparse(T,value))
          else
            control_blocks[:control][flag] = value
          end
        end
      end
    end
  end
  return DFInput(:wannier,control_blocks,Dict{Symbol,String}(),cell_param_dict,atoms,k_points)
end


#---------------------------END WANNIER SECTION ------------------------#
#---------------------------BEGINNING GENERAL SECTION-------------------#
"""
    write_df_input(filename::String,df_input::DFInput)

Writes the input file for a DFInput.
Backend of DFInput decides what writing function is called.
"""
function write_df_input(filename::String,df_input::DFInput)
  if df_input.backend == :QE
    write_qe_input(filename,df_input)
  else
    error("Backend $(df_input.backend) is not yet implemented!")
  end
end

#@TODO after we added defaults and extra job settings to the DFJob type add those instead of hardcoded values
#Possible issue with .in and .win .. maybe we should use the name that is presented when making the job
#@TODO actually having a flow chart wouldn't be a bad idea, now it's always the same order!
#@Incomplete there might be the case that we want multiple different nscf or scf or whatever inputs. Maybe using array of tuple(string,dfinput)
# would be more suitable
"""
    write_job_files(df_job::DFJob)

Writes all the input files and job file that are linked to a DFJob.
"""
function write_job_files(df_job::DFJob)
  files_to_remove = search_dir(df_job.home_dir,".in")
  new_filenames   = String[]

  open(df_job.home_dir*"job.tt","w") do f
    write(f,"#!/bin/bash\n","#SBATCH -N 1\n","#SBATCH --ntasks-per-node=24 \n","#SBATCH --time=24:00:00 \n","#SBATCH -J $(df_job.job_name) \n",
          "#SBATCH -p defpart\n\n","module load open-mpi/gcc/1.10.2\n","module load mkl/2016.1.056\n","\n")
    for (i,(run_command,input)) in enumerate(df_job.flow)
      filename = "$i"*df_job.job_name*"_$input"*".in"
      push!(new_filenames,filename)
      write_df_input(df_job.home_dir*filename,df_job.calculations[input])
      write(f,"mpirun -np 24 $run_command <$filename> $(split(filename,".")[1]).out \n")
    end
  end

  for file in files_to_remove
    if !in(file,new_filenames)
      rm(df_job.home_dir*file)
    end
  end
end

#---------------------------END GENERAL SECTION-------------------#
