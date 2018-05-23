#came from input.jl
"""
    print_flag(input::DFInput, flag)

Prints information of the flag inside the input.
"""
function print_flag(input::DFInput, flag)
    if (:control in fieldnames(input))
        for block in input.control
            if haskey(block.flags, flag)
                s = """
                Filename: $(input.filename)
                Block Name: $(block.name)
                $flag => $(block.flags[flag])\n
                """
                dfprintln(s)
            end
        end
    end

    if (:flags in fieldnames(input))
        if haskey(input.flags, flag)
            s = """
            Filename: $(input.filename)
            $flag => $(input.flags[flag])\n
            """
            dfprintln(s)
        end
    end
end

print_flags(input::DFInput, flags::Array) = print_flag.(input, flags)

"""
    print_block(input::DFInput, block_name::Symbol)

Print the information of a 'Block' inside the input.
"""
function print_block(input::DFInput, block_name::Symbol)
    found        = false
    input_blocks = input.data
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
    input.data |> display
    dfprintln("")
end

"""
    print_filename(input::DFInput)

Prints the filename associated with the input.
"""
function print_filename(input::DFInput)
    dfprintln("Input file: $(input.filename)")
end

"""
    print_info(input::DFInput)

Prints general info of the input.
"""
function print_info(input::DFInput)
    dfprintln("Filename: $(input.filename)")
    if (:control in fieldnames(input))
        dfprintln("  Control Blocks:")
        for (i, block) in enumerate(input.control)
            dfprintln("    $i: $(block.name)")
        end
    end
    dfprintln("  Data Blocks:")
    for (i, block) in enumerate(input.data)
        dfprintln("    $i: $(block.name)")
    end
    dfprintln("  Run command: $(input.runcommand)")
    dfprintln("  Runs: $(input.run)")
end

"""
    print_flags(input::DFInput)

Prints all the flags of the input.
"""
function print_flags(input::DFInput)
    dfprintln("#----------------#")
    dfprintln("Filename: $(input.filename)")
    for (flag, value) in input.flags
        dfprintln("  $flag => $value")
    end
    dfprintln("#----------------#\n")
end

#---------------------------------------------------------
#came from job_control
"""
    print_run_command(job::DFJob, filenames)

Prints the run command of the specified calculations.
"""
function print_run_command(job::DFJob, filenames)
    for calc in inputs(job, filenames)
        dfprintln("Run command of file '$(calc.filename)' is: '$(calc.runcommand)'.")
        dfprintln("")
    end
end

"""
    print_flow(job::DFJob)

Prints the calculation sequence of the job.
"""
function print_flow(job::DFJob)
    for (i, calc) in enumerate(job.calculations)
        dfprintln("$i: $(calc.filename) => runs: $(calc.run)")
    end
end

"""
    print_block(job::DFJob, block_name::Symbol)

Prints information of the specified block name of all the calculations in the job.
"""
function print_block(job::DFJob, block_name::Symbol)
    for calc in job.calculations
        if print_block(calc, block_name) dfprintln("") end
    end
end

"""
    print_block(job::DFJob, filenames, name::Symbol)

Prints the information of the block in a selected file of the job.
"""
function print_block(job::DFJob, filenames, name::Symbol)
    for calc in inputs(job, filenames)
        print_block(calc, name)
    end
end

"""
    print_blocks(job::DFJob, calc_filenames)

Prints information on all the data in the specified calculations.
"""
function print_blocks(job::DFJob, calc_filenames)
    for calc in inputs(job, calc_filenames)
        print_blocks(calc)
    end
end

"""
    print_blocks(job::DFJob)

Prints information of all the data of all the calculations in the job.
"""
function print_blocks(job::DFJob)
    for calc in job.calculations
        print_blocks(calc)
        dfprintln("#------------------------------------#")
    end
end
print_data(job::DFJob)                     = print_blocks(job)
print_data(job::DFJob, calc_filenames)     = print_blocks(job, calc_filenames)
print_data(job::DFJob, block_name::Symbol) = print_block(job, block_name)

"""
    print_info(job::DFJob, filenames::Array{String,1})

Prints general info of the job, and the specified filenames.
"""
function print_info(job::DFJob, filenames::Array{String,1})
    s = """--------------------
    DFJob:      $(job.name)
    Local_dir:  $(job.local_dir)
    Server:     $(job.server)
    Server_dir: $(job.server_dir)
    $(length(job.calculations)) calculations
    --------------------
    """
    dfprintln(s)

    for calc in inputs(job, filenames)
        print_info(calc)
        dfprintln(" ")
    end

end
print_info(job::DFJob)                   = print_info(job, [calc.filename for calc in job.calculations])
print_info(job::DFJob, filename::String) = print_info(job, [filename])

"""
    print_flags(job::DFJob)

Prints flags of all the calculations in the job.
"""
function print_flags(job::DFJob)
    for calc in job.calculations
        print_flags(calc)
    end
end

"""
    print_flags(job::DFJob, calc_filename::String)

Prints flags of the specified calculation.
"""
function print_flags(job::DFJob, calc_filename::String)
    for calc in inputs(job, calc_filename)
        print_flags(calc)
    end
end

"""
    print_flags(job::DFJob, calc_filenames::Array{String,1})

Prints the flags of the specified calculations.
"""
function print_flags(job::DFJob, calc_filenames::Array{String,1})
    for file in calc_filenames
        print_flags(job, file)
    end
end

"""
    print_flags(job::DFJob, flags)

Prints the specified flags running through all the claculations in the job.
"""
function print_flags(job::DFJob, flags::Array{Symbol,1})
    for flag in flags
        print_flag(job, flag)
    end
end

"""
    print_flag(job::DFJob, flag::Symbol)

Prints the specified flag running through all the calculations in the job.
"""
function print_flag(job::DFJob, flag::Symbol)
    for calc in job.calculations
        print_flag(calc, flag)
    end
end

#--------------------------------------
