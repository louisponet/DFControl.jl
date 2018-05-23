"""
    print_info(job::DFJob, filenames::Vector{String})

Prints general info of the job.
"""
function print_info(job::DFJob, filenames::Vector{String})
    s = """--------------------
    DFJob:      $(job.name)
    Local_dir:  $(job.local_dir)
    Server:     $(job.server)
    Server_dir: $(job.server_dir)
    $(length(job.calculations)) calculations
    --------------------
    """
    dfprintln(s)
end
print_info(job::DFJob)                   = print_info(job, [calc.filename for calc in job.calculations])
print_info(job::DFJob, filename::String) = print_info(job, [filename])

"""
Prints input flags.
"""
function print_flags(input::DFInput)
    dfprintln("#----------------#")
    dfprintln("Filename: $(input.filename)")
    for (flag, value) in input.flags
        dfprintln("  $flag => $value")
    end
    dfprintln("#----------------#\n")
end
print_flags(job::DFJob) = print_flags.(job.calculations)
print_flags(job::DFJob, filename::String) = print_flags.(inputs(job, filename))
print_flag(job::DFJob, flag::Symbol) = print_flag.(job.calculations, flag)

#--------------------------------------
