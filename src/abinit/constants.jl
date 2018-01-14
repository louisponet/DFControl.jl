const ABI_UNIT_NAMES = [lowercase(s) for s in [
        "au",
        "Angstr", "Angstrom", "Angstroms", "Bohr", "Bohrs",
        "eV", "Ha", "Hartree", "Hartrees", "K", "Ry", "Rydberg", "Rydbergs",
        "T", "Tesla"]]

const abi_conversions = OrderedDict{Symbol,Any}(:ev => 1 / 27.2113845,
                                         :ha => 1.0, 
                                         :ry => 0.5, 
                                         :ang => 1.889716164632,
                                         :bohr => 1.0,
                                         :au => 1.0,
                                         :angstr => 1.889716164632,
                                         :angstrom => 1.889716164632,
                                         :bohrs => 1.0,
                                         :hartree => 1.0,
                                         :hartrees => 1.0,
                                         :k => 1/3.1577513e5,
                                         :rydberg => 0.5,
                                         :rydbergs => 0.5,
                                         :t =>1/2.35e5,
                                         :tesla => 1/2.35e5)

convert_2abi(value, s::String) = value * abi_conversions[Symbol(s)]
convert_2abi(value, s::Symbol) = value * abi_conversions[s]

function construct_abi_flags()
    out      = OrderedDict{Symbol,Type}()
    open(joinpath(assets_dir,"inputs/abinit/input_variables.txt"), "r") do f
        while !eof(f)
            spl = split(readline(f))
            out[parse(spl[1])] = eval(parse(spl[2]))
        end
    end
    return out
end

const AbinitFlags    = construct_abi_flags()

get_abi_flag_type(flag) = haskey(AbinitFlags, flag) ? AbinitFlags[flag] : Void