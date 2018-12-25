using DFControl: depsdir
const ABI_UNIT_NAMES = [lowercase(s) for s in [
        "au",
        "Angstr", "Angstrom", "Angstroms", "Bohr", "Bohrs",
        "eV", "Ha", "Hartree", "Hartrees", "K", "Ry", "Rydberg", "Rydbergs",
        "T", "Tesla"]]

#These are in String format because we throw them out even before
#parsing them.
const ABI_STRUCTURE_FLAGS = ["acell", "angdeg", "rprim",
                            "ntypat", "natom", "znucl",
                            "typat", "xred", "xcart", "xangst"]

#These are in Symbol format because we use them after parsing.
const ABI_KPOINTS_FLAGS   = [:kptopt, :ngkpt, :kpt, :kptbounds,
                             :kptnrm, :kptrlatt, :kptrlen, :ndivk,
                             :nkpath, :nkpt, :nshiftk, :wtk]

const ABI_KPTOPT_TO_OPTION = Dict( 0 => :direct,
                                   1 => :auto_fullsymm,
                                   2 => :auto_TRsymm,
                                   3 => :auto_nosymm,
                                   4 => :auto_lattsymm,
                                  -1 => :bands)

abi_kptopt_to_option(i::Int) =
    sign(i) == -1 ? ABI_KPTOPT_TO_OPTION[-1] : ABI_KPTOPT_TO_OPTION[i]
    
#TODO: Handle required flags really
function abi_required_flags(s::Symbol)
    if s == :direct
        return [:nkpt, :kpt, :kptnrm, :wtk]
    elseif s == :bands
        return [:kptbounds, :ndivk]
    else
        return [:ngkpt, :kptrlatt, :nshiftk, :shiftk]
    end
end

#REVIEW: I don't really like this, I feel like it should be all eV etc, but ok fine.
#        These are the natural units of abinit.
const abi_conversions = Dict{Symbol,Any}(:ev       => 1 / 27.2113845,
                                         :ha       => 1.0,
                                         :ry       => 0.5,
                                         :ang      => 1.889716164632,
                                         :bohr     => 1.0,
                                         :au       => 1.0,
                                         :angstr   => 1.889716164632,
                                         :angstrom => 1.889716164632,
                                         :bohrs    => 1.0,
                                         :hartree  => 1.0,
                                         :hartrees => 1.0,
                                         :k        => 1/3.1577513e5,
                                         :rydberg  => 0.5,
                                         :rydbergs => 0.5,
                                         :t        => 1/2.35e5,
                                         :tesla    => 1/2.35e5)

convert2abi(value, s::String) = value * abi_conversions[Symbol(s)]
convert2abi(value, s::Symbol) = value * abi_conversions[s]

include(joinpath(depsdir, "abinitflags.jl"))
const AbinitFlags  = _ABINITFLAGS()
#
flagtype(::Type{ABINIT}, flag) = haskey(AbinitFlags, flag) ? AbinitFlags[flag] : Nothing
