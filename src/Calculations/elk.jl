include(joinpath(DFC.DEPS_DIR, "elkflags.jl"))

struct ElkFlagInfo{T}
    name::Symbol
    default::Union{T,Nothing}  #check again v0.7 Some
    description::String
end
ElkFlagInfo() = ElkFlagInfo{Nothing}(:error, "")
Base.eltype(x::ElkFlagInfo{T}) where {T} = T

struct ElkControlBlockInfo <: AbstractBlockInfo
    name::Symbol
    flags::Vector{<:ElkFlagInfo}
    description::String
end

const ELK_CONTROLBLOCKS = _ELK_CONTROLBLOCKS()

function elk_flaginfo(flag::Symbol)
    for b in ELK_CONTROLBLOCKS
        for f in b.flags
            if f.name == flag
                return f
            end
        end
    end
end

elk_block_info(name::Symbol) = getfirst(x -> x.name == name, ELK_CONTROLBLOCKS)

function elk_block_variable(flag_name::Symbol)
    for b in ELK_CONTROLBLOCKS
        for f in b.flags
            if f.name == flag_name
                return b
            end
        end
    end
end

flagtype(::Calculation{Elk}, flag::Symbol) = eltype(elk_flaginfo(flag))

infilename(c::Calculation{Elk}) = "elk.in"
isbands(c::Calculation{Elk}) = c.name == "20"
isnscf(c::Calculation{Elk}) = c.name == "elk2wannier" #nscf == elk2wan??
isscf(c::Calculation{Elk}) = c.name ∈ ["0", "1"]
