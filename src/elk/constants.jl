include(joinpath(depsdir, "elkflags.jl"))

struct ElkFlagInfo{T}
    name::Symbol
    default::Union{T, Nothing}  #check again v0.7 Some
    description::String
end
ElkFlagInfo() = ElkFlagInfo{Nothing}(:error, "")
Base.eltype(x::ElkFlagInfo{T}) where T = T


struct ElkControlBlockInfo <: AbstractBlockInfo
    name::Symbol
    flags::Vector{<:ElkFlagInfo}
    description::String
end

const ELK_CONTROLBLOCKS = _ELK_CONTROLBLOCKS()

const ELK_EXECS = ["elk", "elk-omp"]


elk_block_info(name::Symbol) = getfirst(x->x.name == name, ELK_CONTROLBLOCKS)

function elk_block_variable(flag_name::Symbol)
	for b in ELK_CONTROLBLOCKS
		for f in b.flags
			if f.name == flag_name
				return b
			end
		end
	end
end

