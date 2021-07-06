const ORBITALS = [Orbital("s", 1, 0, 1), Orbital("p", 3, 1, 0), Orbital("px", 1, 1, 2),
                  Orbital("py", 1, 1, 3), Orbital("pz", 1, 1, 1), Orbital("d", 5, 2, 0),
                  Orbital("dz2", 1, 2, 1), Orbital("dxz", 1, 2, 2), Orbital("dyz", 1, 2, 3),
                  Orbital("dx2-y2", 1, 2, 4), Orbital("dxy", 1, 2, 5),
                  Orbital("f", 7, 3, 0), Orbital("fz3", 1, 3, 1), Orbital("fxz2", 1, 3, 2),
                  Orbital("fyz2", 1, 3, 3), Orbital("fz(x2-y2)", 1, 3, 4),
                  Orbital("fxyz", 1, 3, 5), Orbital("fx(x2-3y2)", 1, 3, 6),
                  Orbital("fy(3x2-y2)", 1, 3, 7), Orbital("sp", 2, -1, 0),
                  Orbital("sp2", 3, -2, 0), Orbital("sp3", 4, -3, 0),
                  Orbital("sp3d2", 6, -5, 0)]
"""
    orbital(s::String)

Returns the [Orbital](@ref Orbitals) identified by `s`. 
"""
orbital(s::String) = getfirst(x -> x.name == s, ORBITALS)
orbital(l::Number) = getfirst(x -> x.l == l, ORBITALS)

Base.convert(::Type{String}, x::Orbital) = x.name
orbsize(orb::Orbital) = orb.size
orbsize(orb::String) = orbsize(orbital(orb))

orbital(proj::Projection) = proj.orb
orbsize(proj::Projection) = proj.last - proj.start + 1

function Base.show(io::IO, proj::Projection)
    dfprintln(io, crayon"cyan", "Orbital: ", crayon"reset", "$(proj.orb.name)")
    dfprintln(io, crayon"red", "start index: ", crayon"reset", "$(proj.start)")
    return dfprintln(io, crayon"red", "last index: ", crayon"reset", "$(proj.last)")
end

Base.hash(orb::Orbital, h::UInt) = hash(orb.mr, hash(orb.l, h))
Base.:(==)(o1::Orbital, o2::Orbital) = o1.l == o2.l && o1.mr == o2.mr

"""
Adds projections to atoms.
"""
function addprojections!(atoms, projections_, soc; print = true)
    t_start = 1
    for (proj_at, projs) in projections_
        for proj in projs
            for at in atoms
                if length(projs) <= length(projections(at))
                    continue
                end
                if name(at) == proj_at
                    size   = soc ? 2 * orbsize(proj) : orbsize(proj)
                    t_proj = Projection(orbital(proj), t_start, t_start + size - 1)
                    push!(projections(at), t_proj)
                    t_start += size
                end
            end
        end
    end
    print && @info "Projections are now:"
    for at in atoms
        if !isempty(projections(at))
            print && @info "$(name(at));$(position_cryst(at)) :"
            for p in projections(at)
                print && Base.show(p)
            end
        end
    end
end

Base.zero(::Type{Projection}) = Projection(Orbital("zero", 0, 0, 0), 0, 0)
Base.zero(::Projection) = Projection(Orbital("zero", 0, 0, 0), 0, 0)

Base.range(proj::Projection) = proj.start:proj.last

"Returns all the projections inside the job."
projections(job::DFJob) = projections(structure(job))

"""
    set_projections!(at::AbstractAtom, projections::Vector{Projection}; print=true)

Sets the projections of the atom.

    set_projections!(str::Structure, projs::Pair...; soc=false)
    set_projections!(job::DFJob, projs::Pair...; kwargs...)

Sets the projections of the specified atoms. `projs` has to be of form `:atsym => [:proj]`,
where `proj = "s", "p", "d", "f"`, see [Orbital](@ref Orbitals) to see which projections are allowed. 
If `soc` is set to `true` both up and down projections will be taken into account.
"""
function set_projections!(at::AbstractAtom, projections::Vector{Projection}; print = true)
    print && @info "Setting projections for atom $(name(at)) to $projections"
    return atom(at).projections = projections
end

function set_projections!(job::DFJob, projections...; kwargs...)
    socid = findfirst(issoc, calculations(job))
    return set_projections!(job.structure, projections...; soc = socid !== nothing,
                            kwargs...)
end

function set_projections!(str::Structure, projs::Pair...; soc = false, kwargs...)
    projdict = Dict{Symbol, typeof(projs[1][2])}()
    for (sym, proj) in projs
        ats = atoms(str, sym)
        for a in ats
            projdict[name(a)] = proj
        end
    end
    for at in unique(str.atoms)
        if !haskey(projdict, name(at))
            projdict[name(at)] = [proj.orb for proj in projections(at)]
        end
    end
    empty_projections!(str)
    return addprojections!(atoms(str), projdict, soc; kwargs...)
end

function empty_projections!(str::Structure)
    for at in str.atoms
        empty!(projections(at))
    end
end

function nprojections(structure)
    n = 0
    for at in atoms(structure)
        projs = projections(at)
        if !isempty(projs)
            n += sum(orbsize.(projs))
        end
    end
    return n
end
