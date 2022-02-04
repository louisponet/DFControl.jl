module Database
    using ..DFControl: config_path
    using JSON3
    using HTTP

    abstract type Storable end

    storage_name(s::S) where {S<:Storable} = hasfield(S, :name) ? s.name : error("Please define a name function for type $S.")
    storage_directory(s::S) where {S<:Storable} = error("Please define a storage_directory function for type $S.")
    
    storage_path(s::S) where {S<:Storable} = joinpath("storage", storage_directory(s), storage_name(s) * ".json")
    storage_url(s::Storable) = joinpath("/database", storage_path(s))
    verify(s::Storable) = nothing

    # These are the standard functions where things are saved as simple jsons in the
    # config_path.
    """
        save([server::Server], e::Environment)
        save([server::Server], e::Exec)
        save([server::Server], s::Server)

    Saves an item to the database of `server`. If `server` is not specified the item will be stored in the local database.
    """
    function save(s::Storable)
        p = config_path(storage_path(s))
        mkpath(splitdir(p)[1])
        if ispath(p)
            @warn "Overwriting previously existing item at $p."
        end
        verify(s)
        JSON3.write(p, s)
    end
    function save(server, s::Storable)
        url = storage_url(s) 
        HTTP.post(server, url, s)
    end

    """
        load([server::Server], e::Environment)
        load([server::Server], e::Exec)
        load([server::Server], s::Server)

    Loads a previously stored item from the `server`. If `server` is not specified the local item is loaded.
    """
    function load(s::S) where {S <: Storable}
        p = config_path(storage_path(s))
        if !ispath(p)
            error("No item found at $p.")
        end
        JSON3.read(read(p, String), S)
    end
    function load(server, s::S) where {S <: Storable}
        url = storage_url(s)
        if isempty(s.name)
            if s == S() # asking for all possibilities
                return JSON3.read(HTTP.get(server, splitdir(url)[1]).body, Vector{String})
            else
                n = name(server, s) 
                if n !== nothing # stored item with matching data existed, just replace name
                    s.name = n
                    return s
                else
                    @warn "No exact match found. Returning the closest options..."
                    return JSON3.read(HTTP.get(server, url, s).body, Vector{String})
                end
            end
        elseif exists(server, s) # asking for a stored item
            return JSON3.read(JSON3.read(HTTP.get(server, url).body, String), S)
                
        else
            @warn "No exact match found. Returning the closest options..."
            return JSON3.read(HTTP.get(server, url, s).body, Vector{String})
        end
    end

    exists(s::Storable) = ispath(config_path(storage_path(s)))
    function exists(server, s::Storable)
        url = storage_url(s)
        possibilities = JSON3.read(HTTP.get(server, splitdir(url)[1]).body, Vector{String})
        return storage_name(s) ∈ possibilities
    end

    function replacements(s::S) where {S<:Storable}
        dir = config_path("storage", storage_directory(s))
        default = S() # We only score matches that are not the default
        all = map(x -> JSON3.read(read(joinpath(dir, x), String), S), readdir(dir))
        if isempty(all)
            return S[]
        end
        score = zeros(Int, length(all))
        for (i, t) in enumerate(all)
            for f in fieldnames(S)
                if getfield(t, f) == getfield(s, f) != getfield(default, f)
                    score[i] += 1
                end
            end
        end
        s = sortperm(score, rev=true)
        best = maximum(score)
        return all[s][1:length(findall(isequal(best), score))]
    end

    function name(s::S) where {S<:Storable}
        dir = config_path("storage", storage_directory(s))
        for f in readdir(dir)
            t =  JSON3.read(read(joinpath(dir, f), String), S)
            if t == s
                return t.name
            end
        end
    end
    function name(server, s::Storable)
        resp = HTTP.get(server, "/database/name", s)
        if resp.status == 204
            return nothing
        else
            return JSON3.read(resp.body, String)
        end
    end

    """
        Base.rm([server::Server], e::Environment)
        Base.rm([server::Server], e::Exec)
        Base.rm([server::Server], s::Server)

    Removes an item from the database of `server`. If `server` is not specified the item will removed from the local database.
    """
    function Base.rm(s::Storable)
        p = config_path(storage_path(s))
        if !ispath(p)
            error("No item found at $p.")
        end
        rm(p)
    end
        
    function Base.rm(server, s::Storable)
        url = storage_url(s) 
        HTTP.put(server, url)
    end

    
        
    export Storable
    export save, load, exists
   
    
end
