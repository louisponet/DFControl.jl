using ..Database
function Database.load(req::HTTP.Request)
    p = path(req)
    if !isempty(HTTP.body(req))
        # This part basically exists for when more complicated things need to be done
        # when loading an entity (e.g. a Job).
        typ = Symbol(HTTP.header(req, "Type"))
        val = eval(:(JSON3.read($(req.body), $typ)))
        try
            return Database.load(val)
        catch
            return map(x->Database.storage_name(x), Database.replacements(val))
        end
    else
        cpath = config_path(p) 
        if isempty(splitext(p)[end])
            # Here we return the possibilities
            return map(x->splitext(x)[1], readdir(cpath))
        else
            return read(cpath, String)
        end
    end
end
HTTP.register!(ROUTER, "GET", "/database/storage/**", load)

function Database.save(req::HTTP.Request)
    p = path(req)
    if HTTP.hasheader(req, "Type")
        # This part basically exists for when more complicated things need to be done
        # when storing an entity (e.g. a Job).
        typ = Symbol(HTTP.header(req, "Type"))
        eval(:(save(JSON3.read($(req.body), $typ))))
    else
        mkpath(splitdir(p)[1])
        write(p, req.body)
    end
end
HTTP.register!(ROUTER, "POST", "/database/storage/**", save)

function database_rm(req)
    p = config_path(path(req))
    ispath(p)
    rm(p)
end
HTTP.register!(ROUTER, "PUT", "/database/storage/**", database_rm)

function Database.name(req)
    typ = Symbol(HTTP.header(req, "Type"))
    val = eval(:(JSON3.read($(req.body), $typ)))
    return Database.name(val)
end
    
HTTP.register!(ROUTER, "GET", "/database/name", Database.name)
 
