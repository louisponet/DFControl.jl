using ..Database
function Database.load(req::HTTP.Request)
    p = path(req)
    if !isempty(HTTP.body(req))
        # This part basically exists for when more complicated things need to be done
        # when loading an entity (e.g. a Job).
        typ = Symbol(HTTP.header(req, "Type"))
        val = eval(:(JSON3.read($(req.body), $typ)))
        if Database.exists(val)
            return Database.load(val)
        else
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
HTTP.@register(ROUTER, "GET", "/database/*", load)

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
HTTP.@register(ROUTER, "POST", "/database/*", save)

function database_rm(req)
    ispath(path(req))
    rm(path(req))
end
HTTP.@register(ROUTER, "PUT", "/database/*", database_rm)
 
