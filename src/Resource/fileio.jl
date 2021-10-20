HTTP.@register(ROUTER, "GET", "/fileio/qe/parse/output/pw", req -> FileIO.qe_read_pw_output(JSON3.read(req.body, String)))
HTTP.@register(ROUTER, "GET", "/fileio/qe/parse/output/projwfc", req -> FileIO.qe_read_projwfc_output(JSON3.read(req.body, String)))
HTTP.@register(ROUTER, "GET", "/fileio/qe/parse/output/hp", req -> FileIO.qe_read_hp_output(JSON3.read(req.body, String)))
HTTP.@register(ROUTER, "GET", "/fileio/qe/parse/input", req -> FileIO.qe_read_calculation(JSON3.read(req.body, String)))

HTTP.@register(ROUTER, "GET", "/fileio/wannier90/parse/input", req -> FileIO.wan_read_calculation(JSON3.read(req.body, String)))
HTTP.@register(ROUTER, "GET", "/fileio/wannier90/parse/output", req -> FileIO.wan_read_output(JSON3.read(req.body, String)))
