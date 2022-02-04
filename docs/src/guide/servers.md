# [Servers](@id servers_header)
`DFControl` is structure as a client-server architecture where the communication happens through a rest-API. This means that
on the server-side a small daemon will run that the client will communicate with in order to `load`, `save` and `submit` [`Jobs`](@ref Job).
A [`Server`](@ref) can be setup simply by calling the constructor with a `String` that signifies the name of the [`Server`](@ref), e.g. `"loalhost"`, or an `ssh` string such as `"user@domain"`, which will prompt an interactive setup menu.

Alternatively a [`Server`](@ref) can be set up by filling out the constructor, and saving it as `save(server)`.
A previously saved [`Server`](@ref) can be loaded again through e.g. `Server("localhost")` which will retrieve the previously saved
configuration.
```@docs
Server
start
save(::Server)
```
