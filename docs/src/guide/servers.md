# [Servers](@id servers_header)
`DFControl` is structure as a client-server architecture where the communication happens through a rest-API. This means that
on the server-side a small daemon will run that the client will communicate with in order to `load`, `save` and `submit` [`Jobs`](@ref Job).
There are three ways to set up a server. To setup a local server, simply call [`configure_local()`](@ref) and follow the prompt. A [`Server`](@ref) can also be setup by calling the constructor with a `String` that signifies the name of the [`Server`](@ref), e.g. `"loalhost"`, or an `ssh` string such as `"user@domain"`, which will prompt an interactive setup menu. 

Alternatively, for full cusmotization, a [`Server`](@ref) can be set up by filling out the constructor, and saving it as `save(server)`.
A previously saved [`Server`](@ref) can be loaded again through e.g. `Server("localhost")` which will retrieve the previously saved
configuration.
```@docs
Server
configure_local
start
save(::Server)
```
