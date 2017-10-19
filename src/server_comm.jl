
function read_errors(server::String,server_dir::String; error_fuzzies=["CRASH","*.werr"])
  server_dir = form_directory(server_dir)
  tmp_dir = joinpath(@__DIR__,"tmp")
  if !isdir(tmp_dir)
    mkdir(tmp_dir)
  end
  for fuzzy in error_fuzzies
    run(`scp $(server*":"*server_dir*fuzzy) $local_dir`)
  end

  #for now very dumb!
  crash_readlines = Dict{Symbol,Array{String,1}}()
  for fuzzy in fuzzies
    filenames = search_dir(tmp_dir,strip(fuzzy,'*'))
    if length(filenames)==1
      crash_readlines[filename[1]] = readlines(filenames[1])
    end
  end
  return crash_readlines
end
