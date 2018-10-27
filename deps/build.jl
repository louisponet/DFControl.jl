using Compat
using BinDeps

if !isdir(joinpath(@__DIR__,"../user_defaults"))
    mkdir(joinpath(@__DIR__,"../user_defaults"))
end
if !ispath(joinpath(@__DIR__,"../user_defaults/user_defaults.jl"))
    open(joinpath(@__DIR__,"../user_defaults/user_defaults.jl"), "w") do f
        write(f,"#Generated by DFControl.")
    end
end



relpath = x -> joinpath(@__DIR__, x)
dlpath = relpath("downloads")
if !ispath(dlpath)
    mkdir(dlpath)
end
pythonpath = relpath("python2")
if !ispath(pythonpath)
    mkdir(pythonpath)
end
#From Conda installation
url = "https://repo.continuum.io/miniconda/Miniconda2-latest-"
if Compat.Sys.isapple()
    url *= "MacOSX"
elseif Compat.Sys.islinux()
    url *= "Linux"
elseif Compat.Sys.iswindows()
    url = "https://repo.continuum.io/miniconda/Miniconda2-4.5.4-"
    url *= "Windows"
else
    error("Unsuported OS.")
end
url *= Sys.WORD_SIZE == 64 ? "-x86_64" : "-x86"
url *= Compat.Sys.iswindows() ? ".exe" : ".sh"

if Compat.Sys.isunix()
    installer = joinpath(dlpath, "installer.sh")
end
if Compat.Sys.iswindows()
    installer = joinpath(dlpath, "installer.exe")
end
download(url, installer)
if Compat.Sys.isunix()
    chmod(installer, 33261)  # 33261 corresponds to 755 mode of the 'chmod' program
    run(`$installer -b -f -p $pythonpath`)
end
if Compat.Sys.iswindows()
    run(Cmd(`$installer /S /AddToPath=0 /RegisterPython=0 /D=$pythonpath`, windows_verbatim=true))
end

tarpath = relpath("cif2cell.tar.gz")
download("https://sourceforge.net/projects/cif2cell/files/latest/download", tarpath)
run(unpack_cmd("cif2cell.tar.gz", @__DIR__, ".gz",".tar"))
cif2celldir = relpath("cif2cell-1.2.10")
cd(cif2celldir)
pyex = Compat.Sys.iswindows() ? joinpath(pythonpath, "python") : joinpath(pythonpath, "bin", "python2")
run(`$pyex setup.py install --prefix=$pythonpath`)
cd("..")

rm(tarpath)
rm(relpath("cif2cell-1.2.10"), recursive=true)
rm(relpath("downloads"), recursive=true)
