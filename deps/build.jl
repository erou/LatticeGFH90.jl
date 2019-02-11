using Libdl, Nemo

old_dir = pwd()
deps_dir = dirname(@__FILE__)
nemo_dir = dirname(dirname(pathof(Nemo)))
flint_dir = joinpath(nemo_dir, "local", "lib")

# Make libembed

cd(deps_dir)
run(`make libembed.so ARGS=$flint_dir`)

# Add the path

push!(Libdl.DL_LOAD_PATH, deps_dir)
cd(old_dir)
