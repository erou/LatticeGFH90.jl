using Libdl, Nemo

old_dir = pwd()
deps_dir = dirname(@__FILE__)
nemo_dir = dirname(dirname(pathof(Nemo)))
flint_lib = joinpath(nemo_dir, "deps", "usr", "lib")
flint_headers = joinpath(nemo_dir, "deps", "usr", "include")

# Make libembed

cd(deps_dir)
run(`make libembed.so ARG1=$flint_headers ARG2=$flint_lib`)

# Add the path

push!(Libdl.DL_LOAD_PATH, deps_dir)
cd(old_dir)
