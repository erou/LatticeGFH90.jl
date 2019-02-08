using Libdl

old_dir = pwd()
deps_dir = dirname(@__FILE__)

# Make libembed

cd(deps_dir)
run(`make libembed.so`)

# Add the path

push!(Libdl.DL_LOAD_PATH, deps_dir)
cd(old_dir)
