using Libdl, FLINT_jll, Pkg.Artifacts

old_dir = pwd()
deps_dir = dirname(@__FILE__)
flint_jll_dir = dirname(dirname(pathof(FLINT_jll)))
flint_hash = artifact_hash("FLINT", joinpath(flint_jll_dir, "Artifacts.toml"))
flint_lib = joinpath(artifact_path(flint_hash), "lib")
flint_headers = joinpath(artifact_path(flint_hash), "include")

# Make libembed

cd(deps_dir)
run(`make libembed.so ARG1=$flint_headers ARG2=$flint_lib`)

# Add the path

push!(Libdl.DL_LOAD_PATH, deps_dir)
cd(old_dir)
