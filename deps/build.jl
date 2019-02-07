using Libdl

old_dir = pwd()
deps_dir = dirname(@__FILE__)
pkg_dir = dirname(deps_dir)
cd(pkg_dir)

# Cloning the repo

run(`git clone https://github.com/erou/gf-h90-lattice.git`)
mv(joinpath(pkg_dir, "gf-h90-lattice/implementation"), joinpath(pkg_dir,
                                                                "deps/implementation"))
rm(joinpath(pkg_dir, "gf-h90-lattice"), force = true, recursive = true)

# Make libembed

cd(joinpath(deps_dir, "implementation"))
run(`make libembed.so`)

# Add the path

push!(Libdl.DL_LOAD_PATH, deps_dir)
cd(old_dir)
