using Libdl

old_dir = pwd()
deps_dir = dirname(@__FILE__)
pkg_dir = dirname(deps_dir)
cd(pkg_dir)

# Cloning the repo

run(`git clone https://github.com/erou/gf-h90-lattice.git`)
cd(joinpath(pkg_dir, "gf-h90-lattice/implementation"))
run(`mv * ../../deps`)
cd(pkg_dir)
rm(joinpath(pkg_dir, "gf-h90-lattice"), force = true, recursive = true)

# Make libembed

cd(deps_dir)
run(`make libembed.so`)

# Add the path

push!(Libdl.DL_LOAD_PATH, deps_dir)
cd(old_dir)
