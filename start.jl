# activate the environment
using Pkg
Pkg.activate(joinpath(@__DIR__))

# start Pluto
using Pluto
Pluto.run()