using Pkg
Pkg.activate(".")
Pkg.instantiate()
using FinEtoolsDeforLinear

# Pkg.test()
include("examples/statics/3-d/bend_hex_spectrum_examples.jl"); bend_hex_spectrum_examples.allrun()
include("examples/statics/3-d/comp_hex_spectrum_examples.jl"); comp_hex_spectrum_examples.allrun()

