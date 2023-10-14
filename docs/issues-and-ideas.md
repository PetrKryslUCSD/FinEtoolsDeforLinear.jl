Issues and ideas:


-- Documenter:
using FinEtoolsDeforLinear
using DocumenterTools
Travis.genkeys(user="PetrKryslUCSD", repo="https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl")
DocumenterTools.genkeys(user="PetrKryslUCSD", repo="git@github.com:PetrKryslUCSD/FinEtoolsDeforLinear.jl.git")

– Get rid of lumpedmass(): It is already implemented with an HRZ assembler.

-- Formatting: JuliaFormatter
using JuliaFormatter
format("./examples", SciMLStyle())   
