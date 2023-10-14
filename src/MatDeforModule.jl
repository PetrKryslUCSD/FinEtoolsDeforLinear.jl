"""
    MatDeforModule

Module to support general operations for deformation material models.
"""
module MatDeforModule

__precompile__(true)

using FinEtools.FTypesModule:
    FInt,
    FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.DeforModelRedModule:
    DeforModelRed3D,
    DeforModelRed2DStrain,
    DeforModelRed2DStress,
    DeforModelRed2DAxisymm,
    DeforModelRed1D
using FinEtools.MatModule: AbstractMat
using LinearAlgebra

_RotationMatrix = Union{Array{T, 2}, Transpose{T, Array{T, 2}}} where {T}

"""
    AbstractMatDefor

Abstract type that represents deformable materials.
"""
abstract type AbstractMatDefor <: AbstractMat end

include("genconv.jl")
include("genrot.jl")
include("tens4impl.jl")

end # module
