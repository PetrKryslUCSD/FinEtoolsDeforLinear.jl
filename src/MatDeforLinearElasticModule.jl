module MatDeforLinearElasticModule

__precompile__(true)

using FinEtools.DeforModelRedModule: AbstractDeforModelRed
using FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor

"""
    AbstractMatDeforLinearElastic <: AbstractMatDefor

Abstract Linear Elasticity  material.

"""
abstract type AbstractMatDeforLinearElastic <: AbstractMatDefor end

"""
    tangentmoduli!(self::AbstractMatDeforLinearElastic,  D::Matrix{FT},  t::FT, dt::FT, loc::Matrix{FT}, label::Int) where {FT}

Calculate the material stiffness matrix.

- `D` = matrix of tangent moduli, supplied as a buffer and overwritten. Returned
as output.
"""
function tangentmoduli!(self::AbstractMatDeforLinearElastic,
    D::Matrix{FT},
    t::FT,
    dt::FT,
    loc::Matrix{FT},
    label::Int) where {FT}
    return self.tangentmoduli!(self, D, t, dt, loc, label)
end

"""
    update!(self::AbstractMatDeforLinearElastic,  stress::Vector{FT}, output::Vector{FT},  strain::Vector{FT}, thstrain::Vector{FT}=zeros(6), t::FT= 0.0, dt::FT= 0.0,  loc::Matrix{FT}=zeros(3,1), label::Int=0, quantity=:nothing) where {FT}

Update material state.

- `strain` = strain vector,
- `thstrain` = thermal strain vector,
- `t` = current time,
- `dt` = current time step,
- `loc` = location of the quadrature point in global Cartesian coordinates,
- `label` = label of the finite element in which the quadrature point is found.

# Output
- `stress` = stress vector, allocated by the caller with a size of the number of stress and
strain components, `nstressstrain`. The components of the stress vector are
calculated and stored in the `stress` vector.
- `output` =  array which is (if necessary) allocated  in an appropriate size, filled
  with the output quantity, and returned.
"""
function update!(self::AbstractMatDeforLinearElastic,
    stress::Vector{FT},
    output::Vector{FT},
    strain::Vector{FT},
    thstrain::Vector{FT} = zeros(6),
    t::FT = 0.0,
    dt::FT = 0.0,
    loc::Matrix{FT} = zeros(3, 1),
    label::Int = 0,
    quantity = :nothing) where {FT}
    return self.update!(self, stress, output, strain, thstrain, t, dt, loc, label, quantity)
end

"""
    thermalstrain!(self::AbstractMatDeforLinearElastic, thstrain::Vector{FT}, dT= 0.0) where {FT}

Compute thermal strain from the supplied temperature increment.

- `thstrain` = thermal strain vector, supplied as buffer, returned as output.
"""
function thermalstrain!(self::AbstractMatDeforLinearElastic,
    thstrain::Vector{FT},
    dT = 0.0) where {FT}
    return self.thermalstrain!(self, thstrain, dT)
end

end
