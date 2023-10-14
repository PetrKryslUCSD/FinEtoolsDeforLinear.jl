"""
    MatDeforElastOrthoModule

Module for  orthotropic elastic material.
"""
module MatDeforElastOrthoModule

__precompile__(true)

using FinEtools.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D,
    DeforModelRed2DStrain, DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed1D,
    nstressstrain, nthermstrain
using FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor, stressvtot!
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: AbstractMatDeforLinearElastic
using LinearAlgebra: mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
using LinearAlgebra: eigen, eigvals, rank, dot

"""
    MatDeforElastOrtho{MR<:AbstractDeforModelRed,  MTAN<:Function, MUPD<:Function, MTHS<:Function} <: AbstractMatDeforLinearElastic

Linear orthotropic elasticity  material.
"""
struct MatDeforElastOrtho{
    MR <: AbstractDeforModelRed,
    FT,
    MTAN <: Function,
    MUPD <: Function,
    MTHS <: Function,
} <: AbstractMatDeforLinearElastic
    mr::Type{MR}
    mass_density::FT # mass density
    E1::FT # Young's modulus for material direction 1
    E2::FT # Young's modulus for material direction 2
    E3::FT # Young's modulus for material direction 3
    nu12::FT
    nu13::FT
    nu23::FT
    G12::FT
    G13::FT
    G23::FT
    CTE1::FT # three thermal expansion coefficients
    CTE2::FT # three thermal expansion coefficients
    CTE3::FT # three thermal expansion coefficients
    D::Matrix{FT} # cached matrix of tangent moduli
    tangentmoduli!::MTAN
    update!::MUPD
    thermalstrain!::MTHS
end

function _threedD(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23)
    C = [1.0/E1 -nu12/E1 -nu13/E1 0.0 0.0 0.0;
        -nu12/E1 1.0/E2 -nu23/E2 0.0 0.0 0.0;
        -nu13/E1 -nu23/E2 1.0/E3 0.0 0.0 0.0;
        0.0 0.0 0.0 1/G12 0.0 0.0;
        0.0 0.0 0.0 0.0 1/G13 0.0;
        0.0 0.0 0.0 0.0 0.0 1/G23]
    D = inv(C)
    if (rank(D) < 6)
        error("Non-positive definite D!")
    end
    ev = eigvals(D)
    for e in ev
        #println("$e")
        if (e < 0.0)
            error("Indefinite D!")
        end
    end
    return D
end

"""
    MatDeforElastOrtho(mr::Type{MR}, mass_density, E1, E2, E3,
        nu12, nu13, nu23,
        G12, G13, G23, CTE1, CTE2, CTE3) where {MR}

Create elastic orthotropic material.
"""
function MatDeforElastOrtho(mr::Type{MR}, mass_density, E1, E2, E3,
    nu12, nu13, nu23,
    G12, G13, G23, CTE1, CTE2, CTE3) where {MR}
    return MatDeforElastOrtho(mr,
        float.(promote(mass_density, E1, E2, E3, nu12, nu13, nu23,
            G12, G13, G23, CTE1, CTE2, CTE3)))
end

"""
	MatDeforElastOrtho(mr::Type{MR}, E1, E2, E3,
        nu12, nu13, nu23,
        G12, G13, G23) where {MR}

Create elastic orthotropic material.

Convenience version with only the specification of the elastic properties.
"""
function MatDeforElastOrtho(mr::Type{MR}, E1, E2, E3,
    nu12, nu13, nu23,
    G12, G13, G23) where {MR}
    mass_density = 1.0
    CTE1 = CTE2 = CTE3 = 0.0
    return MatDeforElastOrtho(mr,
        float.(promote(mass_density, E1, E2, E3, nu12, nu13, nu23,
            G12, G13, G23, CTE1, CTE2, CTE3)))
end

"""
	MatDeforElastOrtho(mr::Type{MR}, E, nu) where {MR}

Create elastic orthotropic material which is really isotropic.

Convenience version with only the specification of the elastic properties.
"""
function MatDeforElastOrtho(mr::Type{MR}, E, nu) where {MR}
    mass_density = 1.0
    E1 = E2 = E3 = E
    nu12 = nu13 = nu23 = nu
    CTE1 = CTE2 = CTE3 = 0.0
    G = E / 2.0 / (1 + nu)
    G12 = G13 = G23 = G
    return MatDeforElastOrtho(mr,
        float.(promote(mass_density, E1, E2, E3, nu12, nu13, nu23,
            G12, G13, G23, CTE1, CTE2, CTE3)))
end

"""
	MatDeforElastOrtho(mr::Type{MR}, mass_density::FT,  E::FT, nu::FT, CTE::FT) where {MR}

Create elastic orthotropic material which is really isotropic.

Convenience version with only the specification of the elastic and thermal expansion properties.
"""
function MatDeforElastOrtho(mr::Type{MR}, mass_density, E, nu, CTE) where {MR}
    mass_density = 1.0
    E1 = E2 = E3 = E
    nu12 = nu13 = nu23 = nu
    CTE1 = CTE2 = CTE3 = CTE
    G = E / 2.0 / (1 + nu)
    G12 = G13 = G23 = G
    return MatDeforElastOrtho(mr,
        float.(promote(mass_density, E1, E2, E3, nu12, nu13, nu23,
            G12, G13, G23, CTE1, CTE2, CTE3)))
end

################################################################################
# 3-D solid model
################################################################################

"""
	MatDeforElastOrtho(mr::Type{DeforModelRed3D}, args::NTuple{13, FT}) where FT

Create elastic orthotropic material for 3D models.
"""
function MatDeforElastOrtho(mr::Type{DeforModelRed3D}, args::NTuple{13, FT}) where {FT}
    mass_density, E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, CTE1, CTE2, CTE3 = args
    function tangentmoduli3d!(self::MatDeforElastOrtho,
        D::Matrix{FT},
        t::FT,
        dt::FT,
        loc::Matrix{FT},
        label::Int)
        copyto!(D, self.D)
        return D
    end
    function update3d!(self::MatDeforElastOrtho,
        stress::Vector{FT},
        output::Vector{FT},
        strain::Vector{FT},
        thstrain::Vector{FT} = zeros(6),
        t::FT = 0.0,
        dt::FT = 0.0,
        loc::Matrix{FT} = zeros(3, 1),
        label::Int = 0,
        quantity = :nothing)
        @assert length(stress) == nstressstrain(self.mr)
        A_mul_B!(stress, self.D, strain - thstrain)
        if quantity == :nothing
            #Nothing to be copied to the output array
        elseif quantity == :cauchy || quantity == :Cauchy
            (length(output) >= 6) || (output = zeros(6)) # make sure we can store it
            copyto!(output, stress)
        elseif quantity == :pressure || quantity == :Pressure
            output[1] = -sum(stress[1:3]) / 3.0
        elseif quantity == :princCauchy || quantity == :princcauchy
            t = zeros(FT, 3, 3)
            t = stressvtot!(mr, t, stress)
            ep = eigen(t)
            (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
            copyto!(output, sort(ep.values, rev = true))
        elseif quantity == :vonMises || quantity == :vonmises || quantity == :von_mises ||
               quantity == :vm
            s1 = stress[1]
            s2 = stress[2]
            s3 = stress[3]
            s4 = stress[4]
            s5 = stress[5]
            s6 = stress[6]
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = sqrt(1.0 / 2 * ((s1 - s2)^2 + (s1 - s3)^2 + (s2 - s3)^2 +
                              6 * (s4^2 + s5^2 + s6^2)))
        end
        return output
    end
    function thermalstrain3d!(self::MatDeforElastOrtho, thstrain::Vector{FT}, dT = 0.0)
        thstrain[1] = self.CTE1 * dT
        thstrain[2] = self.CTE2 * dT
        thstrain[3] = self.CTE3 * dT
        thstrain[4] = 0.0
        thstrain[5] = 0.0
        thstrain[6] = 0.0
        return thstrain
    end
    return MatDeforElastOrtho(mr, mass_density, E1, E2, E3, nu12, nu13, nu23,
        G12, G13, G23, CTE1, CTE2, CTE3,
        _threedD(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23),
        tangentmoduli3d!, update3d!, thermalstrain3d!)
end

################################################################################
# 2-D plane stress
################################################################################

"""
	MatDeforElastOrtho(mr::Type{DeforModelRed2DStress}, args::NTuple{13, FT}) where FT

Create elastic orthotropic material for 2D plane stress models.
"""
function MatDeforElastOrtho(mr::Type{DeforModelRed2DStress},
    args::NTuple{13, FT}) where {FT}
    mass_density, E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, CTE1, CTE2, CTE3 = args
    function tangentmoduli2dstrs!(self::MatDeforElastOrtho,
        D::Matrix{FT},
        t::FT,
        dt::FT,
        loc::Matrix{FT},
        label::Int)
        D[1:2, 1:2] = self.D[1:2, 1:2] -
                      (reshape(self.D[1:2, 3], 2, 1) * reshape(self.D[3, 1:2], 1, 2)) /
                      self.D[3, 3]
        ix = [1, 2, 4]
        for i in 1:3
            D[3, i] = D[i, 3] = self.D[4, ix[i]]
        end
        return D
    end
    function update2dstrs!(self::MatDeforElastOrtho,
        stress::Vector{FT},
        output::Vector{FT},
        strain::Vector{FT},
        thstrain::Vector{FT} = zeros(3),
        t::FT = 0.0,
        dt::FT = 0.0,
        loc::Matrix{FT} = zeros(3, 1),
        label::Int = 0,
        quantity = :nothing)
        @assert length(stress) == nstressstrain(self.mr)
        D = zeros(3, 3)
        tangentmoduli2dstrs!(self, D, t, dt, loc, label)
        A_mul_B!(stress, D, strain - thstrain)
        if quantity == :nothing
            #Nothing to be copied to the output array
        elseif quantity == :cauchy || quantity == :Cauchy
            (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
            copyto!(output, stress)
        elseif quantity == :pressure || quantity == :Pressure
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = -sum(stress[1:2]) / 3.0
        elseif quantity == :princCauchy || quantity == :princcauchy
            t = zeros(FT, 2, 2)
            t = stressvtot!(mr, t, stress)
            ep = eigen(t)
            (length(output) >= 2) || (output = zeros(2)) # make sure we can store it
            copyto!(output, sort(ep.values, rev = true))
        elseif quantity == :vonMises || quantity == :vonmises || quantity == :von_mises ||
               quantity == :vm
            s1 = stress[1]
            s2 = stress[2]
            s3 = 0.0
            s4 = stress[3]
            s5 = 0.0
            s6 = 0.0
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = sqrt(1.0 / 2 * ((s1 - s2)^2 + (s1 - s3)^2 + (s2 - s3)^2 +
                              6 * (s4^2 + s5^2 + s6^2)))
        end
        return output
    end
    function thermalstrain2dstrs!(self::MatDeforElastOrtho, thstrain::Vector{FT}, dT = 0.0)
        @assert length(thstrain) == nthermstrain(self.mr)
        thstrain[1] = self.CTE1 * dT
        thstrain[2] = self.CTE2 * dT
        thstrain[3] = 0.0
        return thstrain
    end
    return MatDeforElastOrtho(mr, mass_density, E1, E2, E3, nu12, nu13, nu23,
        G12, G13, G23, CTE1, CTE2, CTE3,
        _threedD(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23),
        tangentmoduli2dstrs!, update2dstrs!, thermalstrain2dstrs!)
end

################################################################################
# 2-D plane strain
################################################################################

"""
	MatDeforElastOrtho(mr::Type{DeforModelRed2DStrain}, args::NTuple{13, FT}) where FT

Create elastic orthotropic material for 2D plane strain models.
"""
function MatDeforElastOrtho(mr::Type{DeforModelRed2DStrain},
    args::NTuple{13, FT}) where {FT}
    mass_density, E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, CTE1, CTE2, CTE3 = args
    function tangentmoduli2dstrn!(self::MatDeforElastOrtho,
        D::Matrix{FT},
        t::FT,
        dt::FT,
        loc::Matrix{FT},
        label::Int)
        ix = [1, 2, 4]
        for i in 1:length(ix)
            for j in 1:length(ix)
                D[j, i] = self.D[ix[j], ix[i]]
            end
        end
        return D
    end
    function update2dstrn!(self::MatDeforElastOrtho,
        stress::Vector{FT},
        output::Vector{FT},
        strain::Vector{FT},
        thstrain::Vector{FT} = zeros(4),
        t::FT = 0.0,
        dt::FT = 0.0,
        loc::Matrix{FT} = zeros(3, 1),
        label::Int = 0,
        quantity = :nothing)
        @assert length(stress) == nstressstrain(self.mr)
        D = zeros(3, 3)
        tangentmoduli2dstrn!(self, D, t, dt, loc, label)
        A_mul_B!(stress, D, strain - thstrain[1:3])
        if quantity == :nothing
            #Nothing to be copied to the output array
        elseif quantity == :cauchy || quantity == :Cauchy
            # sigmax, sigmay, tauxy, sigmaz
            # thstrain[4] =The through the thickness thermal strain
            sz = dot(self.D[3, 1:2], strain[1:2] - thstrain[1:2]) -
                 self.D[3, 3] * thstrain[4]
            (length(output) >= 4) || (output = zeros(4)) # make sure we can store it
            copyto!(output, stress)
            output[4] = sz
        elseif quantity == :pressure || quantity == :Pressure
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            sz = dot(self.D[3, 1:2], strain[1:2] - thstrain[1:2]) -
                 self.D[3, 3] * thstrain[4]
            output[1] = -(sum(stress[[1, 2]]) + sz) / 3.0
        elseif quantity == :princCauchy || quantity == :princcauchy
            (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
            t = zeros(FT, 3, 3)
            sz = dot(self.D[3, 1:2], strain[1:2] - thstrain[1:2]) -
                 self.D[3, 3] * thstrain[4]
            t = stressvtot!(mr, t, vcat(stress[1:3], [sz]))
            ep = eigen(t)
            (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
            copyto!(output, sort(ep.values, rev = true))
        elseif quantity == :vonMises || quantity == :vonmises || quantity == :von_mises ||
               quantity == :vm
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            sz = dot(self.D[3, 1:2], strain[1:2] - thstrain[1:2]) -
                 self.D[3, 3] * thstrain[4]
            s1 = stress[1]
            s2 = stress[2]
            s3 = sz
            s4 = stress[3]
            s5 = 0.0
            s6 = 0.0
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = sqrt(1.0 / 2 * ((s1 - s2)^2 + (s1 - s3)^2 + (s2 - s3)^2 +
                              6 * (s4^2 + s5^2 + s6^2)))
        end
        return output
    end
    function thermalstrain2dstrn!(self::MatDeforElastOrtho, thstrain::Vector{FT}, dT = 0.0)
        @assert length(thstrain) == nthermstrain(self.mr)
        thstrain[1] = self.CTE1 * dT
        thstrain[2] = self.CTE2 * dT
        thstrain[3] = 0.0
        thstrain[4] = self.CTE3 * dT
        return thstrain
    end
    return MatDeforElastOrtho(mr, mass_density, E1, E2, E3, nu12, nu13, nu23,
        G12, G13, G23, CTE1, CTE2, CTE3,
        _threedD(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23),
        tangentmoduli2dstrn!, update2dstrn!, thermalstrain2dstrn!)
end

################################################################################
# 2-D axially symmetric
################################################################################

"""
	MatDeforElastOrtho(mr::Type{DeforModelRed2DAxisymm}, args::NTuple{13, FT}) where FT

Create elastic orthotropic material for 2D axially symmetric models.
"""
function MatDeforElastOrtho(mr::Type{DeforModelRed2DAxisymm},
    args::NTuple{13, FT}) where {FT}
    mass_density, E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, CTE1, CTE2, CTE3 = args
    function tangentmoduli2daxi!(self::MatDeforElastOrtho,
        D::Matrix{FT},
        t::FT,
        dt::FT,
        loc::Matrix{FT},
        label::Int)
        for i in 1:4
            for j in 1:4
                D[i, j] = self.D[i, j]
            end
        end
        return D
    end
    function update2daxi!(self::MatDeforElastOrtho,
        stress::Vector{FT},
        output::Vector{FT},
        strain::Vector{FT},
        thstrain::Vector{FT} = zeros(3),
        t::FT = 0.0,
        dt::FT = 0.0,
        loc::Matrix{FT} = zeros(3, 1),
        label::Int = 0,
        quantity = :nothing)
        @assert length(stress) == nstressstrain(self.mr)
        D = zeros(4, 4)
        tangentmoduli2daxi!(self, D, t, dt, loc, label)
        A_mul_B!(stress, D, strain - thstrain)
        if quantity == :nothing
            #Nothing to be copied to the output array
        elseif quantity == :cauchy || quantity == :Cauchy
            (length(output) >= 4) || (output = zeros(4)) # make sure we can store it
            copyto!(output, stress)
        elseif quantity == :pressure || quantity == :Pressure
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = -sum(stress[[1, 2, 3]]) / 3.0
        elseif quantity == :princCauchy || quantity == :princcauchy
            t = zeros(FT, 3, 3)
            t = stressvtot!(mr, t, stress)
            ep = eigen(t)
            (length(output) >= 3) || (output = zeros(3)) # make sure we can store it
            copyto!(output, sort(ep.values, rev = true))
        elseif quantity == :vonMises || quantity == :vonmises || quantity == :von_mises ||
               quantity == :vm
            s1 = stress[1]
            s2 = stress[2]
            s3 = stress[3]
            s4 = stress[4]
            s5 = 0.0
            s6 = 0.0
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = sqrt(1.0 / 2 * ((s1 - s2)^2 + (s1 - s3)^2 + (s2 - s3)^2 +
                              6 * (s4^2 + s5^2 + s6^2)))
        end
        return output
    end
    function thermalstrain2daxi!(self::MatDeforElastOrtho, thstrain::Vector{FT}, dT = 0.0)
        @assert length(thstrain) == nthermstrain(self.mr)
        thstrain[1] = self.CTE1 * dT
        thstrain[2] = self.CTE2 * dT
        thstrain[3] = self.CTE3 * dT
        thstrain[4] = 0.0
        return thstrain
    end
    return MatDeforElastOrtho(mr, mass_density, E1, E2, E3, nu12, nu13, nu23,
        G12, G13, G23, CTE1, CTE2, CTE3,
        _threedD(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23),
        tangentmoduli2daxi!, update2daxi!, thermalstrain2daxi!)
end

end
