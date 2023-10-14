module MatDeforElastIsoModule

__precompile__(true)

using FinEtools.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D, DeforModelRed2DStrain, DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed1DStrain, DeforModelRed1DStress, nstressstrain, nthermstrain
using FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor, stressvtot!
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: AbstractMatDeforLinearElastic
using LinearAlgebra: Transpose, Diagonal, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
using LinearAlgebra: eigen, eigvals, norm, cholesky, cross, dot


"""
	MatDeforElastIso{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function, MTHS<:Function} <: AbstractMatDeforLinearElastic

Linear isotropic elasticity  material.

"""
struct  MatDeforElastIso{MR<:AbstractDeforModelRed, FT, MTAN<:Function, MUPD<:Function, MTHS<:Function} <: AbstractMatDeforLinearElastic
	mr::Type{MR} # model reduction type
	mass_density::FT # mass density
	E::FT # Young's modulus
	nu::FT # Poisson ratio
	CTE::FT # Coefficient of Thermal Expansion
	D::Matrix{FT} # cached matrix of 3D tangent moduli
	tangentmoduli!::MTAN # Function to return the tangent moduli matrix
	update!::MUPD # Function to update the material state
	thermalstrain!::MTHS # Function to calculate the thermal strains
end

function _threedD(E, nu)
	mI = Matrix(Diagonal([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]))
	m1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
	lambda = E * nu / (1 + nu) / (1 - 2*(nu));
	mu = E / 2. / (1+nu);
	D = lambda * m1 * m1' + 2. * mu * mI;
	return D
end

"""
    MatDeforElastIso(mr::Type{MR}, mass_density, E, nu, CTE) where {MR<:AbstractDeforModelRed}

Create an isotropic elastic material providing all material parameters.
"""
function MatDeforElastIso(mr::Type{MR}, mass_density, E, nu, CTE) where {MR<:AbstractDeforModelRed}
    # Accept input data of any numerical type, then promote
    return MatDeforElastIso(mr, float.(promote(mass_density, E, nu, CTE)))
end

"""
    MatDeforElastIso(mr::Type{MR}, E, nu) where {MR}

Create isotropic elastic material with default mass density and thermal expansion.
"""
function MatDeforElastIso(mr::Type{MR}, E, nu) where {MR<:AbstractDeforModelRed}
    mass_density = 1.0
    CTE = 0.0
    return MatDeforElastIso(mr, float.(promote(mass_density, E, nu, CTE)))
end


################################################################################
# 3-D solid model
################################################################################

"""
    MatDeforElastIso(mr::Type{DeforModelRed3D}, args::NTuple{4, FT}) where FT

Create elastic isotropic material for 3D stress models.
"""
function MatDeforElastIso(mr::Type{DeforModelRed3D}, args::NTuple{4, FT}) where FT
    mass_density, E, nu, CTE = args
    function tangentmoduli3d!(self::MatDeforElastIso,  D::Matrix{FT},  t::FT, dt::FT, loc::Matrix{FT}, label::Int)
		copyto!(D,self.D);
		return D
	end
	function update3d!(self::MatDeforElastIso,  stress::Vector{FT}, output::Vector{FT},  strain::Vector{FT}, thstrain::Vector{FT}=zeros(6), t= 0.0, dt= 0.0,  loc::Matrix{FT}=zeros(3,1), label::Int=0, quantity=:nothing)
		@assert length(stress) == nstressstrain(self.mr)
		A_mul_B!(stress, self.D, strain-thstrain);
		if quantity == :nothing
			#Nothing to be copied to the output array
		elseif quantity == :cauchy || quantity == :Cauchy
			(length(output) >= 6) || (output = zeros(6)) # make sure we can store it
			copyto!(output, stress);
		elseif quantity == :pressure || quantity == :Pressure
			output[1]  =  -sum(stress[1:3])/3.
		elseif quantity == :princCauchy || quantity == :princcauchy
			t = zeros(FT,3,3)
			t = stressvtot!(mr, t, stress);
			ep = eigen(t);
			(length(output) >= 3) || (output = zeros(3)) # make sure we can store it
			copyto!(output,  sort(ep.values, rev=true));
		elseif quantity == :maxshear 
			t = zeros(FT,3,3)
			t = stressvtot!(mr, t, stress);
			ep = eigen(t);
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			s = sort(ep.values, rev=true)
			copyto!(output, s[1] - s[3]);
		elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
			s1=stress[1]; s2=stress[2]; s3=stress[3];
			s4=stress[4]; s5=stress[5]; s6=stress[6];
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
		end
		return output
	end
	function thermalstrain3d!(self::MatDeforElastIso, thstrain::Vector{FT}, dT= 0.0)
		@assert length(thstrain) == nthermstrain(self.mr)
		thstrain[1] = self.CTE*dT
		thstrain[2] = self.CTE*dT
		thstrain[3] = self.CTE*dT
		thstrain[4] = 0.0
		thstrain[5] = 0.0
		thstrain[6] = 0.0
		return thstrain
	end

	return MatDeforElastIso(mr, mass_density, E, nu, CTE, _threedD(E, nu),
		tangentmoduli3d!, update3d!, thermalstrain3d!)
end



################################################################################
# 2-D plane stress
################################################################################

"""
    MatDeforElastIso(mr::Type{DeforModelRed2DStress}, args::NTuple{4, FT}) where FT

Create elastic isotropic material for 2D plane stress models.
"""
function MatDeforElastIso(mr::Type{DeforModelRed2DStress}, args::NTuple{4, FT}) where FT
    mass_density, E, nu, CTE = args
	function tangentmoduli2dstrs!(self::MatDeforElastIso, D::Matrix{FT},  t::FT, dt::FT, loc::Matrix{FT}, label::Int)
		D[1:2, 1:2] = self.D[1:2, 1:2] -  (reshape(self.D[1:2,3], 2, 1) * reshape(self.D[3,1:2], 1, 2))/self.D[3, 3]
		ix=[1, 2, 4];
		for i = 1:3
			D[3,i] = D[i,3] = self.D[4, ix[i]];
		end
		return D
	end
	function update2dstrs!(self::MatDeforElastIso, stress::Vector{FT}, output::Vector{FT},  strain::Vector{FT}, thstrain::Vector{FT}=zeros(3), t= 0.0, dt= 0.0,  loc::Matrix{FT}=zeros(3,1), label::Int=0, quantity=:nothing)
		@assert length(stress) == nstressstrain(self.mr)
		D = zeros(3, 3)
		tangentmoduli2dstrs!(self, D, t, dt, loc, label)
		A_mul_B!(stress, D, strain-thstrain);
		if quantity == :nothing
			#Nothing to be copied to the output array
		elseif quantity == :cauchy || quantity == :Cauchy
			(length(output) >= 3) || (output = zeros(3)) # make sure we can store it
			copyto!(output, stress);
		elseif quantity == :pressure || quantity == :Pressure
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			output[1] = -sum(stress[1:2])/3.
		elseif quantity == :princCauchy || quantity == :princcauchy
			t = zeros(FT,2,2)
			t = stressvtot!(mr, t, stress);
			ep = eigen(t);
			(length(output) >= 2) || (output = zeros(2)) # make sure we can store it
			copyto!(output,  sort(ep.values, rev=true));
		elseif quantity == :maxshear 
			t = zeros(FT,2,2)
			t = stressvtot!(mr, t, stress);
			ep = eigen(t);
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			s = sort(ep.values, rev=true)
			copyto!(output, s[1] - s[2]);
		elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
			s1=stress[1]; s2=stress[2]; s3=0.0;
			s4=stress[3]; s5=0.0; s6=0.0;
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
		end
		return output
	end
	function thermalstrain2dstrs!(self::MatDeforElastIso, thstrain::Vector{FT}, dT= 0.0)
		@assert length(thstrain) == nthermstrain(self.mr)
		thstrain[1] = self.CTE*dT
		thstrain[2] = self.CTE*dT
		thstrain[3] = 0.0
		return thstrain
	end
	return MatDeforElastIso(mr, mass_density, E, nu, CTE, _threedD(E, nu),
		tangentmoduli2dstrs!, update2dstrs!, thermalstrain2dstrs!)
end


################################################################################
# 2-D plane strain
################################################################################

"""
    MatDeforElastIso(mr::Type{DeforModelRed2DStrain}, args::NTuple{4, FT}) where FT

Create elastic isotropic material for 2D plane strain models.
"""
function MatDeforElastIso(mr::Type{DeforModelRed2DStrain}, args::NTuple{4, FT}) where FT
    mass_density, E, nu, CTE = args
	function tangentmoduli2dstrn!(self::MatDeforElastIso,  D::Matrix{FT},  t, dt, loc::Matrix{FT}, label::Int)
		ix = [1, 2, 4];
		for i = 1:length(ix)
			for j = 1:length(ix)
				D[j,i] = self.D[ix[j], ix[i]];
			end
		end
		return D
	end
	# Note on the principal stresses: The principal stresses are calculated for the
	# fully three-dimensional  stress state, that is not the "in-plane" maximum and
	# minimum,  but rather  the  three-dimensional maximum (1) and minimum (3).
	# The intermediate principal stress is (2).
	function update2dstrn!(self::MatDeforElastIso,  stress::Vector{FT}, output::Vector{FT},  strain::Vector{FT}, thstrain::Vector{FT}=zeros(4), t= 0.0, dt= 0.0,  loc::Matrix{FT}=zeros(3,1), label::Int=0, quantity=:nothing)
		@assert length(stress) == nstressstrain(self.mr)
		D = zeros(3, 3)
		tangentmoduli2dstrn!(self, D, t, dt, loc, label)
		A_mul_B!(stress, D, strain-thstrain[1:3]);
		if quantity == :nothing
			#Nothing to be copied to the output array
		elseif quantity == :cauchy || quantity == :Cauchy
			# sigmax, sigmay, tauxy, sigmaz
			# thstrain[4] =The through the thickness thermal strain
			sz = dot(self.D[3, 1:2], strain[1:2]-thstrain[1:2])-self.D[3,3]*thstrain[4];
			(length(output) >= 4) || (output = zeros(4)) # make sure we can store it
			copyto!(output, stress); output[4] = sz
		elseif quantity == :pressure || quantity == :Pressure
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			sz = dot(self.D[3, 1:2], strain[1:2]-thstrain[1:2])-self.D[3,3]*thstrain[4];
			output[1] = -(sum(stress[[1,2]]) + sz)/3.
		elseif quantity == :princCauchy || quantity == :princcauchy
			(length(output) >= 3) || (output = zeros(3)) # make sure we can store it
			t = zeros(FT, 3,3)
			sz = dot(self.D[3, 1:2], strain[1:2]-thstrain[1:2])-self.D[3,3]*thstrain[4];
			t = stressvtot!(mr, t, vcat(stress[1:3], [sz]));
			ep = eigen(t);
			(length(output) >= 3) || (output = zeros(3)) # make sure we can store it
			copyto!(output,  sort(ep.values, rev=true));
		elseif quantity == :maxshear 
			(length(output) >= 3) || (output = zeros(3)) # make sure we can store it
			t = zeros(FT, 3,3)
			sz = dot(self.D[3, 1:2], strain[1:2]-thstrain[1:2])-self.D[3,3]*thstrain[4];
			t = stressvtot!(mr, t, vcat(stress[1:3], [sz]));
			ep = eigen(t);
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			s = sort(ep.values, rev=true)
			copyto!(output, s[1] - s[3]);
		elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			sz = dot(self.D[3, 1:2], strain[1:2]-thstrain[1:2])-self.D[3,3]*thstrain[4];
			s1=stress[1]; s2=stress[2]; s3=sz;
			s4=stress[3]; s5=0.0; s6=0.0;
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
		end
		return output
	end
	function thermalstrain2dstrn!(self::MatDeforElastIso, thstrain::Vector{FT}, dT= 0.0)
		@assert length(thstrain) == nthermstrain(self.mr)
		thstrain[1] = self.CTE*dT
		thstrain[2] = self.CTE*dT
		thstrain[3] = 0.0
		thstrain[4] = self.CTE*dT
		return thstrain
	end
	return MatDeforElastIso(mr, mass_density, E, nu, CTE, _threedD(E, nu),
		tangentmoduli2dstrn!, update2dstrn!, thermalstrain2dstrn!)
end


################################################################################
# 2-D axially symmetric
################################################################################

"""
    MatDeforElastIso(mr::Type{DeforModelRed2DAxisymm}, args::NTuple{4, FT}) where FT

Create elastic isotropic material for 2D axially symmetric models.
"""
function MatDeforElastIso(mr::Type{DeforModelRed2DAxisymm}, args::NTuple{4, FT}) where FT
    mass_density, E, nu, CTE = args
	function tangentmoduli2daxi!(self::MatDeforElastIso,  D::Matrix{FT},  t, dt, loc::Matrix{FT}, label::Int)
		for i = 1:4
			for j = 1:4
				D[j,i] = self.D[i, j];
			end
		end
		return D
	end
	function update2daxi!(self::MatDeforElastIso,  stress::Vector{FT}, output::Vector{FT},  strain::Vector{FT}, thstrain::Vector{FT}=zeros(3), t= 0.0, dt= 0.0,  loc::Matrix{FT}=zeros(3,1), label::Int=0, quantity=:nothing)
		@assert length(stress) == nstressstrain(self.mr)
		D = zeros(4, 4)
		tangentmoduli2daxi!(self, D, t, dt, loc, label)
		A_mul_B!(stress, D, strain-thstrain);
		if quantity == :nothing
			#Nothing to be copied to the output array
		elseif quantity == :cauchy || quantity == :Cauchy
			(length(output) >= 4) || (output = zeros(4)) # make sure we can store it
			copyto!(output, stress);
		elseif quantity == :pressure || quantity == :Pressure
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			output[1] = -sum(stress[[1,2,3]])/3.
		elseif quantity == :princCauchy || quantity == :princcauchy
			t = zeros(FT,3,3)
			t = stressvtot!(mr, t, stress);
			ep = eigen(t);
			(length(output) >= 3) || (output = zeros(3)) # make sure we can store it
			copyto!(output,  sort(ep.values, rev=true));
		elseif quantity == :maxshear 
			t = zeros(FT,3,3)
			t = stressvtot!(mr, t, stress);
			ep = eigen(t);
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			s = sort(ep.values, rev=true)
			copyto!(output, s[1] - s[3]);
		elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
			s1=stress[1]; s2=stress[2]; s3=stress[3];
			s4=stress[4]; s5=0.0; s6=0.0;
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
		end
		return output
	end
	function thermalstrain2daxi!(self::MatDeforElastIso, thstrain::Vector{FT}, dT= 0.0)
		@assert length(thstrain) == nthermstrain(self.mr)
		thstrain[1] = self.CTE*dT
		thstrain[2] = self.CTE*dT
		thstrain[3] = self.CTE*dT
		thstrain[4] = 0.0
		return thstrain
	end
	return MatDeforElastIso(mr, mass_density, E, nu, CTE, _threedD(E, nu),
		tangentmoduli2daxi!, update2daxi!, thermalstrain2daxi!)
end

################################################################################
# 1-D zero transverse strain
################################################################################

"""
    MatDeforElastIso(mr::Type{DeforModelRed1DStrain}, args::NTuple{4, FT}) where FT

Create elastic isotropic material for 1D models.
"""
function MatDeforElastIso(mr::Type{DeforModelRed1DStrain}, args::NTuple{4, FT}) where FT
    mass_density, E, nu, CTE = args
    function tangentmoduli1d!(self::MatDeforElastIso, D::Matrix{FT}, t, dt, loc::Matrix{FT}, label::Int)
        lambda = E * nu / (1 + nu) / (1 - 2*(nu));
        mu = E / 2. / (1+nu);
        D[1, 1] = lambda + 2*mu;
        return D
    end
    function update1d!(self::MatDeforElastIso,  stress::Vector{FT}, output::Vector{FT},  strain::Vector{FT}, thstrain::Vector{FT}=zeros(1), t= 0.0, dt= 0.0,  loc::Matrix{FT}=zeros(3,1), label::Int=0, quantity=:nothing)
        @assert length(stress) == nstressstrain(self.mr)
        D = zeros(1, 1)
        tangentmoduli1d!(self, D, t, dt, loc, label)
        A_mul_B!(stress, D, strain-thstrain);
        if quantity == :nothing
            #Nothing to be copied to the output array
        elseif quantity == :cauchy || quantity == :Cauchy
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            copyto!(output, stress);
        elseif quantity == :pressure || quantity == :Pressure
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = -sum(stress[[1]])/3.
        elseif quantity == :princCauchy || quantity == :princcauchy
            copyto!(output,  stress[1]);
        elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
            s1=stress[1]; s2=0.0; s3=0.0;
            s4=0.0; s5=0.0; s6=0.0;
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
        end
        return output
    end
    function thermalstrain1d!(self::MatDeforElastIso, thstrain::Vector{FT}, dT= 0.0)
        @assert length(thstrain) == nthermstrain(self.mr)
        thstrain[1] = self.CTE*dT
        return thstrain
    end
    return MatDeforElastIso(mr, mass_density, E, nu, CTE, _threedD(E, nu),
        tangentmoduli1d!, update1d!, thermalstrain1d!)
end

################################################################################
# 1-D zero transverse stress
################################################################################

"""
    MatDeforElastIso(mr::Type{DeforModelRed1DStress}, args::NTuple{4, FT}) where FT

Create elastic isotropic material for 1D models.
"""
function MatDeforElastIso(mr::Type{DeforModelRed1DStress}, args::NTuple{4, FT}) where FT
    mass_density, E, nu, CTE = args
    function tangentmoduli1d!(self::MatDeforElastIso, D::Matrix{FT}, t, dt, loc::Matrix{FT}, label::Int)
        D3d = _threedD(E, nu)
        D[1, 1] = D3d[1, 1] - transpose(D3d[1, 2:3])*(D3d[2:3, 2:3]\D3d[2:3, 1]);
        return D
    end
    function update1d!(self::MatDeforElastIso,  stress::Vector{FT}, output::Vector{FT},  strain::Vector{FT}, thstrain::Vector{FT}=zeros(1), t= 0.0, dt= 0.0,  loc::Matrix{FT}=zeros(3,1), label::Int=0, quantity=:nothing)
        @assert length(stress) == nstressstrain(self.mr)
        D = zeros(1, 1)
        tangentmoduli1d!(self, D, t, dt, loc, label)
        A_mul_B!(stress, D, strain-thstrain);
        if quantity == :nothing
            #Nothing to be copied to the output array
        elseif quantity == :cauchy || quantity == :Cauchy
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            copyto!(output, stress);
        elseif quantity == :pressure || quantity == :Pressure
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = -sum(stress[[1]])/3.
        elseif quantity == :princCauchy || quantity == :princcauchy
            copyto!(output,  stress[1]);
        elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
            s1=stress[1]; s2=0.0; s3=0.0;
            s4=0.0; s5=0.0; s6=0.0;
            (length(output) >= 1) || (output = zeros(1)) # make sure we can store it
            output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
        end
        return output
    end
    function thermalstrain1d!(self::MatDeforElastIso, thstrain::Vector{FT}, dT= 0.0)
        @assert length(thstrain) == nthermstrain(self.mr)
        thstrain[1] = self.CTE*dT
        return thstrain
    end
    return MatDeforElastIso(mr, mass_density, E, nu, CTE, _threedD(E, nu),
        tangentmoduli1d!, update1d!, thermalstrain1d!)
end

end
