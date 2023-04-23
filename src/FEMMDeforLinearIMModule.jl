"""
    FEMMDeforLinearIMModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models: incompatible-mode formulation.
"""
module FEMMDeforLinearIMModule

__precompile__(true)

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, FESetH8, FESetT10, manifdim, nodesperelem, gradN!, bfun, bfundpar
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.IntegRuleModule: GaussRule
import FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: AbstractFEMMDeforLinear
import FinEtoolsDeforLinear.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D
import FinEtoolsDeforLinear.MatDeforLinearElasticModule: AbstractMatDeforLinearElastic, tangentmoduli!, update!, thermalstrain!
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
import FinEtools.NodalFieldModule: NodalField
import FinEtools.CSysModule: CSys, updatecsmat!, csmat
import FinEtoolsDeforLinear.DeforModelRedModule: nstressstrain, nthermstrain, Blmat!, divmat, vgradmat
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, loc!, jac!, locjac!
import FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: stiffness, nzebcloadsstiffness, mass, thermalstrainloads, inspectintegpoints
import FinEtools.FEMMBaseModule: associategeometry!
import FinEtoolsDeforLinear.MatDeforModule: rotstressvec!
import LinearAlgebra: mul!, Transpose, UpperTriangular
using LinearAlgebra: Symmetric, cholesky, eigen
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: norm, qr, diag, dot, cond
import Statistics: mean

"""
    FEMMDeforLinearIMH8{MR<:AbstractDeforModelRed, S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic} 

Type for mean-strain linear deformation FEMM based on eight-node hexahedral elements.
"""
mutable struct FEMMDeforLinearIMH8{MR<:AbstractDeforModelRed, S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic}  <: AbstractFEMMDeforLinear 
    mr::Type{MR}
    integdomain::IntegDomain{S, F} # geometry data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
    nmodes::Int64
end

function FEMMDeforLinearIMH8(mr::Type{MR}, integdomain::IntegDomain{S, F}, mcsys::CSys, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    return FEMMDeforLinearIMH8(mr, integdomain, mcsys, material, 12)
end

function FEMMDeforLinearIMH8(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    return FEMMDeforLinearIMH8(mr, integdomain, CSys(manifdim(integdomain.fes)), material, 12)
end

function FEMMDeforLinearIMH8(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M, nmodes::Int64) where {MR<:AbstractDeforModelRed,  S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    return FEMMDeforLinearIMH8(mr, integdomain, CSys(manifdim(integdomain.fes)), material, nmodes)
end

function  centroidintegrationdata(self) 
	integration_rule = GaussRule(3, 1)
    pc::FFltMat = integration_rule.param_coords;
    w::FFltMat  =  integration_rule.weights ;
    npts::FInt = integration_rule.npts;
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    Ns = FFltMat[];
    gradNparams = FFltMat[];
    for j=1:npts
        push!(Ns, bfun(self.fes,vec(pc[j,:])));
        push!(gradNparams, bfundpar(self.fes,vec(pc[j,:])));
    end
    @assert npts == 1
    return npts, reshape(Ns, 1,npts), reshape(gradNparams, 1,npts), w, pc
end

function imintegrationdata(nmodes, integration_rule) 
    pc::FFltMat = integration_rule.param_coords;
    w::FFltMat  =  integration_rule.weights ;
    npts::FInt = integration_rule.npts;
    bfun(nmodes, pc) = begin
        if nmodes == 12 # Simo basis functions
        	N = [0.5 * (pc[1]^2 - 1);
        		 0.5 * (pc[2]^2 - 1);
        		 0.5 * (pc[3]^2 - 1);
        		 pc[1] * pc[2] * pc[3]]
        else # Wilson basis functions
        	N = [0.5 * (pc[1]^2 - 1);
        		 0.5 * (pc[2]^2 - 1);
        		 0.5 * (pc[3]^2 - 1)]
        end
        return reshape(N, length(N), 1)
    end
    bfundpar(nmodes, pc) = begin
        if nmodes == 12 # Simo basis functions
        	gradN = [pc[1] 0 0;
	        		 0 pc[2] 0;
	        		 0 0 pc[3];
	        		 (pc[2] * pc[3]) (pc[1] * pc[3]) (pc[1] * pc[2])]
        else # Wilson basis functions
        	gradN = [pc[1] 0 0;
	        		 0 pc[2] 0;
	        		 0 0 pc[3];]
        end
        return gradN
    end
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    Ns = FFltMat[];
    gradNparams = FFltMat[];
    for j=1:npts
        push!(Ns, bfun(nmodes, vec(pc[j,:])));
        push!(gradNparams, bfundpar(nmodes, vec(pc[j,:])));
    end
    return reshape(Ns, 1, npts), reshape(gradNparams, 1, npts)
end

function imBlmat!(mr, imB, imNs, imgradN, loc0, csmat, nmodes)
	Blmat!(mr, imB, imNs, imgradN, loc0, csmat);
	if nmodes == 12
		imB[1:3, 10] .= imgradN[4, 1]
		imB[1:3, 11] .= imgradN[4, 2]
		imB[1:3, 12] .= imgradN[4, 3]
		imB[4:6, 10:12] .= 0.0
	end
end

function _buffers2(self, geom::NodalField, u::NodalField)
	fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    if self.nmodes == 12
    	nimne = 4
    else
    	nimne = 3
    end
    elmatdim = ndn*(nne + nimne); # dimension of the element matrix
    elmatcdim = ndn*nne;  # dimension of the condensed element matrix
    # Prepare buffers
    ecoords = fill(zero(FFlt), nne, ndofs(geom)); # array of Element coordinates
    elmat = fill(zero(FFlt), elmatdim, elmatdim);      # element matrix -- buffer
    elmatc = fill(zero(FFlt), elmatcdim, elmatcdim);      # element matrix -- buffer
    dofnums = zeros(FInt, 1, elmatcdim); # degree of freedom array -- buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    loc0 = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J0 = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ0 = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    gradN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    gradN0 = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    imgradN = fill(zero(FFlt), nimne, mdim); # intermediate result -- buffer
    D = fill(zero(FFlt), nstrs, nstrs); # material stiffness matrix -- buffer
    B = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    DB = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    Bc = fill(zero(FFlt), nstrs, elmatcdim); # strain-displacement matrix -- buffer
    imB = fill(zero(FFlt), nstrs, elmatdim - elmatcdim); # strain-displacement matrix -- buffer
    return ecoords, dofnums, loc, J, csmatTJ, loc0, J0, csmatTJ0, gradN, gradN0, imgradN, D, B, DB, Bc, imB, elmatc, elmat
end

"""
    associategeometry!(self::FEMMAbstractBase,  geom::NodalField{FFlt})

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(self::F,  geom::NodalField{FFlt}) where {F<:FEMMDeforLinearIMH8}
    # Nothing needs to be done
    return self
end

"""
stiffness(self::AbstractFEMMDeforLinearIM, assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute and assemble  stiffness matrix.
"""
function stiffness(self::FEMMDeforLinearIMH8, assembler::A,
    geom::NodalField{FFlt},
    u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    npts0,  Ns0,  gradNparams0,  w0,  pc0 = centroidintegrationdata(self.integdomain);
    imNs, imgradNparams = imintegrationdata(self.nmodes, self.integdomain.integration_rule);
    ecoords, dofnums, loc, J, csmatTJ, loc0, J0, csmatTJ0, gradN, gradN0, imgradN, D, B, DB, Bc, imB, elmatc, elmat = _buffers2(self, geom, u)
    tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)
    startassembly!(assembler, size(elmatc, 1), size(elmatc, 2), count(fes), u.nfreedofs, u.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i]);
        # NOTE: the coordinate system should be evaluated at a single point within the
        # element in order for the derivatives to be consistent at all quadrature points
        # Centroid data
        locjac!(loc0, J0, ecoords, Ns0[1], gradNparams0[1])
        updatecsmat!(self.mcsys, loc0, J0, fes.label[i]);
        Jac0 = Jacobianvolume(self.integdomain, J0, loc0, fes.conn[i], Ns0[1]);
        At_mul_B!(csmatTJ0, csmat(self.mcsys), J0); # local Jacobian matrix
        gradN!(fes, gradN0, gradNparams0[1], csmatTJ0);
        fill!(elmat,  0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
        	locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
        	Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
        	At_mul_B!(csmatTJ, csmat(self.mcsys), J); # local Jacobian matrix
        	gradN!(fes, gradN, gradNparams[j], csmatTJ);
        	Blmat!(self.mr, Bc, Ns[j], gradN, loc, csmat(self.mcsys));
        	B[:, 1:24] .= sqrt(Jac * w[j]) .* Bc
            gradN!(fes, imgradN, imgradNparams[j], csmatTJ0);
        	imBlmat!(self.mr, imB, imNs[j], imgradN, loc0, csmat(self.mcsys), self.nmodes);
        	B[:, 25:end] .= sqrt(Jac0 * w[j]) .* imB
        	add_btdb_ut_only!(elmat, B, 1.0, D, DB)
        end # Loop over quadrature points
        complete_lt!(elmat)
        # Static condensation
        elmatc .= elmat[1:24, 1:24] - elmat[1:24, 25:end] * (elmat[25:end, 25:end] \ elmat[25:end, 1:24])
        gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, elmatc, dofnums, dofnums); # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end


# """
# nzebcloadsstiffness(self::AbstractFEMMDeforLinear,  assembler::A,
#   geom::NodalField{FFlt},
#   u::NodalField{T}) where {A<:AbstractSysvecAssembler, T<:Number}

# Compute load vector for nonzero EBC for fixed displacement.
# """
# function nzebcloadsstiffness(self::FEMMDeforLinearIMH8,  assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysvecAssembler, T<:Number}
#     fes = self.integdomain.fes
#     npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
#     ecoords, dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac, D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix = _buffers2(self, geom, u, npts)
#     realmat = self.material
#     stabmat = self.stabilization_material
#     tangentmoduli!(realmat, D, 0.0, 0.0, loc, 0)
#     tangentmoduli!(stabmat, Dstab, 0.0, 0.0, loc, 0)
#     startassembly!(assembler,  u.nfreedofs);
#     for i = 1:count(fes) # Loop over elements
#         gatherfixedvalues_asvec!(u, elvecfix, fes.conn[i]);# retrieve element displacement vector
#         if norm(elvecfix, Inf) != 0.0   # Is the load nonzero?
#             gathervalues_asmat!(geom, ecoords, fes.conn[i]);
#             # NOTE: the coordinate system should be evaluated at a single point within the
#             # element in order for the derivatives to be consistent at all quadrature points
#             loc = centroid!(self,  loc, ecoords)
#             updatecsmat!(self.mcsys, loc, J, fes.label[i]);
#             vol = 0.0; # volume of the element
#             fill!(MeangradN, 0.0) # mean basis function gradients
#             for j = 1:npts # Loop over quadrature points
#                 jac!(J, ecoords, gradNparams[j])
#                 Jac[j] = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
#                 At_mul_B!(csmatTJ, csmat(self.mcsys), J); # local Jacobian matrix
#                 gradN!(fes, AllgradN[j], gradNparams[j], csmatTJ);
#                 dvol = Jac[j]*w[j]
#                 MeangradN .= MeangradN .+ AllgradN[j]*dvol
#                 vol = vol + dvol
#             end # Loop over quadrature points
#             MeangradN .= MeangradN/vol
#             Blmat!(self.mr, Bbar, Ns[1], MeangradN, loc, csmat(self.mcsys));
#             fill!(elmat,  0.0); # Initialize element matrix
#             add_btdb_ut_only!(elmat, Bbar, vol, D, DB)
#             add_btdb_ut_only!(elmat, Bbar, -self.phis[i]*vol, Dstab, DB)
#             for j = 1:npts # Loop over quadrature points
#                 Blmat!(self.mr, B, Ns[j], AllgradN[j], loc, csmat(self.mcsys));
#                 add_btdb_ut_only!(elmat, B, self.phis[i]*Jac[j]*w[j], Dstab, DB)
#             end # Loop over quadrature points
#             complete_lt!(elmat)
#             gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
#             A_mul_B!(elvec, elmat, elvecfix)
#             assemble!(assembler,  -elvec,  dofnums); # assemble element load vector
#         end
#     end # Loop over elements
#     return makevector!(assembler);
# end

# function nzebcloadsstiffness(self::FEMMDeforLinearIMH8, geom::NodalField{FFlt}, u::NodalField{T}) where {T<:Number}
#     assembler = SysvecAssembler()
#     return  nzebcloadsstiffness(self, assembler, geom, u);
# end

end
