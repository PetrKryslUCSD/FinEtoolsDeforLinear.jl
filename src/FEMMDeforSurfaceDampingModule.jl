"""
    FEMMDeforSurfaceDampingModule

Module for operations on the damping associated with absorbing boundary 
conditions (ABC) representation of the effect of infinite extent 
of inviscid fluid next to the surface.
"""
module FEMMDeforSurfaceDampingModule

__precompile__(true)

using FinEtools.FENodeSetModule: FENodeSet
using FinEtools.FESetModule: AbstractFESet, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobiansurface
using FinEtools.FieldModule: ndofs, gatherdofnums!, gathervalues_asmat!, nalldofs
using FinEtools.NodalFieldModule: NodalField
using FinEtools.FEMMBaseModule: AbstractFEMM
using FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!
using FinEtools.MatrixUtilityModule: add_nnt_ut_only!, complete_lt!, locjac!
using FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
using LinearAlgebra: norm, cross

"""
    FEMMDeforSurfaceDamping{S<:AbstractFESet, F<:Function} <: AbstractFEMM

Type for surface damping model.
"""
mutable struct FEMMDeforSurfaceDamping{ID<:IntegDomain} <: AbstractFEMM
    integdomain::ID # geometry data
end

"""
    dampingABC(self::FEMMDeforSurfaceDamping, assembler::A, geom::NodalField{GFT}, u::NodalField{UFT}, impedance::UFT, surfacenormal::SurfaceNormal) where {A<:AbstractSysmatAssembler, GFT<:Number, UFT<:Number}

Compute the damping matrix associated with absorbing boundary conditions (ABC).

Compute the damping matrix associated with absorbing boundary conditions (ABC)
representation of the effect of infinite extent of inviscid fluid next to
the surface.
"""
function dampingABC(self::FEMMDeforSurfaceDamping, assembler::A, geom::NodalField{GFT}, u::NodalField{UFT}, impedance::FT, surfacenormal::SurfaceNormal) where {A<:AbstractSysmatAssembler, GFT<:Number, UFT<:Number, FT<:Number}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes per finite element
    sdim = ndofs(geom); # spatial dimension
    mdim = manifdim(fes); # manifold dimension of the finite elements
    Cedim = ndn * nne; # size of damping element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    loc = zeros(GFT, 1, sdim); # quadrature point coordinate -- used as a buffer
    J = zeros(GFT, sdim, mdim); # Jacobian matrix -- used as a buffer
    ecoords = zeros(GFT, nne, sdim)
    Ce = zeros(FT, Cedim, Cedim); # element damping matrix -- used as a buffer
    Nn = zeros(GFT, Cedim); # column vector
    dofnums = zeros(eltype(u.dofnums), Cedim); # degree of freedom array -- used as a buffer
    # Prepare assembler and temporaries
    startassembly!(assembler, Cedim^2 * nfes, nalldofs(u), nalldofs(u))
    for i  in eachindex(fes)  # loop over finite elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i]);
        fill!(Ce, 0.0); # Initialize element damping matrix
        for j  in  1:npts # loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            n = updatenormal!(surfacenormal, loc, J, i, j)
            for k = 1:nne
                Nn[(k-1)*ndn+1:k*ndn] = n * Ns[j][k]
            end
            add_nnt_ut_only!(Ce, Nn, (+1)*impedance*Jac*w[j]);
        end # end loop over quadrature points
        complete_lt!(Ce);
        gatherdofnums!(u, dofnums, self.integdomain.fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, Ce, dofnums, dofnums); # assemble the element damping matrix
    end # end loop over finite elements
    return makematrix!(assembler);
end

function dampingABC(self::FEMMDeforSurfaceDamping, geom::NodalField{GFT}, u::NodalField{UFT}, impedance::FT, surfacenormal::SurfaceNormal) where {GFT<:Number, UFT<:Number, FT<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return dampingABC(self, assembler, geom, u, impedance, surfacenormal);
end

end
