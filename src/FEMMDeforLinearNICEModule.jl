"""
    FEMMDeforLinearNICEModule

Formulation for the small displacement, small strain deformation
model for Nodally-Integrated Continuum Elements (NICE).

The approximation is  originally from Dohrmann et al IJNME 47 (2000).
The formulation was subsequently developed in Krysl, P. and Zhu, B.
Locking-free continuum displacement finite elements with nodal
integration, International Journal for Numerical Methods in Engineering,
76,7,1020-1043,2008.
"""
module FEMMDeforLinearNICEModule

__precompile__(true)

using FinEtools.FENodeSetModule: FENodeSet
using FinEtools.FESetModule: AbstractFESet, FESetH8, FESetT4, manifdim, nodesperelem, gradN!
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: AbstractFEMMDeforLinear
using FinEtools.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D
using FinEtoolsDeforLinear.MatDeforLinearElasticModule:
    AbstractMatDeforLinearElastic, tangentmoduli!, update!, thermalstrain!
using FinEtools.FieldModule:
    ndofs,
    gatherdofnums!,
    gatherfixedvalues_asvec!,
    gathervalues_asvec!,
    gathervalues_asmat!,
    nalldofs
using FinEtools.NodalFieldModule: NodalField, nnodes
using FinEtools.CSysModule: CSys, updatecsmat!, csmat
using FinEtools.FENodeToFEMapModule: FENodeToFEMap
using FinEtools.DeforModelRedModule: nstressstrain, nthermstrain, blmat!
using FinEtools.AssemblyModule:
    AbstractSysvecAssembler,
    AbstractSysmatAssembler,
    SysmatAssemblerSparseSymm,
    startassembly!,
    assemble!,
    makematrix!,
    makevector!,
    SysvecAssembler
using FinEtools.MatrixUtilityModule:
    add_btdb_ut_only!, complete_lt!, loc!, jac!, locjac!, adjugate3!
import FinEtoolsDeforLinear.FEMMDeforLinearBaseModule:
    stiffness, mass, thermalstrainloads, inspectintegpoints
import FinEtools.FEMMBaseModule: associategeometry!
using FinEtoolsDeforLinear.MatDeforModule: rotstressvec!
using LinearAlgebra: mul!, Transpose, UpperTriangular, eigvals
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
using LinearAlgebra: norm, qr, diag, dot, cond, I
using Statistics: mean

const StabParamFloat = Float64

"""
    AbstractFEMMDeforLinearNICE <: AbstractFEMMDeforLinear

Abstract FEMM type for Nodally Integrated Continuum Elements (NICE).
"""
abstract type AbstractFEMMDeforLinearNICE <: AbstractFEMMDeforLinear end

mutable struct _NodalBasisFunctionGradients{FT,IT}
    gradN::Matrix{FT}
    patchconn::Vector{IT}
    Vpatch::FT
end

"""
    mutable struct FEMMDeforLinearNICEH8{
        MR <: AbstractDeforModelRed,
        S <: FESetH8,
        F <: Function,
        M <: AbstractMatDeforLinearElastic,
    } <: AbstractFEMMDeforLinearNICE

FEMM type for Nodally Integrated Continuum Elements (NICE) based on the eight-node hexahedron.
"""
mutable struct FEMMDeforLinearNICEH8{
    MR<:AbstractDeforModelRed,
    ID<:IntegDomain{S,F} where {S<:FESetH8,F<:Function},
    CS<:CSys,
    M<:AbstractMatDeforLinearElastic,
} <: AbstractFEMMDeforLinearNICE
    mr::Type{MR}
    integdomain::ID # integration domain
    mcsys::CS # updater of the material orientation matrix
    material::M # material object
    stabfact::StabParamFloat
    nodalbasisfunctiongrad::Vector{_NodalBasisFunctionGradients}
end

function FEMMDeforLinearNICEH8(
    mr::Type{MR},
    integdomain::IntegDomain{S,F},
    mcsys::CSys,
    material::M,
) where {MR<:AbstractDeforModelRed,S<:FESetH8,F<:Function,M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabfact = 0.05
    return FEMMDeforLinearNICEH8(
        mr,
        integdomain,
        mcsys,
        material,
        stabfact,
        _NodalBasisFunctionGradients[],
    )
end

function FEMMDeforLinearNICEH8(
    mr::Type{MR},
    integdomain::IntegDomain{S,F},
    material::M,
) where {MR<:AbstractDeforModelRed,S<:FESetH8,F<:Function,M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabfact = 0.05
    return FEMMDeforLinearNICEH8(
        mr,
        integdomain,
        CSys(manifdim(integdomain.fes)),
        material,
        stabfact,
        _NodalBasisFunctionGradients[],
    )
end

function FEMMDeforLinearNICEH8(
    mr::Type{MR},
    integdomain::IntegDomain{S,F},
    material::M,
    stabfact,
) where {MR<:AbstractDeforModelRed,S<:FESetH8,F<:Function,M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    return FEMMDeforLinearNICEH8(
        mr,
        integdomain,
        CSys(manifdim(integdomain.fes)),
        material,
        stabfact,
        _NodalBasisFunctionGradients[],
    )
end

"""
    mutable struct FEMMDeforLinearNICET4{
        MR <: AbstractDeforModelRed,
        S <: FESetT4,
        F <: Function,
        M <: AbstractMatDeforLinearElastic,
    } <: AbstractFEMMDeforLinearNICE

FEMM type for Nodally Integrated Continuum Elements (NICE) based on the 4-node tetrahedron.
"""
mutable struct FEMMDeforLinearNICET4{
    MR<:AbstractDeforModelRed,
    ID<:IntegDomain{S,F} where {S<:FESetT4,F<:Function},
    CS<:CSys,
    M<:AbstractMatDeforLinearElastic,
} <: AbstractFEMMDeforLinearNICE
    mr::Type{MR}
    integdomain::ID # geometry data
    mcsys::CS # updater of the material orientation matrix
    material::M # material object
    stabfact::StabParamFloat
    nodalbasisfunctiongrad::Vector{_NodalBasisFunctionGradients}
end

function FEMMDeforLinearNICET4(
    mr::Type{MR},
    integdomain::IntegDomain{S,F},
    mcsys::CSys,
    material::M,
) where {MR<:AbstractDeforModelRed,S<:FESetT4,F<:Function,M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabfact = 0.015
    return FEMMDeforLinearNICET4(
        mr,
        integdomain,
        mcsys,
        material,
        stabfact,
        _NodalBasisFunctionGradients[],
    )
end

function FEMMDeforLinearNICET4(
    mr::Type{MR},
    integdomain::IntegDomain{S,F},
    material::M,
) where {MR<:AbstractDeforModelRed,S<:FESetT4,F<:Function,M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabfact = 0.015
    return FEMMDeforLinearNICET4(
        mr,
        integdomain,
        CSys(manifdim(integdomain.fes)),
        material,
        stabfact,
        _NodalBasisFunctionGradients[],
    )
end

function FEMMDeforLinearNICET4(
    mr::Type{MR},
    integdomain::IntegDomain{S,F},
    material::M,
    stabfact,
) where {MR<:AbstractDeforModelRed,S<:FESetT4,F<:Function,M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    return FEMMDeforLinearNICET4(
        mr,
        integdomain,
        CSys(manifdim(integdomain.fes)),
        material,
        stabfact,
        _NodalBasisFunctionGradients[],
    )
end

function _buffers1(
    self::AbstractFEMMDeforLinearNICE,
    geom::NodalField{GFT},
    npts::Int,
) where {GFT}
    fes = self.integdomain.fes
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    # Prepare buffers
    loc = fill(zero(GFT), 1, sdim) # quadrature point location -- buffer
    J = fill(zero(GFT), sdim, mdim) # Jacobian matrix -- buffer
    adjJ = fill(zero(GFT), sdim, mdim) # Jacobian matrix -- buffer
    csmatTJ = fill(zero(GFT), mdim, mdim) # intermediate result -- buffer
    gradN = fill(zero(GFT), nne, mdim)
    xl = fill(zero(GFT), nne, mdim)
    return loc, J, adjJ, csmatTJ, gradN, xl
end

function _buffers2(
    self::AbstractFEMMDeforLinearNICE,
    geom::NodalField{GFT},
    u::NodalField,
    npts::Int,
) where {GFT}
    fes = self.integdomain.fes
    ndn = ndofs(u) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    nstrs = nstressstrain(self.mr)  # number of stresses
    loc = fill(zero(GFT), 1, sdim) # quadrature point location -- buffer
    J = fill(zero(GFT), sdim, mdim) # Jacobian matrix -- buffer
    csmatTJ = fill(zero(GFT), mdim, mdim) # intermediate result -- buffer
    Jac = fill(zero(GFT), npts)
    D = fill(zero(GFT), nstrs, nstrs) # material stiffness matrix -- buffer
    return loc, J, csmatTJ, Jac, D
end

function patchconn(fes, gl, thisnn)
    # Generate patch connectivity for a given node (thisnn)
    # from the connectivities of the finite elements attached to it.
    return vcat(
        collect(setdiff(Set([i for j in eachindex(gl) for i in fes.conn[gl[j]]]), thisnn)),
        [thisnn],
    )
end

function computenodalbfungrads(self, geom)
    # # Compute the nodal basis function gradients.
    # # Return the cell array of structures with attributes
    # %		 bfun_gradients{nix}.Nspd= basis function gradient matrix
    # #        bfun_gradients{nix}.Vpatch= nodal patch volume
    # #        bfun_gradients{nix}.patchconn= nodal patch connectivity

    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    loc, J, adjJ, csmatTJ, gradN, xl = _buffers1(self, geom, npts)

    # Get the inverse map from finite element nodes to geometric cells
    fen2fe = FENodeToFEMap(fes.conn, nnodes(geom))
    # Initialize the nodal gradients, nodal patch, and patch connectivity
    bfungrads =
        fill(_NodalBasisFunctionGradients(fill(0.0, 0, 0), fill(0, 0), 0.0), nnodes(geom))
    # Now loop over all finite element nodes in the map
    lnmap = fill(0, length(fen2fe.map)) # Local node map: buffer to speed up operations
    for nix in eachindex(fen2fe.map)
        gl = fen2fe.map[nix]
        thisnn = nix # We are at this node
        if !isempty(gl) # This node has an element patch in this block
            # establish local numbering of all nodes of the patch @ node thisnn
            p = patchconn(fes, gl, thisnn)
            np = length(p)
            lnmap[p] .= 1:np# now store the local numbers
            c = reshape(geom.values[thisnn, :], 1, ndofs(geom))
            updatecsmat!(self.mcsys, c, J, nix, 0)
            gradNavg = fill(0.0, np, ndofs(geom))# preallocate strain-displacement matrix
            Vpatch = 0.0
            for k in eachindex(gl)
                i = gl[k]
                kconn = collect(fes.conn[i])
                pci = findfirst(cx -> cx == thisnn, kconn)# at which node in the element are we with this quadrature point?
                @assert 1 <= pci <= nodesperelem(fes)
                # centered coordinates of nodes in the material coordinate system
                for cn in eachindex(kconn)
                    xl[cn, :] =
                        (reshape(geom.values[kconn[cn], :], 1, ndofs(geom)) - c) *
                        csmat(self.mcsys)
                end
                jac!(J, xl, gradNparams[pci])
                At_mul_B!(csmatTJ, csmat(self.mcsys), J) # local Jacobian matrix
                Jac = Jacobianvolume(self.integdomain, J, c, fes.conn[i], Ns[pci])
                Vpatch += Jac * w[pci]
                sgradN = gradNparams[pci] * adjugate3!(adjJ, J)
                gradNavg[lnmap[kconn], :] += (w[pci] .* sgradN)
            end
            @assert Vpatch != 0
            gradNavg ./= Vpatch
            bfungrads[nix] = _NodalBasisFunctionGradients(gradNavg, p, Vpatch)
            lnmap[p] .= 0 # Restore the buffer to pristine condition
        end
    end
    self.nodalbasisfunctiongrad = bfungrads
    return self
end

"""
    associategeometry!(self::F,
        geom::NodalField{GFT}) where {F <: AbstractFEMMDeforLinearNICE, GFT}

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(
    self::F,
    geom::NodalField{GFT},
) where {F<:AbstractFEMMDeforLinearNICE,GFT}
    return computenodalbfungrads(self, geom)
end

function Phi3(dim, np, lx, c)
    lx[:, 1] = lx[:, 1] .- c[1]
    lx[:, 2] = lx[:, 2] .- c[2]
    lx[:, 3] = lx[:, 3] .- c[3]
    Phi = fill(0.0, dim * np, 12)
    Phi[1:dim:(end-1), 1] .= 1
    Phi[1:dim:(end-1), 2] = lx[:, 1]'
    Phi[1:dim:(end-1), 3] = lx[:, 2]'
    Phi[1:dim:(end-1), 4] = lx[:, 3]'
    Phi[2:dim:end, 5] .= 1
    Phi[2:dim:end, 6] = lx[:, 1]'
    Phi[2:dim:end, 7] = lx[:, 2]'
    Phi[2:dim:end, 8] = lx[:, 3]'
    Phi[3:dim:end, 9] .= 1
    Phi[3:dim:end, 10] = lx[:, 1]'
    Phi[3:dim:end, 11] = lx[:, 2]'
    Phi[3:dim:end, 12] = lx[:, 3]'
    return Phi
end

"""
    stiffness(self::AbstractFEMMDeforLinearNICE,
        assembler::A,
        geom::NodalField{GFT},
        u::NodalField{UFT}) where {A <: AbstractSysmatAssembler, GFT <: Number, UFT <: Number}

Compute and assemble  stiffness matrix.
"""
function stiffness(
    self::AbstractFEMMDeforLinearNICE,
    assembler::A,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {A<:AbstractSysmatAssembler,GFT<:Number,UFT<:Number}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    loc, J, csmatTJ, Jac, D = _buffers2(self, geom, u, npts)
    tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)
    Dmod = sort(eigvals(D))
    stabDmod = mean(Dmod[1:2]; dims = 1)
    elmatsizeguess = 4 * nodesperelem(fes) * ndofs(u)
    startassembly!(assembler, elmatsizeguess, elmatsizeguess, nnodes(u), nalldofs(u), nalldofs(u))
    for nix in eachindex(self.nodalbasisfunctiongrad)
        gradN = self.nodalbasisfunctiongrad[nix].gradN
        patchconn = self.nodalbasisfunctiongrad[nix].patchconn
        Vpatch = self.nodalbasisfunctiongrad[nix].Vpatch
        c = reshape(geom.values[nix, :], 1, ndofs(geom))
        updatecsmat!(self.mcsys, c, J, nix, 0)
        nd = length(patchconn) * ndofs(u)
        Bnodal = fill(0.0, size(D, 1), nd)
        blmat!(self.mr, Bnodal, Ns[1], gradN, c, csmat(self.mcsys))
        elmat = fill(0.0, nd, nd) # Can we SPEED it UP?
        DB = fill(0.0, size(D, 1), nd)
        add_btdb_ut_only!(elmat, Bnodal, Vpatch, D, DB)
        complete_lt!(elmat)
        if (self.stabfact > 0)
            Phi = Phi3(ndofs(u), length(patchconn), geom.values[patchconn, :], c)
            A1 = Phi * ((Phi' * Phi) \ Phi')
            elmat += (self.stabfact * stabDmod) .* (I - A1)
        end
        dofnums = fill(0, nd)
        gatherdofnums!(u, dofnums, patchconn) # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums) # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

"""
    stiffness(self::AbstractFEMMDeforLinearNICE,
        geom::NodalField{GFT},
        u::NodalField{UFT}) where {GFT <: Number, UFT <: Number}

Compute and assemble  stiffness matrix.
"""
function stiffness(
    self::AbstractFEMMDeforLinearNICE,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {GFT<:Number,UFT<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return stiffness(self, assembler, geom, u)
end

end
