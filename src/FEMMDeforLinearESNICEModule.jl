"""
Formulation for the small displacement, small strain deformation
model for Nodally-Integrated Continuum Elements (NICE).

The approximation is  originally from Dohrmann et al IJNME 47 (2000).
The formulation was subsequently developed in Krysl, P. and Zhu, B.
Locking-free continuum displacement finite elements with nodal
integration, International Journal for Numerical Methods in Engineering,
76,7,1020-1043,2008.

The stabilization scheme comes from papers on energy-sampling stabilization
for mean-strain elements (Krysl and coauthors).
"""
module FEMMDeforLinearESNICEModule

__precompile__(true)

using FinEtools.FENodeSetModule: FENodeSet
using FinEtools.FESetModule: AbstractFESet, FESetH8, FESetT4, manifdim, nodesperelem, gradN!
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: AbstractFEMMDeforLinear
using FinEtools.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D
using FinEtoolsDeforLinear.MatDeforLinearElasticModule:
    AbstractMatDeforLinearElastic, tangentmoduli!, update!, thermalstrain!
using FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso
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
using FinEtools.DeforModelRedModule: nstressstrain, nthermstrain, blmat!, divmat, vgradmat
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
using FinEtools.MatModule: massdensity
using LinearAlgebra: mul!, Transpose, UpperTriangular, eigvals, det
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
using LinearAlgebra: norm, qr, diag, dot, cond, I, cross
using Statistics: mean

const StabParamFloat = Float64

"""
    AbstractFEMMDeforLinearESNICE <: AbstractFEMMDeforLinear

Abstract FEMM type for Nodally Integrated Continuum Elements (NICE) with
energy-sampling stabilization (ESNICE).
"""
abstract type AbstractFEMMDeforLinearESNICE <: AbstractFEMMDeforLinear end

# Tetrahedron
# The coefficient set below was obtained by fitting the ratio of energies
# true/approximate for the finite element model of six tetrahedra arranged
# into a rectangular block and subject to pure bending
# Fitting for a small aspect-ratio range (1.0 to 10)
_T4_stabilization_parameters = (2.101588423297799, 1.311321055432958)

mutable struct _NodalBasisFunctionGradients{FT,IT}
    gradN::Matrix{FT}
    patchconn::Vector{IT}
    Vpatch::FT
end

function _make_stabilization_material(material::M) where {M}
    ns = fieldnames(typeof(material))
    E = 0.0
    nu = 0.0
    if :E in ns
        E = material.E
        if material.nu < 0.3
            nu = material.nu
        else
            nu = 0.3 + (material.nu - 0.3) / 2.0
        end
    else
        if :E1 in ns
            E = mean([material.E1, material.E2, material.E3])
            nu = min(material.nu12, material.nu13, material.nu23)
        else
            error("No clues on how to construct the stabilization material")
        end
    end
    return MatDeforElastIso(material.mr, 0.0, E, nu, 0.0)
end

"""
    mutable struct FEMMDeforLinearESNICET4{
        MR<:AbstractDeforModelRed,
        ID<:IntegDomain{S} where {S<:FESetT4},
        CS<:CSys,
        M<:AbstractMatDeforLinearElastic,
        MS<:MatDeforElastIso,
    } <: AbstractFEMMDeforLinearESNICE

FEMM type for Energy-sampling Stabilized Nodally Integrated Continuum Elements
(ESNICE) based on the 4-node tetrahedron.
"""
mutable struct FEMMDeforLinearESNICET4{
    MR<:AbstractDeforModelRed,
    ID<:IntegDomain{S} where {S<:FESetT4},
    CS<:CSys,
    M<:AbstractMatDeforLinearElastic,
    MS<:MatDeforElastIso,
} <: AbstractFEMMDeforLinearESNICE
    mr::Type{MR}
    integdomain::ID # geometry data
    mcsys::CS # updater of the material orientation matrix
    material::M # material object
    stabilization_material::MS
    nodalbasisfunctiongrad::Vector{_NodalBasisFunctionGradients}
    ephis::Vector{StabParamFloat}
    nphis::Vector{StabParamFloat}
end

"""
    mutable struct FEMMDeforLinearESNICEH8{
        MR<:AbstractDeforModelRed,
        ID<:IntegDomain{S} where {S<:FESetH8},
        CS<:CSys,
        M<:AbstractMatDeforLinearElastic,
        MS<:MatDeforElastIso,
    } <: AbstractFEMMDeforLinearESNICE

FEMM type for Energy-sampling Stabilized Nodally Integrated Continuum Elements
(ESNICE) based on the 8-node hexahedron.
"""
mutable struct FEMMDeforLinearESNICEH8{
    MR<:AbstractDeforModelRed,
    ID<:IntegDomain{S} where {S<:FESetH8},
    CS<:CSys,
    M<:AbstractMatDeforLinearElastic,
    MS<:MatDeforElastIso,
} <: AbstractFEMMDeforLinearESNICE
    mr::Type{MR}
    integdomain::ID # geometry data
    mcsys::CS # updater of the material orientation matrix
    material::M # material object
    stabilization_material::MS
    nodalbasisfunctiongrad::Vector{_NodalBasisFunctionGradients}
    ephis::Vector{StabParamFloat}
    nphis::Vector{StabParamFloat}
end

"""
    FEMMDeforLinearESNICET4(
        mr::Type{MR},
        integdomain::ID,
        mcsys::CS,
        material::M,
    ) where {
        MR<:AbstractDeforModelRed,
        ID<:IntegDomain{S} where {S<:FESetT4},
        CS<:CSys,
        M<:AbstractMatDeforLinearElastic,
    }

Constructor.
"""
function FEMMDeforLinearESNICET4(
    mr::Type{MR},
    integdomain::ID,
    mcsys::CS,
    material::M,
) where {
    MR<:AbstractDeforModelRed,
    ID<:IntegDomain{S} where {S<:FESetT4},
    CS<:CSys,
    M<:AbstractMatDeforLinearElastic,
}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearESNICET4(
        mr,
        integdomain,
        mcsys,
        material,
        stabilization_material,
        _NodalBasisFunctionGradients[],
        fill(zero(StabParamFloat), 1),
        fill(zero(StabParamFloat), 1),
    )
end

"""
    FEMMDeforLinearESNICET4(
        mr::Type{MR},
        integdomain::ID,
        material::M,
    ) where {
        MR<:AbstractDeforModelRed,
        ID<:IntegDomain{S} where {S<:FESetT4},
        M<:AbstractMatDeforLinearElastic,
    }

Constructor.
"""
function FEMMDeforLinearESNICET4(
    mr::Type{MR},
    integdomain::ID,
    material::M,
) where {
    MR<:AbstractDeforModelRed,
    ID<:IntegDomain{S} where {S<:FESetT4},
    M<:AbstractMatDeforLinearElastic,
}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearESNICET4(
        mr,
        integdomain,
        CSys(manifdim(integdomain.fes)),
        material,
        stabilization_material,
        _NodalBasisFunctionGradients[],
        fill(zero(StabParamFloat), 1),
        fill(zero(StabParamFloat), 1),
    )
end

"""
    FEMMDeforLinearESNICEH8(
        mr::Type{MR},
        integdomain::ID,
        mcsys::CS,
        material::M,
    ) where {
        MR<:AbstractDeforModelRed,
        ID<:IntegDomain{S} where {S<:FESetH8},
        CS<:CSys,
        M<:AbstractMatDeforLinearElastic,
    }

Constructor.
"""
function FEMMDeforLinearESNICEH8(
    mr::Type{MR},
    integdomain::ID,
    mcsys::CS,
    material::M,
) where {
    MR<:AbstractDeforModelRed,
    ID<:IntegDomain{S} where {S<:FESetH8},
    CS<:CSys,
    M<:AbstractMatDeforLinearElastic,
}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearESNICEH8(
        mr,
        integdomain,
        mcsys,
        material,
        stabilization_material,
        _NodalBasisFunctionGradients[],
        fill(zero(StabParamFloat), 1),
        fill(zero(StabParamFloat), 1),
    )
end

"""
    FEMMDeforLinearESNICEH8(
        mr::Type{MR},
        integdomain::ID,
        material::M,
    ) where {
        MR<:AbstractDeforModelRed,
        ID<:IntegDomain{S} where {S<:FESetH8},
        M<:AbstractMatDeforLinearElastic,
    }

Constructor.
"""
function FEMMDeforLinearESNICEH8(
    mr::Type{MR},
    integdomain::ID,
    material::M,
) where {
    MR<:AbstractDeforModelRed,
    ID<:IntegDomain{S} where {S<:FESetH8},
    M<:AbstractMatDeforLinearElastic,
}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearESNICEH8(
        mr,
        integdomain,
        CSys(manifdim(integdomain.fes)),
        material,
        stabilization_material,
        _NodalBasisFunctionGradients[],
        fill(zero(StabParamFloat), 1),
        fill(zero(StabParamFloat), 1),
    )
end

function _buffers1(self::FEMM, geom::NodalField{GFT}) where {FEMM<:AbstractFEMMDeforLinearESNICE,GFT}
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
    self::FEMM,
    geom::NodalField{GFT},
    u::NodalField,
    npts::Int,
) where {FEMM<:AbstractFEMMDeforLinearESNICE,GFT}
    fes = self.integdomain.fes
    ndn = ndofs(u) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    nstrs = nstressstrain(self.mr)  # number of stresses
    elmatdim = ndn * nne             # dimension of the element matrix
    # Prepare buffers
    ecoords = fill(zero(GFT), nne, ndofs(geom)) # array of Element coordinates
    dofnums = zeros(eltype(u.dofnums), elmatdim) # degree of freedom array -- buffer
    loc = fill(zero(GFT), 1, sdim) # quadrature point location -- buffer
    J = fill(zero(GFT), sdim, mdim) # Jacobian matrix -- buffer
    csmatTJ = fill(zero(GFT), mdim, mdim) # intermediate result -- buffer
    Jac = fill(zero(GFT), npts)
    D = fill(zero(GFT), nstrs, nstrs) # material stiffness matrix -- buffer
    Dstab = fill(zero(GFT), nstrs, nstrs) # material stiffness matrix -- buffer
    elmat = fill(zero(GFT), elmatdim, elmatdim)      # element matrix -- buffer
    B = fill(zero(GFT), nstrs, elmatdim) # strain-displacement matrix -- buffer
    return ecoords, dofnums, loc, J, csmatTJ, Jac, D, Dstab, elmat, B
end

function _buffers3(
    self::FEMM,
    geom::NodalField{GFT},
    u::NodalField,
) where {FEMM<:AbstractFEMMDeforLinearESNICE,GFT}
    fes = self.integdomain.fes
    ndn = ndofs(u) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    nstrs = nstressstrain(self.mr)  # number of stresses
    elmatdim = ndn * nne             # dimension of the element matrix
    # Prepare buffers
    elmat = fill(zero(GFT), elmatdim, elmatdim)      # element matrix -- buffer
    dofnums = zeros(eltype(u.dofnums), elmatdim) # degree of freedom array -- buffer
    B = fill(zero(GFT), nstrs, elmatdim) # strain-displacement matrix -- buffer
    DB = fill(zero(GFT), nstrs, elmatdim) # strain-displacement matrix -- buffer
    elvecfix = fill(zero(GFT), elmatdim) # vector of prescribed displ. -- buffer
    elvec = fill(zero(GFT), elmatdim) # element vector -- buffer
    gradN = fill(zero(GFT), nne, mdim) # intermediate result -- buffer
    return dofnums, B, DB, elmat, elvec, elvecfix, gradN
end

function _patchconn(fes, gl, thisnn)
    # Generate patch connectivity for a given node (thisnn)
    # from the connectivities of the finite elements attached to it.
    return vcat(
        collect(setdiff(Set([i for j in eachindex(gl) for i in fes.conn[gl[j]]]), thisnn)),
        [thisnn],
    )
end

function _computenodalbfungrads(self, geom)
    # # Compute the nodal basis function gradients.
    # # Return the cell array of structures with attributes
    # %      bfun_gradients{nix}.Nspd= basis function gradient matrix
    # #        bfun_gradients{nix}.Vpatch= nodal patch volume
    # #        bfun_gradients{nix}.patchconn= nodal patch connectivity

    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    loc, J, adjJ, csmatTJ, gradN, xl = _buffers1(self, geom)

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
            p = _patchconn(fes, gl, thisnn)
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

function _tetaspectratiovol(X)
    edge1 = vec(X[2, :] - X[1, :])
    edge2 = vec(X[3, :] - X[1, :])
    edge3 = vec(X[4, :] - X[1, :])
    edge4 = vec(X[3, :] - X[2, :])
    edge5 = vec(X[4, :] - X[3, :])
    edge6 = vec(X[4, :] - X[2, :])
    V = dot(edge3, cross(edge1, edge2))
    A1 = norm(cross(edge1, edge2)) # This is twice the area of the triangle
    A2 = norm(cross(edge2, edge3))
    A3 = norm(cross(edge3, edge1))
    A4 = norm(cross(edge4, edge6))
    #h1, h2, h3, h4 = V/A1, V/A2, V/A3, V/A4
    L1, L2, L3, L4, L5, L6 =
        norm(edge1), norm(edge2), norm(edge3), norm(edge4), norm(edge5), norm(edge6)
    # the heights and the edge lengths will now be used to derive a composite
    # index: aspect ratio index. If this cannot be computed without generating
    # not-a-number or infinity, this number is assumed to be better represented
    # with 1.0 (perfect ratio), so that the element shape is then governed by
    # the other aspect ratio values.
    f = maximum
    f([L1, L2, L4]) != 0.0 && A1 != 0.0 ? ar1 = V / A1 / f([L1, L2, L4]) : ar1 = 1.0
    f([L3, L2, L5]) != 0.0 && A2 != 0.0 ? ar2 = V / A2 / f([L3, L2, L5]) : ar2 = 1.0
    f([L1, L3, L6]) != 0.0 && A3 != 0.0 ? ar3 = V / A3 / f([L1, L3, L6]) : ar3 = 1.0
    f([L6, L5, L4]) != 0.0 && A4 != 0.0 ? ar4 = V / A4 / f([L6, L5, L4]) : ar4 = 1.0
    return ar1, ar2, ar3, ar4, V / 6
end

"""
    associategeometry!(
        self::FEMM,
        geom::NodalField{GFT};
        stabilization_parameters = _T4_stabilization_parameters,
    ) where {FEMM<:FEMMDeforLinearESNICET4,GFT}

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(
    self::FEMM,
    geom::NodalField{GFT};
    stabilization_parameters = _T4_stabilization_parameters,
) where {FEMM<:FEMMDeforLinearESNICET4,GFT}
    (a, b) = stabilization_parameters
    fes = self.integdomain.fes
    self.ephis = fill(zero(StabParamFloat), count(fes))
    evols = fill(zero(StabParamFloat), count(fes))
    self.nphis = fill(zero(StabParamFloat), nnodes(geom))
    nvols = fill(zero(StabParamFloat), nnodes(geom))
    for i = 1:count(fes) # Loop over elements
        ar1, ar2, ar3, ar4, V = _tetaspectratiovol(geom.values[collect(fes.conn[i]), :])
        evols[i] = V
        # If the aspect ratios are not reasonable, such as when the element is a
        # sliver or inverted, we turn off the stabilization for the element
        # by setting its stabilization factor to zero.
        if min(ar1, ar2, ar3, ar4) <= 0
            self.ephis[i] = 0.0
        else
            self.ephis[i] = (1.0 / (b * min(ar1, ar2, ar3, ar4)^a) + 1.0)^(-1)
        end
        # Accumulate: the stabilization factor at the node is the weighted mean
        # of the stabilization factors of the elements at that node
        for k = 1:nodesperelem(fes)
            nvols[fes.conn[i][k]] += evols[i]
            self.nphis[fes.conn[i][k]] += self.ephis[i] * evols[i]
        end
    end # Loop over elements
    @assert any(isnan.(self.ephis)) == false
    # Now scale the values at the nodes with the nodal volumes
    for k in eachindex(nvols)
        if nvols[k] != 0.0
            self.nphis[k] /= nvols[k]
        end
    end
    @assert any(isnan.(self.nphis)) == false
    # Now calculate the nodal basis function gradients
    return _computenodalbfungrads(self, geom)
end

"""
    associategeometry!(
        self::FEMM,
        geom::NodalField{GFT},
    ) where {FEMM<:FEMMDeforLinearESNICEH8,GFT}

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(
    self::FEMM,
    geom::NodalField{GFT},
) where {FEMM<:FEMMDeforLinearESNICEH8,GFT}
    fes = self.integdomain.fes
    self.ephis = fill(zero(StabParamFloat), count(fes))
    evols = fill(zero(StabParamFloat), count(fes))
    self.nphis = fill(zero(StabParamFloat), nnodes(geom))
    nvols = fill(zero(StabParamFloat), nnodes(geom))
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)

    for i = 1:count(fes) # Loop over elements
        X = geom.values[collect(fes.conn[i]), :]
        V = 0.0
        for j = 1:npts
            J = X' * gradNparams[j]
            Jac = det(J)
            @assert Jac > 0 "Nonpositive Jacobian"
            V = V + Jac * w[j]
            hs = sum(J .* J; dims = 1)
            phi = 2 * (1.0 + self.stabilization_material.nu) * minimum(hs) / maximum(hs)
            self.ephis[i] = max(self.ephis[i], phi / (1 + phi))
        end
        evols[i] = V
        # Accumulate: the stabilization factor at the node is the weighted mean of the stabilization factors of the elements at that node
        for k = 1:nodesperelem(fes)
            nvols[fes.conn[i][k]] += evols[i]
            self.nphis[fes.conn[i][k]] += self.ephis[i] * evols[i]
        end
    end # Loop over elements
    # Now scale the values at the nodes with the nodal volumes
    for k in eachindex(nvols)
        self.nphis[k] /= nvols[k]
    end
    # Now calculate the nodal basis function gradients
    return _computenodalbfungrads(self, geom)
end

"""
    stiffness(
        self::FEMM,
        assembler::A,
        geom::NodalField{GFT},
        u::NodalField{T},
    ) where {FEMM<:AbstractFEMMDeforLinearESNICE,A<:AbstractSysmatAssembler,GFT<:Number,T<:Number}

Compute and assemble  stiffness matrix.
"""
function stiffness(
    self::FEMM,
    assembler::A,
    geom::NodalField{GFT},
    u::NodalField{T},
) where {FEMM<:AbstractFEMMDeforLinearESNICE,A<:AbstractSysmatAssembler,GFT<:Number,T<:Number}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    ecoords, dofnums, loc, J, csmatTJ, Jac, D, Dstab = _buffers2(self, geom, u, npts)
    tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)
    tangentmoduli!(self.stabilization_material, Dstab, 0.0, 0.0, loc, 0)
    elmatsizeguess = 4 * nodesperelem(fes) * ndofs(u)
    startassembly!(assembler, elmatsizeguess, elmatsizeguess, nnodes(u), nalldofs(u), nalldofs(u))
    for nix = 1:length(self.nodalbasisfunctiongrad)
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
        add_btdb_ut_only!(elmat, Bnodal, -self.nphis[nix] * Vpatch, Dstab, DB)
        complete_lt!(elmat)
        dofnums = fill(0, nd)
        gatherdofnums!(u, dofnums, patchconn) # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums) # assemble symmetric matrix
    end # Loop over elements
    Kn = makematrix!(assembler)
    dofnums, B, DB, elmat, elvec, elvecfix, gradN = _buffers3(self, geom, u)
    # OPTIMIZATION: switch to a single-point quadrature rule here
    startassembly!(
        assembler,
        nodesperelem(fes) * ndofs(u),
        nodesperelem(fes) * ndofs(u),
        count(fes),
        nalldofs(u),
        nalldofs(u),
    )
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elmat, 0.0) # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            # Do the  following only if the element is well shaped and the
            # stabilization factor is positive; if the element is so distorted
            # that its Jacobian is non-positive, skip the following step.
            if self.ephis[i] > 0 && Jac != 0.0
                updatecsmat!(self.mcsys, loc, J, i, j)
                At_mul_B!(csmatTJ, csmat(self.mcsys), J) # local Jacobian matrix
                gradN!(fes, gradN, gradNparams[j], csmatTJ)
                blmat!(self.mr, B, Ns[j], gradN, loc, csmat(self.mcsys))
                add_btdb_ut_only!(elmat, B, self.ephis[i] * Jac * w[j], Dstab, DB)
            end
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(u, dofnums, fes.conn[i]) # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums) # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler) + Kn
end

function stiffness(
    self::FEMM,
    geom::NodalField{GFT},
    u::NodalField{T},
) where {FEMM<:AbstractFEMMDeforLinearESNICE,GFT<:Number,T<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return stiffness(self, assembler, geom, u)
end

"""
    inspectintegpoints(
        self::FEMM,
        geom::NodalField{GFT},
        u::NodalField{UFT},
        dT::NodalField{TFT},
        felist::Vector{IT},
        inspector::F,
        idat,
        quantity = :Cauchy;
        context...,
    ) where {FEMM<:AbstractFEMMDeforLinearESNICE,GFT<:Number,UFT<:Number,TFT<:Number,IT<:Integer,F<:Function}

Inspect integration point quantities.

# Arguments

  - `geom` - reference geometry field
  - `u` - displacement field
  - `dT` - temperature difference field
  - `felist` - indexes of the finite elements that are to be inspected:
    The fes to be included are: `fes[felist]`.
  - `context`    - structure: see the update!() method of the material.
  - `inspector` - functionwith the signature
    idat = inspector(idat, j, conn, x, out, loc);
    where
    `idat` - a structure or an array that the inspector may
    use to maintain some state,  for instance minimum or maximum of
    stress, `j` is the element number, `conn` is the element connectivity,
    `out` is the output of the update!() method,  `loc` is the location
    of the integration point in the *reference* configuration.

# Return

The updated inspector data is returned.
"""
function inspectintegpoints(
    self::FEMM,
    geom::NodalField{GFT},
    u::NodalField{UFT},
    dT::NodalField{TFT},
    felist::Vector{IT},
    inspector::F,
    idat,
    quantity = :Cauchy;
    context...,
) where {FEMM<:AbstractFEMMDeforLinearESNICE,GFT<:Number,UFT<:Number,TFT<:Number,IT<:Integer,F<:Function}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    ecoords, dofnums, loc, J, csmatTJ, Jac, D, Dstab = _buffers2(self, geom, u, npts)
    # Sort out  the output requirements
    outputcsys = self.mcsys # default: report the stresses in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = val
        end
    end
    t = 0.0
    dt = 0.0
    qpdT = 0.0 # node temperature increment
    qpstrain = fill(zero(GFT), nstressstrain(self.mr), 1) # total strain -- buffer
    qpthstrain = fill(zero(GFT), nthermstrain(self.mr)) # thermal strain -- buffer
    qpstress = fill(zero(GFT), nstressstrain(self.mr)) # stress -- buffer
    out1 = fill(zero(GFT), nstressstrain(self.mr)) # stress -- buffer
    out = fill(zero(GFT), nstressstrain(self.mr))# output -- buffer
    outtot = fill(zero(GFT), nstressstrain(self.mr))# output -- buffer
    xe = fill(0.0, nodesperelem(fes), ndofs(geom))
    eue = fill(zero(UFT), nodesperelem(fes) * ndofs(u))
    gradN = fill(zero(GFT), nodesperelem(fes), ndofs(u))
    B = fill(zero(GFT), nstressstrain(self.mr), nodesperelem(fes) * ndofs(u))
    # Loop over  all the elements and all the quadrature points within them
    for ilist = 1:length(felist) # Loop over elements
        i = felist[ilist]
        gathervalues_asmat!(geom, ecoords, fes.conn[i])# retrieve element coords
        for nix in fes.conn[i] # For all nodes connected by this element
            nodalgradN = self.nodalbasisfunctiongrad[nix].gradN
            patchconn = self.nodalbasisfunctiongrad[nix].patchconn
            Vpatch = self.nodalbasisfunctiongrad[nix].Vpatch
            ue = fill(zero(UFT), length(patchconn) * ndofs(u))
            gathervalues_asvec!(u, ue, patchconn)# retrieve element displacements
            loc = reshape(geom.values[nix, :], 1, ndofs(geom))
            updatecsmat!(self.mcsys, loc, J, nix, 1)
            nd = length(patchconn) * ndofs(u)
            Bnodal = fill(0.0, size(D, 1), nd)
            blmat!(self.mr, Bnodal, Ns[1], nodalgradN, loc, csmat(self.mcsys))
            updatecsmat!(outputcsys, loc, J, nix, 1) # Update output coordinate system
            # Quadrature point quantities
            A_mul_B!(qpstrain, Bnodal, ue) # strain in material coordinates
            qpdT = dT.values[nix] # Quadrature point temperature increment
            thermalstrain!(self.material, qpthstrain, qpdT)
            # Real Material: update the state and return the output
            out = update!(
                self.material,
                qpstress,
                out,
                vec(qpstrain),
                qpthstrain,
                t,
                dt,
                loc,
                nix,
                quantity,
            )
            copyto!(outtot, out)
            # Stabilization Material: update the state and return the output
            out = update!(
                self.stabilization_material,
                qpstress,
                out,
                vec(qpstrain),
                qpthstrain,
                t,
                dt,
                loc,
                nix,
                quantity,
            )
            outtot .+= -self.nphis[nix] .* out
            pci = findfirst(cx -> cx == nix, fes.conn[i])# at which node are we?
            locjac!(loc, J, ecoords, Ns[pci], gradNparams[pci])
            updatecsmat!(self.mcsys, loc, J, i, 1)
            At_mul_B!(csmatTJ, csmat(self.mcsys), J) # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[pci], csmatTJ)
            blmat!(self.mr, B, Ns[pci], gradN, loc, csmat(self.mcsys))
            gathervalues_asvec!(u, eue, fes.conn[i])# retrieve element displacements
            A_mul_B!(qpstrain, B, eue) # strain in material coordinates
            out = update!(
                self.stabilization_material,
                qpstress,
                out,
                vec(qpstrain),
                qpthstrain,
                t,
                dt,
                loc,
                nix,
                quantity,
            )
            outtot .+= +self.ephis[i] .* out
            out, outtot = outtot, out # swap in the total output
            if (quantity == :Cauchy)   # Transform stress tensor,  if that is "out"
                (length(out1) >= length(out)) || (out1 = zeros(length(out)))
                rotstressvec!(self.mr, out1, out, transpose(csmat(self.mcsys)))# To global coord sys
                rotstressvec!(self.mr, out, out1, csmat(outputcsys))# To output coord sys
            end
            # Call the inspector
            idat = inspector(idat, i, fes.conn[i], ecoords, out, loc)
        end # Loop over nodes
    end # Loop over elements
    return idat
end

function inspectintegpoints(
    self::FEMM,
    geom::NodalField{GFT},
    u::NodalField{UFT},
    felist::Vector{IT},
    inspector::F,
    idat,
    quantity = :Cauchy;
    context...,
) where {FEMM<:AbstractFEMMDeforLinearESNICE,GFT<:Number,UFT<:Number,IT<:Integer,F<:Function}
    dT = NodalField(fill(zero(GFT), nnodes(geom), 1)) # zero difference in temperature
    return inspectintegpoints(
        self,
        geom,
        u,
        dT,
        felist,
        inspector,
        idat,
        quantity;
        context...,
    )
end

"""
    infsup_gh(
        self::FEMM,
        assembler::A,
        geom::NodalField{GFT},
        u::NodalField{UFT},
    ) where {FEMM<:AbstractFEMMDeforLinearESNICE,A<:AbstractSysmatAssembler,GFT<:Number,UFT<:Number}

Compute the matrix to produce the norm of the divergence of the displacement.

This matrix is used in the numerical infsup test (Klaus-Jurgen Bathe, The
inf-sup condition and its evaluation for mixed finite element methods, Computers
and Structures 79 (2001) 243-252.)

!!! note 
    This computation has not been optimized in any way. It can be expected
    to be inefficient.
"""
function infsup_gh(
    self::FEMM,
    assembler::A,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {FEMM<:AbstractFEMMDeforLinearESNICE,A<:AbstractSysmatAssembler,GFT<:Number,UFT<:Number}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    elmatsizeguess = 4 * nodesperelem(fes) * ndofs(u)
    startassembly!(
        assembler,
        elmatsizeguess,
        elmatsizeguess,
        nnodes(u) + count(fes),
        u.nfreedofs,
        u.nfreedofs,
    )
    for nix = 1:length(self.nodalbasisfunctiongrad)
        gradN = self.nodalbasisfunctiongrad[nix].gradN
        patchconn = self.nodalbasisfunctiongrad[nix].patchconn
        Vpatch = self.nodalbasisfunctiongrad[nix].Vpatch
        c = reshape(geom.values[nix, :], 1, ndofs(geom))
        nd = length(patchconn) * ndofs(u)
        divm = divmat(self.mr, Ns[1], gradN, c)
        elmat = (transpose(divm) * divm) * Vpatch
        dofnums = fill(0, nd)
        gatherdofnums!(u, dofnums, patchconn) # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums) # assemble symmetric matrix
    end # Loop over nodes
    return makematrix!(assembler)
end

function infsup_gh(
    self::AbstractFEMMDeforLinearESNICE,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {GFT<:Number,UFT<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return infsup_gh(self, assembler, geom, u)
end

"""
    infsup_sh(
        self::AbstractFEMMDeforLinearESNICE,
        assembler::A,
        geom::NodalField{GFT},
        u::NodalField{UFT},
    ) where {A<:AbstractSysmatAssembler,GFT<:Number,UFT<:Number}

Compute the matrix to produce the seminorm of the displacement (square root of
the sum of the squares of the derivatives of the components of displacement).

This matrix is used in the numerical infsup test (Klaus-Jurgen Bathe, The
inf-sup condition and its evaluation for mixed finite element methods,
Computers and Structures 79 (2001) 243-252.)

!!! note
    This computation has not been optimized in any way. It can be expected to be inefficient.
"""
function infsup_sh(
    self::FEMM,
    assembler::A,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {FEMM<:AbstractFEMMDeforLinearESNICE,A<:AbstractSysmatAssembler,GFT<:Number,UFT<:Number}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    elmatsizeguess = 4 * nodesperelem(fes) * ndofs(u)
    startassembly!(
        assembler,
        elmatsizeguess,
        elmatsizeguess,
        nnodes(u) + count(fes),
        u.nfreedofs,
        u.nfreedofs,
    )
    for nix = 1:length(self.nodalbasisfunctiongrad)
        gradN = self.nodalbasisfunctiongrad[nix].gradN
        patchconn = self.nodalbasisfunctiongrad[nix].patchconn
        Vpatch = self.nodalbasisfunctiongrad[nix].Vpatch
        c = reshape(geom.values[nix, :], 1, ndofs(geom))
        nd = length(patchconn) * ndofs(u)
        vgradm = vgradmat(self.mr, Ns[1], gradN, c)
        elmat = (transpose(vgradm) * vgradm) * Vpatch
        dofnums = fill(0, nd)
        gatherdofnums!(u, dofnums, patchconn) # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums) # assemble symmetric matrix
    end # Loop over nodes
    return makematrix!(assembler)
end

function infsup_sh(
    self::FEMM,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {FEMM<:AbstractFEMMDeforLinearESNICE,GFT<:Number,UFT<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return infsup_sh(self, assembler, geom, u)
end

end
