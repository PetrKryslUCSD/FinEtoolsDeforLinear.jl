"""
    AbstractFEMMDeforLinearBaseModule

Base module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models.
"""
module FEMMDeforLinearBaseModule

__precompile__(true)

using FinEtools.FENodeSetModule: FENodeSet
using FinEtools.DataCacheModule: DataCache
using FinEtools.FESetModule: gradN!, nodesperelem, manifdim
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using FinEtools.FieldModule:
    ndofs,
    gatherdofnums!,
    gatherfixedvalues_asvec!,
    gathervalues_asvec!,
    gathervalues_asmat!,
    nalldofs
using FinEtools.NodalFieldModule: NodalField, nnodes
using FinEtools.AssemblyModule:
    AbstractSysvecAssembler,
    AbstractSysmatAssembler,
    SysmatAssemblerSparseSymm,
    startassembly!,
    assemble!,
    makematrix!,
    makevector!,
    SysvecAssembler
using FinEtools.FEMMBaseModule: AbstractFEMM
import FinEtools.FEMMBaseModule: inspectintegpoints, bilform_dot, bilform_lin_elastic
using FinEtools.CSysModule: CSys, updatecsmat!, csmat
using FinEtools.DeforModelRedModule: nstressstrain, nthermstrain, blmat!, divmat, vgradmat
using FinEtools.MatrixUtilityModule:
    add_btdb_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, add_btsigma!
using FinEtoolsDeforLinear.MatDeforModule: rotstressvec!
using FinEtools.MatModule: massdensity
using FinEtoolsDeforLinear.MatDeforLinearElasticModule:
    tangentmoduli!, update!, thermalstrain!
using FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
using LinearAlgebra: Transpose, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
using LinearAlgebra: norm, dot, I
using LinearAlgebra

"""
    AbstractFEMMDeforLinear <: AbstractFEMMBase

Abstract type of FEMM for linear deformation.
"""
abstract type AbstractFEMMDeforLinear <: AbstractFEMM end

function _buffers(
    self::AbstractFEMMDeforLinear,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {GFT,UFT}
    fes = self.integdomain.fes
    ndn = ndofs(u) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    nstrs = nstressstrain(self.mr)  # number of stresses
    elmatdim = ndn * nne             # dimension of the element matrix
    # Prepare _buffers
    ecoords = fill(zero(GFT), nne, ndofs(geom)) # array of Element coordinates
    elmat = fill(zero(GFT), elmatdim, elmatdim)      # element matrix -- buffer
    dofnums = zeros(eltype(u.dofnums), elmatdim) # degree of freedom array -- buffer
    loc = fill(zero(GFT), 1, sdim) # quadrature point location -- buffer
    J = fill(zero(GFT), sdim, mdim) # Jacobian matrix -- buffer
    csmatTJ = fill(zero(GFT), mdim, mdim) # intermediate result -- buffer
    gradN = fill(zero(GFT), nne, mdim) # intermediate result -- buffer
    D = fill(zero(GFT), nstrs, nstrs) # material stiffness matrix -- buffer
    B = fill(zero(GFT), nstrs, elmatdim) # strain-displacement matrix -- buffer
    DB = fill(zero(GFT), nstrs, elmatdim) # strain-displacement matrix -- buffer
    elvecfix = fill(zero(UFT), elmatdim) # vector of prescribed displ. -- buffer
    elvec = fill(zero(UFT), elmatdim) # element vector -- buffer
    return ecoords, dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix
end

"""
    mass(self::AbstractFEMMDeforLinear,  assembler::A,  geom::NodalField{GFT}, u::NodalField{UFT}) where {A<:AbstractSysmatAssembler, GFT<:Number, UFT<:Number}

Compute the consistent mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function mass(
    self::AbstractFEMMDeforLinear,
    assembler::A,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {A<:AbstractSysmatAssembler,GFT<:Number,UFT<:Number}
    cf = DataCache(massdensity(self.material) * LinearAlgebra.I(ndofs(u)))
    return bilform_dot(self, assembler, geom, u, cf)
end

function mass(
    self::AbstractFEMMDeforLinear,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {GFT<:Number,UFT<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return mass(self, assembler, geom, u)
end

"""
    stiffness(self::FEMM, assembler::A, geom::NodalField{GFT}, u::NodalField{UFT}) where {FEMM<:AbstractFEMMDeforLinear, A<:AbstractSysmatAssembler, GFT<:Number, UFT<:Number}

Compute and assemble  stiffness matrix.

!!! note

    The material stiffness matrix is assumed to be the same at all the points of
    the domain (homogeneous material).
"""
function stiffness(
    self::FEMM,
    assembler::A,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {FEMM<:AbstractFEMMDeforLinear,A<:AbstractSysmatAssembler,GFT<:Number,UFT<:Number}
    sdim = ndofs(geom)
    loc = fill(zero(GFT), 1, sdim)
    nstrs = nstressstrain(self.mr)
    D = fill(zero(GFT), nstrs, nstrs)
    tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)
    return bilform_lin_elastic(self, assembler, geom, u, self.mr, DataCache(D))
end

function stiffness(
    self::FEMM,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {FEMM<:AbstractFEMMDeforLinear,GFT<:Number,UFT<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return stiffness(self, assembler, geom, u)
end

"""
    thermalstrainloads(self::AbstractFEMMDeforLinear, assembler::A, geom::NodalField{GFT}, u::NodalField{UFT}, dT::NodalField{TFT}) where {A<:AbstractSysvecAssembler, GFT<:Number, UFT<:Number, TFT<:Number}

Compute the thermal-strain load vector.
"""
function thermalstrainloads(
    self::AbstractFEMMDeforLinear,
    assembler::A,
    geom::NodalField{GFT},
    u::NodalField{UFT},
    dT::NodalField{TFT},
) where {A<:AbstractSysvecAssembler,GFT<:Number,UFT<:Number,TFT<:Number}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    ecoords, dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix =
        _buffers(self, geom, u)
    t = 0.0
    dt = 0.0
    DeltaT = fill(zero(GFT), nodesperelem(fes))
    strain = fill(zero(GFT), nstressstrain(self.mr)) # total strain -- buffer
    thstrain = fill(zero(GFT), nthermstrain(self.mr)) # thermal strain -- buffer
    thstress = fill(zero(GFT), nstressstrain(self.mr)) # thermal stress -- buffer
    startassembly!(assembler, nalldofs(u))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asvec!(dT, DeltaT, fes.conn[i])# retrieve element temperatures
        if norm(DeltaT, Inf) != 0     # Is the thermal increment nonzero?
            gathervalues_asmat!(geom, ecoords, fes.conn[i])
            fill!(elvec, 0.0) # Initialize element matrix
            for j = 1:npts # Loop over quadrature points
                locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
                updatecsmat!(self.mcsys, loc, J, i, j)
                At_mul_B!(csmatTJ, csmat(self.mcsys), J) # local Jacobian matrix
                gradN!(fes, gradN, gradNparams[j], csmatTJ)#Do: gradN = gradNparams[j]/csmatTJ;
                blmat!(self.mr, B, Ns[j], gradN, loc, csmat(self.mcsys))# strains in mater cs, displ in global cs
                thermalstrain!(self.material, thstrain, dot(vec(Ns[j]), DeltaT))
                thstress = update!(
                    self.material,
                    thstress,
                    thstress,
                    strain,
                    thstrain,
                    t,
                    dt,
                    loc,
                    fes.label[i],
                    :nothing,
                )
                add_btsigma!(elvec, B, (-1) * (Jac * w[j]), thstress)
            end
            gatherdofnums!(u, dofnums, fes.conn[i]) # retrieve degrees of freedom
            assemble!(assembler, elvec, dofnums) # assemble element load vector
        end
    end # Loop over elements
    return makevector!(assembler)
end

function thermalstrainloads(
    self::AbstractFEMMDeforLinear,
    geom::NodalField{GFT},
    u::NodalField{UFT},
    dT::NodalField{TFT},
) where {GFT<:Number,UFT<:Number,TFT<:Number}
    assembler = SysvecAssembler()
    return thermalstrainloads(self, assembler, geom, u, dT)
end

"""
    inspectintegpoints(self::FEMM, geom::NodalField{GFT},  u::NodalField{UFT}, dT::NodalField{TFT}, felist::AbstractVector{IT}, inspector::F, idat, quantity=:Cauchy; context...) where {FEMM<:AbstractFEMMDeforLinear, GFT<:Number, UFT<:Number, TFT<:Number, IT<:Integer, F<:Function}

Inspect integration point quantities.

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

### Return

The updated inspector data is returned.
"""
function inspectintegpoints(
    self::FEMM,
    geom::NodalField{GFT},
    u::NodalField{UFT},
    dT::NodalField{TFT},
    felist::AbstractVector{IT},
    inspector::F,
    idat,
    quantity = :Cauchy;
    context...,
) where {
    FEMM<:AbstractFEMMDeforLinear,
    GFT<:Number,
    UFT<:Number,
    TFT<:Number,
    IT<:Integer,
    F<:Function,
}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    ecoords, dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix =
        _buffers(self, geom, u)
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
    dTe = fill(zero(TFT), nodesperelem(fes)) # nodal temperatures -- buffer
    ue = fill(zero(GFT), size(elmat, 1)) # array of node displacements -- buffer
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    xe = fill(zero(GFT), nne, sdim) # array of node coordinates -- buffer
    qpdT = 0.0 # node temperature increment
    qpstrain = fill(zero(GFT), nstressstrain(self.mr), 1) # total strain -- buffer
    qpthstrain = fill(zero(GFT), nthermstrain(self.mr)) # thermal strain -- buffer
    qpstress = fill(zero(GFT), nstressstrain(self.mr)) # stress -- buffer
    out1 = fill(zero(GFT), nstressstrain(self.mr)) # stress -- buffer
    out = fill(zero(GFT), nstressstrain(self.mr))# output -- buffer
    # Loop over  all the elements and all the quadrature points within them
    for ilist in eachindex(felist) # Loop over elements
        i = felist[ilist]
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        gathervalues_asvec!(u, ue, fes.conn[i])# retrieve element displacements
        gathervalues_asvec!(dT, dTe, fes.conn[i])# retrieve element temp. increments
        for j  in 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            updatecsmat!(self.mcsys, loc, J, i, j)
            At_mul_B!(csmatTJ, csmat(self.mcsys), J) # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], csmatTJ)
            blmat!(self.mr, B, Ns[j], gradN, loc, csmat(self.mcsys))
            updatecsmat!(outputcsys, loc, J, i, j)
            # Quadrature point quantities
            A_mul_B!(qpstrain, B, ue) # strain in material coordinates
            qpdT = dot(vec(dTe), vec(Ns[j]))# Quadrature point temperature increment
            thermalstrain!(self.material, qpthstrain, qpdT)
            # Material updates the state and returns the output
            out = update!(
                self.material,
                qpstress,
                out,
                vec(qpstrain),
                qpthstrain,
                t,
                dt,
                loc,
                fes.label[i],
                quantity,
            )
            if (quantity == :Cauchy)   # Transform stress tensor,  if that is "out"
                (length(out1) >= length(out)) || (out1 = zeros(length(out)))
                rotstressvec!(self.mr, out1, out, transpose(csmat(self.mcsys)))# To global coord sys
                rotstressvec!(self.mr, out, out1, csmat(outputcsys))# To output coord sys
            end
            # Call the inspector
            idat = inspector(idat, i, fes.conn[i], ecoords, out, loc)
        end # Loop over quadrature points
    end # Loop over elements
    return idat # return the updated inspector data
end

function inspectintegpoints(
    self::FEMM,
    geom::NodalField{GFT},
    u::NodalField{UFT},
    felist::AbstractVector{IT},
    inspector::F,
    idat,
    quantity = :Cauchy;
    context...,
) where {FEMM<:AbstractFEMMDeforLinear,GFT<:Number,UFT<:Number,IT,F<:Function}
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

function _buffers2(
    self::AbstractFEMMDeforLinear,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {GFT,UFT}
    fes = self.integdomain.fes
    ndn = ndofs(u) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    nstrs = nstressstrain(self.mr)  # number of stresses
    elmatdim = ndn * nne             # dimension of the element matrix
    # Prepare buffers
    ecoords = fill(zero(GFT), nne, ndofs(geom)) # array of Element coordinates
    elmat = fill(zero(GFT), elmatdim, elmatdim)      # element matrix -- buffer
    dofnums = zeros(eltype(u.dofnums), elmatdim) # degree of freedom array -- buffer
    loc = fill(zero(GFT), 1, sdim) # quadrature point location -- buffer
    J = fill(zero(GFT), sdim, mdim) # Jacobian matrix -- buffer
    csmatTJ = fill(zero(GFT), mdim, mdim) # intermediate result -- buffer
    gradN = fill(zero(GFT), nne, mdim) # intermediate result -- buffer
    return ecoords, dofnums, loc, J, csmatTJ, gradN, elmat
end

"""
    infsup_gh(self::AbstractFEMMDeforLinear, assembler::A, geom::NodalField{GFT}, u::NodalField{UFT}) where {A<:AbstractSysmatAssembler, GFT, UFT}

Compute the matrix to produce the norm of the divergence of the displacement.

This matrix is used in the numerical infsup test (Klaus-Jurgen Bathe, The
inf-sup condition and its evaluation for mixed finite element methods,
Computers and Structures 79 (2001) 243-252.)

!!! note


This computation has not been optimized in any way. It can be expected to be
inefficient.
"""
function infsup_gh(
    self::AbstractFEMMDeforLinear,
    assembler::A,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {A<:AbstractSysmatAssembler,GFT,UFT}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    ecoords, dofnums, loc, J, csmatTJ, gradN, elmat = _buffers2(self, geom, u)
    startassembly!(assembler, size(elmat)..., count(fes), nalldofs(u), nalldofs(u))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elmat, 0.0) # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            updatecsmat!(self.mcsys, loc, J, i, j)
            At_mul_B!(csmatTJ, csmat(self.mcsys), J) # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], csmatTJ)
            divm = divmat(self.mr, Ns[j], gradN, loc)
            elmat += (transpose(divm) * divm) * (Jac * w[j])
        end # Loop over quadrature points
        gatherdofnums!(u, dofnums, fes.conn[i]) # retrieve degrees of freedom
        assemble!(assembler, (elmat + elmat') / 2, dofnums, dofnums) # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function infsup_gh(
    self::AbstractFEMMDeforLinear,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {GFT,UFT}
    assembler = SysmatAssemblerSparseSymm()
    return infsup_gh(self, assembler, geom, u)
end

function _buffers3(
    self::AbstractFEMMDeforLinear,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {GFT,UFT}
    fes = self.integdomain.fes
    ndn = ndofs(u) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    nstrs = nstressstrain(self.mr)  # number of stresses
    elmatdim = ndn * nne             # dimension of the element matrix
    # Prepare buffers
    ecoords = fill(zero(GFT), nne, ndofs(geom)) # array of Element coordinates
    elmat = fill(zero(GFT), elmatdim, elmatdim)      # element matrix -- buffer
    dofnums = zeros(eltype(u.dofnums), elmatdim) # degree of freedom array -- buffer
    loc = fill(zero(GFT), 1, sdim) # quadrature point location -- buffer
    J = fill(zero(GFT), sdim, mdim) # Jacobian matrix -- buffer
    csmatTJ = fill(zero(GFT), mdim, mdim) # intermediate result -- buffer
    gradN = fill(zero(GFT), nne, mdim) # intermediate result -- buffer
    return ecoords, dofnums, loc, J, csmatTJ, gradN, elmat
end

"""
    infsup_sh(self::AbstractFEMMDeforLinear, assembler::A, geom::NodalField{GFT}, u::NodalField{UFT}) where {A<:AbstractSysmatAssembler, GFT<:Number, UFT<:Number}

Compute the matrix to produce the seminorm of the displacement (square root of
the sum of the squares of the derivatives of the components of displacement).

This matrix is used in the numerical infsup test (Klaus-Jurgen Bathe, The
inf-sup condition and its evaluation for mixed finite element methods,
Computers and Structures 79 (2001) 243-252.)

!!! note


This computation has not been optimized in any way. It can be expected to be
inefficient.
"""
function infsup_sh(
    self::AbstractFEMMDeforLinear,
    assembler::A,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {A<:AbstractSysmatAssembler,GFT<:Number,UFT<:Number}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    ecoords, dofnums, loc, J, csmatTJ, gradN, elmat = _buffers3(self, geom, u)
    startassembly!(assembler, size(elmat)..., count(fes), nalldofs(u), nalldofs(u))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elmat, 0.0) # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            updatecsmat!(self.mcsys, loc, J, i, j)
            At_mul_B!(csmatTJ, csmat(self.mcsys), J) # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], csmatTJ)
            vgradm = vgradmat(self.mr, Ns[j], gradN, loc)
            elmat += (transpose(vgradm) * vgradm) * (Jac * w[j])
        end # Loop over quadrature points
        gatherdofnums!(u, dofnums, fes.conn[i]) # retrieve degrees of freedom
        assemble!(assembler, (elmat + elmat') / 2, dofnums, dofnums) # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function infsup_sh(
    self::AbstractFEMMDeforLinear,
    geom::NodalField{GFT},
    u::NodalField{UFT},
) where {GFT<:Number,UFT<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return infsup_sh(self, assembler, geom, u)
end

end
