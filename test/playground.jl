
module mholestr1
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 100.0;
    nu = 1.0/3;
    cte = 0.15
    xradius = 1.0
    yradius = 1.0
    L = 3.0
    H = 3.0
    # nL = 50
    # nH = 50
    # nW = 70
    nL = 100
    nH = 100
    nW = 120
    tolerance = min(xradius, yradius, L, H)/min(nL, nH, nW)/1000.;#Geometrical tolerance

    fens, fes = Q4elliphole(xradius, yradius, L, H,    nL, nH, nW)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    l1 = selectnode(fens; box=[0, 0, -Inf,  Inf],  inflate = tolerance)
    setebc!(u, l1, 1, 0.0)
    l2 = selectnode(fens; box=[-Inf,  Inf, 0, 0],  inflate = tolerance)
    setebc!(u, l2, 2, 0.0)
    applyebc!(u)
    numberdofs!(u)

    boundaryfes =  meshboundary(fes);
    lce = selectelem(fens, boundaryfes, facing=true, direction = x -> -x,  inflate=  tolerance);
    lc = connectednodes(subset(boundaryfes, lce))
    
    le = selectelem(fens, boundaryfes,  box= [-Inf,  Inf, H, H],  inflate=  tolerance);
    el1femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([0.0, 1.0]);
    Fm = distribloads(el1femm,  geom,  u,  fi,  2);
    le = selectelem(fens, boundaryfes,  box= [L, L, -Inf,  Inf],  inflate=  tolerance);
    el2femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([1.0, 0.0]);
    Fm += distribloads(el2femm,  geom,  u,  fi,  2);

    MR = DeforModelRed2DStress
    material = MatDeforElastIso(MR,  0.0, E, nu, cte)

    femm = FEMMDeforLinear(MR, IntegDomain(fes,  GaussRule(2, 2)),  material)
    
    K = stiffness(femm,  geom,  u)
    U=  K\(Fm)
    scattersysvec!(u, U[:])

    
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # @test norm(maximum(fld.values) - 8.190847372888073e7)/8.190847372888073e7 <= 1.0e-4
    # File =  "mholestr1-s1.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigma_1", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    fld= fieldfromintegpoints(femm, geom, u, :maxshear, 1)
    @test norm(maximum(fld.values) - 2.210557410276065)/2.210557410276065 <= 1.0e-4
    # File =  "mholestr1-maxshear.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("maxshear", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    

end
end
using .mholestr1
mholestr1.test()


module mholestr2
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 100.0;
    nu = 1.0/3;
    cte = 0.15
    xradius = 1.0
    yradius = 1.0
    L = 3.0
    H = 3.0
    # nL = 50
    # nH = 50
    # nW = 70
    nL = 100
    nH = 100
    nW = 120
    tolerance = min(xradius, yradius, L, H)/min(nL, nH, nW)/1000.;#Geometrical tolerance

    fens, fes = Q4elliphole(xradius, yradius, L, H,    nL, nH, nW)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    l1 = selectnode(fens; box=[0, 0, -Inf,  Inf],  inflate = tolerance)
    setebc!(u, l1, 1, 0.0)
    l2 = selectnode(fens; box=[-Inf,  Inf, 0, 0],  inflate = tolerance)
    setebc!(u, l2, 2, 0.0)
    applyebc!(u)
    numberdofs!(u)

    boundaryfes =  meshboundary(fes);
    lce = selectelem(fens, boundaryfes, facing=true, direction = x -> -x,  inflate=  tolerance);
    lc = connectednodes(subset(boundaryfes, lce))
    
    le = selectelem(fens, boundaryfes,  box= [-Inf,  Inf, H, H],  inflate=  tolerance);
    el1femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([0.0, 1.0]);
    Fm = distribloads(el1femm,  geom,  u,  fi,  2);
    le = selectelem(fens, boundaryfes,  box= [L, L, -Inf,  Inf],  inflate=  tolerance);
    el2femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([1.0, 0.0]);
    Fm += distribloads(el2femm,  geom,  u,  fi,  2);

    MR = DeforModelRed2DStrain
    material = MatDeforElastIso(MR,  0.0, E, nu, cte)

    femm = FEMMDeforLinear(MR, IntegDomain(fes,  GaussRule(2, 2)),  material)
    
    K = stiffness(femm,  geom,  u)
    U=  K\(Fm)
    scattersysvec!(u, U[:])

    
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # @test norm(maximum(fld.values) - 8.190847372888073e7)/8.190847372888073e7 <= 1.0e-4
    # File =  "mholestr2-s1.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigma_1", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    fld= fieldfromintegpoints(femm, geom, u, :maxshear, 1)
    @test norm(maximum(fld.values) - 2.2102738887214257)/2.2102738887214257 <= 1.0e-4
    # File =  "mholestr2-maxshear.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("maxshear", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    

end
end
using .mholestr2
mholestr2.test()

module mholestr3
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 100.0;
    nu = 1.0/3;
    cte = 0.15
    xradius = 1.0
    yradius = 1.0
    L = 3.0
    H = 3.0
    # nL = 50
    # nH = 50
    # nW = 70
    nL = 100
    nH = 100
    nW = 120
    tolerance = min(xradius, yradius, L, H)/min(nL, nH, nW)/1000.;#Geometrical tolerance

    fens, fes = Q4elliphole(xradius, yradius, L, H,    nL, nH, nW)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    l1 = selectnode(fens; box=[0, 0, -Inf,  Inf],  inflate = tolerance)
    setebc!(u, l1, 1, 0.0)
    l2 = selectnode(fens; box=[-Inf,  Inf, 0, 0],  inflate = tolerance)
    setebc!(u, l2, 2, 0.0)
    applyebc!(u)
    numberdofs!(u)

    boundaryfes =  meshboundary(fes);
    lce = selectelem(fens, boundaryfes, facing=true, direction = x -> -x,  inflate=  tolerance);
    lc = connectednodes(subset(boundaryfes, lce))
    
    le = selectelem(fens, boundaryfes,  box= [-Inf,  Inf, H, H],  inflate=  tolerance);
    el1femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([0.0, -1.0]);
    Fm = distribloads(el1femm,  geom,  u,  fi,  2);
    le = selectelem(fens, boundaryfes,  box= [L, L, -Inf,  Inf],  inflate=  tolerance);
    el2femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([1.0, 0.0]);
    Fm += distribloads(el2femm,  geom,  u,  fi,  2);

    MR = DeforModelRed2DStress
    material = MatDeforElastIso(MR,  0.0, E, nu, cte)

    femm = FEMMDeforLinear(MR, IntegDomain(fes,  GaussRule(2, 2)),  material)
    
    K = stiffness(femm,  geom,  u)
    U=  K\(Fm)
    scattersysvec!(u, U[:])

    
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # @test norm(maximum(fld.values) - 8.190847372888073e7)/8.190847372888073e7 <= 1.0e-4
    File =  "mholestr3-s1.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigma_1", fld.values)], vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)
    fld= fieldfromintegpoints(femm, geom, u, :maxshear, 1)
    @test norm(maximum(fld.values) - 5.921999943843146)/5.921999943843146 <= 1.0e-4
    # File =  "mholestr3-maxshear.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("maxshear", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    

end
end
using .mholestr3
mholestr3.test()

