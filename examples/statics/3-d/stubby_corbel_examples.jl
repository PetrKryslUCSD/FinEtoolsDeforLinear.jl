module stubby_corbel_examples
using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule: linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky
using SparseArrays
using Printf
using SymRCM

# Isotropic material
E=1000.0;
nu=0.4999; #Taylor data
W=25.0;
H=50.0;
L= 50.0;
htol = minimum([L,H,W])/1000;
uzex =-12.6;
magn = 0.2*uzex/4;
Force =magn*W*H*2;
CTE = 0.0
n = 5 #

function getfrcL!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    copyto!(forceout, [0.0; 0.0; magn])
end

function stubby_corbel_H8_by_hand()
    elementtag = "H8"
    println("""
    Stubby corbel example. Element: $(elementtag)
    """)

    fens,fes = H8block(W, L, H, n, 2*n, 2*n)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing=true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing=true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing=true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]

    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(csmatout, csmat)
    end

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u,lx0,true,1,0.0)
    setebc!(u,lx0,true,2,0.0)
    setebc!(u,lx0,true,3,0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u,lx1,true,1,0.0)
    applyebc!(u)
    numberdofs!(u)
    println("u.nfreedofs = $(u.nfreedofs)")

    fi = ForceIntensity(FFlt, 3, getfrcL!);
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2);
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    @time K=cholesky(K)
    @time U = K\(F2)
    scattersysvec!(u,U[:])
    @show length(U)
    Tipl = selectnode(fens, box=[0 W L L 0 H], inflate=htol)
    utip = mean(u.values[Tipl, 3], dims=1)
    println("Deflection: $(utip), compared to $(uzex)")

    File =  "stubby_corbel_H8_by_hand.vtk"
    vtkexportmesh(File, fens, fes;  vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # stubby_corbel_H8_by_hand

function stubby_corbel_H8_big(n = 30)
    elementtag = "H8"
    println("""
    Stubby corbel example. Element: $(elementtag)
    """)

    fremem = Base.Sys.free_memory()
    @printf "Initial free memory = %.1f [MB]\n" round(fremem/1024^2, digits = 1)

    fens,fes = H8block(W, L, H, n, 2*n, 2*n)
    println("Number of elements: $(count(fes))")
    remfremem = fremem - Base.Sys.free_memory()
    @printf "After mesh, free memory decrease = %.1f [MB]\n" round(remfremem/1024^2, digits = 1)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing=true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing=true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing=true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]

    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(csmatout, csmat)
    end


    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u,lx0,true,1,0.0)
    setebc!(u,lx0,true,2,0.0)
    setebc!(u,lx0,true,3,0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u,lx1,true,1,0.0)
    applyebc!(u)
    numberdofs!(u, perm)
    println("u.nfreedofs = $(u.nfreedofs)")

    fi = ForceIntensity(FFlt, 3, getfrcL!);
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2);
    associategeometry!(femm, geom)
    remfremem = fremem - Base.Sys.free_memory()
    @printf "Before stiffness, free memory decrease = %.1f [MB]\n" round(remfremem/1024^2, digits = 1)
    K = stiffness(femm, geom, u)
    remfremem = fremem - Base.Sys.free_memory()
    @printf "After stiffness, free memory decrease = %.1f [MB]\n" round(remfremem/1024^2, digits = 1)
    println("Stiffness: number of non zeros = $(nnz(K)) [ND]")

    @time K = cholesky(Symmetric(K))
    remfremem = fremem - Base.Sys.free_memory()
    @printf "After cholesky, free memory decrease = %.1f [MB]\n" round(remfremem/1024^2, digits = 1)
    @time U = K\(F2)
    remfremem = fremem - Base.Sys.free_memory()
    @printf "After solution, free memory decrease = %.1f [MB]\n" round(remfremem/1024^2, digits = 1)
    scattersysvec!(u,U[:])
    # @show length(U)
    Tipl = selectnode(fens, box=[0 W L L 0 H], inflate=htol)
    utip = mean(u.values[Tipl, 3], dims=1)
    # println("Deflection: $(utip), compared to $(uzex)")

    # File =  "stubby_corbel_H8_by_hand.vtk"
    # vtkexportmesh(File, fens, fes;  vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # stubby_corbel_H8_big

function allrun()
    println("#####################################################")
    println("# stubby_corbel_H8_by_hand ")
    stubby_corbel_H8_by_hand()
    return true
end # function allrun

end # module stubby_corbel_examples

    stubby_corbel_examples.stubby_corbel_H8_big()
