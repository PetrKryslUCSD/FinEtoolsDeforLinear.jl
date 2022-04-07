module stubby_corbel_examples
using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule: linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra
using Krylov, LinearOperators, IncompleteLU
using SparseArrays
using SuiteSparse
using Printf
using SymRCM
using UnicodePlots
using Infiltrator
using Random
using ILUZero

# Isotropic material
E=1000.0;
nu=0.4999; #Taylor data
W=25.0;
H=50.0;
L= 50.0;
htol = minimum([L,H,W])/1000;
uzex =-0.16;
magn = 0.2*(-12.6)/4;
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

function stubby_corbel_H8_big(n = 10, solver = :suitesparse)
    elementtag = "H8"
    println("""
    Stubby corbel example. Element: $(elementtag)
    """)

    fens,fes = H8block(W, L, H, n, 2*n, 2*n)
    println("Number of elements: $(count(fes))")
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
    K = stiffness(femm, geom, u)
    println("Stiffness: number of non zeros = $(nnz(K)) [ND]")
    println("Sparsity = $(nnz(K)/size(K, 1)/size(K, 2)) [ND]")
    display(spy(K, canvas = DotCanvas))

    if solver == :suitesparse
    # @show methods(SuiteSparse.CHOLMOD.ldlt, (typeof(K), ))
        # @time K = SuiteSparse.CHOLMOD.ldlt(K)
        @time K = SuiteSparse.CHOLMOD.cholesky(K)
    # @time K = SparseArrays.ldlt(K)
    # @time K = cholesky(K)
    # @infiltrate 
        @time U = K\(F2)
    else
        n = size(K, 1)
        mKd = mean(diag(K))
        # @time factor = ilu(K, τ = mKd / 1000.0) # This may work for compressible materials
        @time factor = ilu(K, τ = mKd / 1000000.0) # This may work for incompressible materials
        # factor = ilu0(K)
        @show nnz(factor) / nnz(K)
        opM = LinearOperator(Float64, n, n, false, false, (y, v) -> ldiv!(y, factor, v))
        @time (U, stats) = Krylov.cg(K, F2;
                        M=opM, 
                        itmax=Int(round(n/2)), 
                        verbose=1)
    end
    scattersysvec!(u,U[:])


    # @show length(U)
    Tipl = selectnode(fens, box=[0 W L L 0 H], inflate=htol)
    utip = mean(u.values[Tipl, 3], dims=1)
    println("Deflection: $(utip), compared to $(uzex)")

    File =  "stubby_corbel_H8_big.vtk"
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
end # stubby_corbel_H8_big

function allrun()
    println("#####################################################")
    println("# stubby_corbel_H8_by_hand ")
    stubby_corbel_H8_by_hand()
    return true
end # function allrun


@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module stubby_corbel_examples
nothing
    
