module uxo_mode_examples
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule

using LinearAlgebra
using Arpack
using CSV

function uxo_mode_esnice_t4()
    
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0*2*pi)^2;
    
    MR = DeforModelRed3D
    output = import_NASTRAN(joinpath(@__DIR__, "UXO.nas"))
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    @show count(fens)
     
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, which=:SM, explicittransform=:none)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    true
    
end # uxo_modes


function uxo_mode_swept_h8()
    
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (1.0*2*pi)^2;
    
    MR = DeforModelRed3D
    xyzrows = CSV.File(joinpath(@__DIR__, "UXO-swept-mesh-xyz.csv"), header=0)
    xyz = fill(0.0, length(xyzrows), 3)
    for i in 1:size(xyz, 1)
        xyz[i, :] .= xyzrows[i]
    end
    connrows = CSV.File(joinpath(@__DIR__, "UXO-swept-mesh-conn.csv"), header=0)
    conn = fill(0, length(connrows), 8)
    for i in 1:size(conn, 1)
        conn[i, :] .= connrows[i]
    end
    fens = FENodeSet(xyz)
    fes = FESetH8(conn)
    
    # fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
     
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    # femm = FEMMDeforLinearESNICEH8(MR, IntegDomain(fes, NodalTensorProductRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev = eigs(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, which=:SM, explicittransform=:none)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")

    File =  "uxo_mode_swept_h8.vtk"
    vtkexportmesh(File, fens, fes)

    for mode = 1:7
    scattersysvec!(u, v[:,mode])
    File =  "alum_cyl_mode-$(mode).vtk"
    vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
end
    
    true
    
end # uxo_modes

function uxo_modes_algo()
    println("""
    % Vibration modes of unit cube  of almost incompressible material.
    %
    % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    % tetrahedral. International Journal for Numerical Methods in
    % Engineering 67: 841-867.""")
    t0 = time()
    
    
    E = 1*phun("PA");
    nu = 0.499;
    rho= 1*phun("KG/M^3");
    a=1*phun("M"); b=a; h= a;
    n1=2;# How many element edges per side?
    na= n1; nb= n1; nh =n1;
    neigvs=20                   # how many eigenvalues
    omega_shift=(0.1*2*pi)^2;
    
    fens,fes =H20block(a,b,h, na,nb,nh)
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)),
    material), "femm_mass"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)),
    material))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "omega_shift"=>omega_shift, "neigvs"=>neigvs)
    
    # Solve
    modeldata =  AlgoDeforLinearModule.modal(modeldata)
    
    fs = modeldata["omega"]/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    modeldata["postprocessing"] = FDataDict("file"=>"uxo_mode",
    "mode"=>10)
    modeldata= AlgoDeforLinearModule.exportmode(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
    
    true
    
end # uxo_modes_algo


function uxo_modes_export()
    println("""
    Vibration modes of unit cube  of almost incompressible material.
    
    This example EXPORTS the model to Abaqus.
    
    Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    tetrahedral. International Journal for Numerical Methods in
    Engineering 67: 841-867.
    """)
    t0 = time()
    
    
    E = 1*phun("PA");
    nu = 0.499;
    rho = 1*phun("KG/M^3");
    a = 1*phun("M"); b = a; h =  a;
    n1 = 5;# How many element edges per side?
    na =  n1; nb =  n1; nh  = n1;
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (0.01*2*pi)^2;
    
    MR = DeforModelRed3D
    fens,fes  = H20block(a,b,h, na,nb,nh)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material=MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)), material)
    
    K =stiffness(femm, geom, u)
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
    M =mass(femm, geom, u)
    d,v,nev,nconv = eigs(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, which=:SM, explicittransform=:none)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    mode = 7
    scattersysvec!(u, v[:,mode])
    File =  "uxo_modes.vtk"
    vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
    @async run(`"paraview.exe" $File`)
    
    AE = AbaqusExporter("uxo_modes_h20");
    # AE.ios = STDOUT;
    HEADING(AE, "Vibration modes of unit cube  of almost incompressible material.");
    COMMENT(AE, "The  first six frequencies are rigid body modes.");
    COMMENT(AE, "The  first nonzero frequency (7) should be around 0.26 Hz");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    COMMENT(AE, "The hybrid form of the serendipity hexahedron is chosen because");
    COMMENT(AE, "the material is  nearly incompressible.");
    ELEMENT(AE, "c3d20rh", "AllElements", 1, connasarray(fes))
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements");
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    DENSITY(AE, rho)
    STEP_FREQUENCY(AE, neigvs)
    END_STEP(AE)
    close(AE)
    
    true
    
end # uxo_modes_export


function uxo_modes_msh8_algo()
    println("""
    % Vibration modes of unit cube  of almost incompressible material.
    % Mean-strain hexahedron.
    % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    % tetrahedral. International Journal for Numerical Methods in
    % Engineering 67: 841-867.""")
    t0 = time()
    
    
    E = 1*phun("PA");
    nu = 0.499;
    rho= 1*phun("KG/M^3");
    a=1*phun("M"); b=a; h= a;
    n1=8 # How many element edges per side?
    na= n1; nb= n1; nh =n1;
    neigvs=20                   # how many eigenvalues
    omega_shift=(0.1*2*pi)^2;
    
    fens,fes = H8block(a,b,h, na,nb,nh)
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),
    material), "femm_mass"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,3)),
    material))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "omega_shift"=>omega_shift, "neigvs"=>neigvs)
    
    # Solve
    modeldata =  AlgoDeforLinearModule.modal(modeldata)
    
    fs = modeldata["omega"]/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    modeldata["postprocessing"] = FDataDict("file"=>"uxo_mode",
    "mode"=>10)
    modeldata= AlgoDeforLinearModule.exportmode(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
    
    true
    
end # uxo_modes_msh8_algo

function allrun()
    println("#####################################################") 
    println("# uxo_mode_nice_t4 omitted due to significant size")
    # uxo_mode_esnice_t4()
    println("#####################################################") 
    println("# uxo_mode_esnice_t4 ")
    uxo_mode_swept_h8()
    uxo_modes_algo()
    uxo_modes_export()
    uxo_modes_msh8_algo()
    return true
end # function allrun


@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module uxo_mode_examples
nothing
