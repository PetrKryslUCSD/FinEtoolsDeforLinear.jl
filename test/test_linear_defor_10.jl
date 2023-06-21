module mbaruniax1
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.FENodeSetModule
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 30e6*phun("psi") # Young's modulus
    alpha = 5.5e-6*phun("in")/phun("in")/phun("F")
    for nu in [0.0 0.3 0.45]

        MR = DeforModelRed1DStress
        material = MatDeforElastIso(MR,  0.0, E, nu, alpha)
        D = fill(0.0, 1, 1)
        t::FFlt, dt::FFlt, loc::FFltMat, label::FInt = 0.0, 0.0, [0.0 0.0 0.0], 0
        tangentmoduli!(material,  D,  t, dt, loc, label)
        @test D[1, 1] ≈ E
    end
    
    true
end
end
using .mbaruniax1
mbaruniax1.test()

module mbaruniax2
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.FENodeSetModule
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 30e6*phun("psi") # Young's modulus
    alpha = 5.5e-6*phun("in")/phun("in")/phun("F")
    for nu in [0.0 0.3 0.45]

        MR = DeforModelRed1DStrain
        material = MatDeforElastIso(MR,  0.0, E, nu, alpha)
        D = fill(0.0, 1, 1)
        t::FFlt, dt::FFlt, loc::FFltMat, label::FInt = 0.0, 0.0, [0.0 0.0 0.0], 0
        tangentmoduli!(material,  D,  t, dt, loc, label)
        lambda = E * nu / (1 + nu) / (1 - 2*(nu));
        mu = E / 2. / (1+nu);
        @test D[1, 1] ≈ lambda + 2*mu
    end
    
    true
end
end
using .mbaruniax2
mbaruniax2.test()

module mTEST13H_in_fluid
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsDeforLinear
using LinearAlgebra
using Arpack
using Test

function TEST13H_hva()
    # Harmonic forced vibration problem is solved for a homogeneous square plate,
    # simply-supported on the circumference.
    # This is the TEST 13H from the Abaqus v 6.12 Benchmarks manual.
    # The test is recommended by the National Agency for Finite Element Methods and Standards (U.K.):
    # Test 13 from NAFEMS “Selected Benchmarks for Forced Vibration,” R0016, March 1993.
    #
    #
    # The plate is discretized with hexahedral solid elements. The simple support
    # condition is approximated by distributed rollers on the boundary.
    # Because only the out of plane displacements are prevented, the structure
    # has three rigid body modes in the plane of the plate.
    #
    #
    # The nonzero benchmark frequencies are (in hertz): 2.377, 5.961, 5.961,
    # 9.483, 12.133, 12.133, 15.468, 15.468 [Hz].
    
    # println(""" Homogeneous square plate, simply-supported on the circumference,
    # from the test 13 from NAFEMS “Selected Benchmarks for Forced Vibration,”
    # R0016, March 1993. The nonzero benchmark frequencies are (in hertz): 2.377,
    # 5.961, 5.961, 9.483, 12.133, 12.133, 15.468, 15.468 [Hz].

    # This problem is extended by including fluid-induced damping by the
    # surrounding air using a matrix expressing the ABC with dampers along the
    # boundary.
    # """)
    
    # t0 = time()
    
    E = 200*phun("GPa");# Young's modulus
    nu = 0.3;# Poisson ratio
    rho = 8000*phun("KG*M^-3");# mass density
    qmagn = 100.0*phun("Pa")
    L = 10.0*phun("M"); # side of the square plate
    t = 0.05*phun("M"); # thickness of the square plate
    nL = 16; nt = 4;
    tolerance = t/nt/100;
    # neigvs = 11;
    # OmegaShift = (2*pi*0.5) ^ 2; # to resolve rigid body modes
    frequencies = vcat(linearspace(0.0,2.377,20), linearspace(2.377,15.0,70))
    rho_fluid = 1.3*phun("kg*m^3")
    c_fluid = 341*phun("m/s")
    
    # Compute the parameters of Rayleigh damping. For the two selected
    # frequencies we have the relationship between the damping ratio and
    # the Rayleigh parameters
    # $\xi_m=a_0/\omega_m+a_1\omega_m$
    # where $m=1,2$.  Solving for the Rayleigh parameters $a_0,a_1$ yields:
    zeta1= 0.02; zeta2  =0.02;
    o1 =2*pi*2.377;  o2 =2*pi*15.468;
    Rayleigh_mass = 2*(o1*o2)/(o2^2-o1^2)*(o2*zeta1-o1*zeta2);# a0
    Rayleigh_stiffness = 2*(o1*o2)/(o2^2-o1^2)*(-1/o2*zeta1+1/o1*zeta2);# a1
    
    Rayleigh_mass = Rayleigh_mass;
    Rayleigh_stiffness = Rayleigh_stiffness;
    
    MR = DeforModelRed3D
    fens,fes  = H8block(L, L, t, nL, nL, nt)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(FCplxFlt, size(fens.xyz,1), 3)) # displacement field
    nl = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    nl = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    nl = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    nl = selectnode(fens, box=[-Inf Inf L L -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    applyebc!(u)
    numberdofs!(u)
    # println("nfreedofs = $(u.nfreedofs)")
    
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)), material)
    femm = associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
    M = mass(femm, geom, u)
    C = Rayleigh_mass*M + Rayleigh_stiffness*K
    bfes = meshboundary(fes)
    sfemm = FEMMDeforSurfaceDamping(IntegDomain(bfes, GaussRule(2, 3)))
    impedance = rho_fluid * c_fluid
    D = dampingABC(sfemm, geom, u, impedance, SurfaceNormal(3))
    
    # if true
    #     t0 = time()
    #     d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    #     d = d - OmegaShift;
    #     fs = real(sqrt.(complex(d)))/(2*pi)
    #     println("Reference Eigenvalues: $fs [Hz]")
    #     println("eigs solution ($(time() - t0) sec)")
    # end
    
    bdryfes = meshboundary(fes)
    topbfl = selectelem(fens, bdryfes, facing=true, direction=[0.0 0.0 1.0])
    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,topbfl), GaussRule(2,2)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, feid::FInt, qpid::FInt) where {T}
        forceout .=  [0.0, 0.0, -qmagn]
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F = distribloads(el1femm, geom, u, fi, 2);
    

    K_ff = matrix_blocked(K, nfreedofs(u), nfreedofs(u))[:ff]
    M_ff = matrix_blocked(M, nfreedofs(u), nfreedofs(u))[:ff]
    D_ff = matrix_blocked(D, nfreedofs(u), nfreedofs(u))[:ff]
    C_ff = matrix_blocked(C, nfreedofs(u), nfreedofs(u))[:ff]

    F_f = vector_blocked(F, nfreedofs(u))[:f]
    U_d = gathersysvec(u, :d)

    U1 = zeros(FCplxFlt, nfreedofs(u), length(frequencies))
    for k = 1:length(frequencies)
        frequency = frequencies[k];
        omega = 2*pi*frequency;
        U1[:, k] = (-omega^2*M_ff + 1im*omega*(C_ff+D_ff) + K_ff)\F_f;
    end
    
    midpoint = selectnode(fens, box=[L/2 L/2 L/2 L/2 0 0], inflate=tolerance);
    midpointdof = u.dofnums[midpoint, 3]
    
    umidAmpl = abs.(U1[midpointdof, :])/phun("MM")
    @test abs(maximum(umidAmpl) - 9.56807e+00) < 1.0e-3 * 9.56807e+00

    umidReal = real.(U1[midpointdof, :])/phun("MM")
    umidImag = imag.(U1[midpointdof, :])/phun("MM")
    
    
    umidPhase = atan.(umidImag,umidReal)/pi*180 
    
 
    true
end # TEST13H_hva

TEST13H_hva()

end # module 


module mdistorted_block_infsup1
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
# using UnicodePlots


function distorted_block_infsup_T4()
    lambdatol = sqrt(1e8*eps(1.0));
    E=1000.0;
    nu=0.24;
    parshiftmult= 0.002;
    A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

    lambdamin = Float64[]
    h = Float64[]
    for ne = [2, 3, 4]
        Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

        fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
        # fens, fes = T4toT10(fens, fes)
        # @show connasarray(fes)

        for i = 1:count(fens)
            fens.xyz[i,:] = fens.xyz[i,:] + vec(reshape(fens.xyz[i,:], 1, 3)*A);
        end
        # @show fens.xyz

        # File =  "minfsuptest1.vtk"
        # vtkexportmesh(File, fens, fes)
        # try rm(File); catch end

        MR  =  DeforModelRed3D

        material = MatDeforElastIso(MR, E, nu)


        geom = NodalField(fens.xyz)
        u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
        bfes = meshboundary(fes)
        l1 = connectednodes(bfes)
        setebc!(u, l1, true, 1, 0.0)
        setebc!(u, l1, true, 2, 0.0)
        setebc!(u, l1, true, 3, 0.0)
        numberdofs!(u)

        femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
        Gh = infsup_gh(femm, geom, u);
        femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
        Sh = infsup_sh(femm, geom, u);


        G_ff = matrix_blocked(Gh, nfreedofs(u), nfreedofs(u))[:ff]
        S_ff = matrix_blocked(Sh, nfreedofs(u), nfreedofs(u))[:ff]

        lambda, modes = eigen(Matrix(G_ff), Matrix(S_ff));

        # @show lambda
        abslambda = real.(filter(y -> !isnan(y), lambda));
        ix = findall(y  -> y < 0.0, abslambda);
        if !isempty(ix)
            abslambda[ix] .= 0;
        end

        abslambda = sqrt.(sort(abslambda));
        ix = findall(y  -> y > 0.0, abslambda);
        # a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
        # display(a)

        ix = findall(y  -> y >= lambdatol, abslambda);
        if isempty(ix)
            @error "Bad guess of the number of eigenvalues"
        end
        push!(lambdamin, abslambda[ix[1]])
        push!(h, 1.0/(count(fens))^(1/3))
    end

    # @show lambdamin
    @test norm(lambdamin - [2.70777e-01, 1.79116e-01, 1.32893e-01]) < 1.0e-3
    # a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
    # display(a)

    # @test norm(lambdamin - [0.262065, 0.1709, 0.126159, 0.100228, 0.0828139]) / norm(lambdamin) <= 1.0e-4

end

distorted_block_infsup_T4()

end # module 
