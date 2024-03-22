module stubby_corbel_examples
using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy, solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule:
    linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra
using Krylov, LinearOperators, IncompleteLU
using SparseArrays
using SuiteSparse
using Printf
using SymRCM
using UnicodePlots
using PlotlyJS
# using Infiltrator
using Random
using ILUZero
# using SkylineSolvers
using LDLFactorizations
using LimitedLDLFactorizations, LinearOperators, Krylov
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: domain_decomposition, parallel_matrix_assembly!, SysmatAssemblerSparsePatt
using DataDrop
import CoNCMOR: CoNCData, transfmatrix, LegendreBasis, SineCosineBasis

# Isotropic material
E = 1000.0
nu = 0.4999 # Taylor data: nearly incompressible material
nu = 0.3 # Compressible material
W = 25.0
H = 50.0
L = 50.0
htol = minimum([L, H, W]) / 1000
uzex = -0.16
magn = 0.2 * (-12.6) / 4
Force = magn * W * H * 2
CTE = 0.0
n = 5 #

function getfrcL!(forceout, XYZ, tangents, feid, qpid)
    copyto!(forceout, [0.0; 0.0; magn])
end

transfm(m, t, tT) = (tT * m * t)
transfv(v, t, tT) = (tT * v)
function morprecon2d(y, v, Phi, Kd, K, mixprop)
    Kr = transfm(K, Phi, Phi')
    y .= (1 - mixprop) * (Kd \ v) + mixprop * (Phi * (Kr \ (Phi' * v)))
    return y
end
function morprecond3(y, v, Phi, Kd, Kr, mixprop)
    y .= (1 - mixprop) * (Kd \ v) + mixprop * (Phi * (Kr \ (Phi' * v)))
    return y
end
function morprecond4(y, v, Phi, Kd, Kr, mixprop)
    y .= (Phi * (((1 - mixprop) * (Phi' * Kd * Phi) + mixprop * Kr) \ (Phi' * v)))
    return y
end
# morprecond3nomix(y, v, Phi, Kd, Kr, mixprop) = begin
#     y .= (Kd \ v) + (Phi * (Kr \ (Phi' * v)))
#     return y
# end

function ildl(K::AbstractMatrix{T}; kwargs...) where {T}
    F = lldl(K; kwargs...)
    F.D .= abs.(F.D)
    n = length(F.D)
    return LinearOperator(T, n, n, true, true, (y, v) -> ldiv!(y, F, v))
end

function stubby_corbel_H8_by_hand()
    elementtag = "H8"
    println("""
    Stubby corbel example. Element: $(elementtag)
    """)

    fens, fes = H8block(W, L, H, n, 2 * n, 2 * n)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(csmatout, XYZ, tangents, feid, qpid)
        copyto!(csmatout, csmat)
    end

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)
    println("nfreedofs(u) = $(nfreedofs(u))")

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F = distribloads(el2femm, geom, u, fi, 2)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    @time solve_blocked!(u, K, F)
    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)
    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "stubby_corbel_H8_by_hand.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
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

function stubby_corbel_H8_big_iso(n = 10, solver = :suitesparse)
    elementtag = "H8"
    println("""
    Stubby corbel example. Element: $(elementtag)
    """)

    fens, fes = H8block(W, L, H, n, 2 * n, 2 * n)
    println("Number of elements: $(count(fes))")
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(csmatout, XYZ, tangents, feid, qpid)
        copyto!(csmatout, csmat)
    end

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    C = connectionmatrix(femm, count(fens))
    # display(spy(C, canvas = DotCanvas))
    # I, J, V = findnz(C)
    # @show bw = maximum(I .- J) + 1
    perm = symrcm(C)
    # C = C[perm, perm]
    # display(spy(C, canvas = DotCanvas))
    # I, J, V = findnz(C)
    # @show bw = maximum(I .- J) + 1 

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u, perm)
    # numberdofs!(u)
    println("nfreedofs(u) = $(nfreedofs(u))")

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F = distribloads(el2femm, geom, u, fi, 2)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    println("Stiffness: number of non zeros = $(nnz(K)) [ND]")
    println("Sparsity = $(nnz(K)/size(K, 1)/size(K, 2)) [ND]")
    display(spy(K, canvas = DotCanvas))
    # DataDrop.store_matrix("K", K)

    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)

    if solver == :suitesparse
        # @show methods(SuiteSparse.CHOLMOD.ldlt, (typeof(K), ))
        # @time K = SuiteSparse.CHOLMOD.ldlt(K)
        @time K = SuiteSparse.CHOLMOD.cholesky(K)
        # @time K = SparseArrays.ldlt(K)
        # @time K = cholesky(K)
        @time U = K \ (F)
    elseif solver == :cg
        n = size(K, 1)
        mKd = mean(diag(K))
        # @time factor = ilu(K, τ = mKd / 100.0) # This may work for compressible materials
        @time factor = ilu(K, τ = mKd / 1000000.0) # This may work for incompressible materials
        # factor = ilu0(K)
        @show nnz(factor) / nnz(K)
        opM = LinearOperator(Float64, n, n, false, false, (y, v) -> ldiv!(y, factor, v))
        @time (U, stats) = Krylov.cg(K, F; M = opM, itmax = Int(round(n / 2)), verbose = 1)
    elseif solver == :scaledcg
        n = size(K, 1)
        idKs = Diagonal(1.0 ./ sqrt.(diag(K)))
        sK = idKs * K * idKs
        @show mKd = mean(diag(sK))
        # @time factor = ilu(sK, τ = 0.01) # This may work for compressible materials
        @time factor = ilu(sK, τ = 0.000001) # This may work for incompressible materials
        # @time factor = ilu0(sK)
        @show nnz(factor) / nnz(K)
        opM = LinearOperator(Float64, n, n, false, false, (y, v) -> ldiv!(y, factor, v))
        @time (U, stats) =
            Krylov.cg(sK, idKs * F; M = opM, itmax = Int(round(n / 2)), verbose = 1)
        U = Vector(idKs * U)
    elseif solver == :skyline
        I, J, V = findnz(K)
        @show bw = maximum(abs.(I .- J)) + 1
        M = size(K, 1)
        K = nothing
        GC.gc()
        sky = SkylineSolvers.Ldlt.SkylineMatrix(I, J, V, M)
        I = nothing
        J = nothing
        V = nothing
        GC.gc()
        @show SkylineSolvers.Ldlt.nnz(sky)
        @time SkylineSolvers.Ldlt.factorize!(sky)
        @time U = SkylineSolvers.Ldlt.solve(sky, F)
    elseif solver == :mor
        Nc = 16
        nbf1max = 3
        mixprop = 0.5
        partitioning = nodepartitioning(fens, Nc)
        mor = CoNCData(fens, partitioning)
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u)
        transfm(m, t, tT) = (tT * m * t)
        transfv(v, t, tT) = (tT * v)
        PhiT = Phi'
        Kr = transfm(K, Phi, PhiT)
        @show size(Kr)
        Krfactor = lu(Kr)
        Ur = Phi * (Krfactor \ (PhiT * F))
        scattersysvec!(u, Ur[:])
        utip = mean(u.values[Tipl, 3], dims = 1)
        println("First Guess of Deflection: $(utip), compared to $(uzex)")
        R = F - K * Ur
        n = size(K, 1)
        Kd = Diagonal(diag(K))
        Krd = fill(0.0, size(K, 1))
        for i = 1:size(K, 1)
            Krd = Phi[i, :] * Kr * PhiT[i, :]
        end
        trace1 = scatter(
            x = 1:length(Krd),
            y = Krd,
            mode = "lines",
            line_width = 1.5,
            line_color = "RoyalBlue",
        )
        trace2 = scatter(
            x = 1:length(Krd),
            y = diag(K),
            mode = "lines",
            line_width = 1.5,
            line_color = "red",
        )
        data = [trace1, trace2]
        layout = Layout(; title = "Data Labels Hover")

        display(plot(data, layout))
        # @show  mKd = mean(diag(K))
        # @time factor0 = ilu(K, τ = mKd / 10.0) # This may work for incompressible materials
        # itr = 0
        function morprecond(y, v)
            y .= (1 - mixprop) * (Kd \ v) + mixprop * (Phi * (Krfactor \ (PhiT * v)))
        end
        opM = LinearOperator(Float64, n, n, false, false, morprecond)
        @time (DU, stats) = Krylov.cg(K, R; M = opM, itmax = 500, verbose = 1)
        U = Ur + DU
    else
        @error "Solver not recognized"
    end
    scattersysvec!(u, U[:])

    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "stubby_corbel_H8_big.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # stubby_corbel_H8_big_iso

function _cg(A, b, x0, maxiter)
    x = deepcopy(x0)
    g = x' * A - b'
    d = -g'
    for iter in 1:maxiter
        Ad = A * d
        rho = (d' * Ad)
        alpha = (-g * d) / rho
        x = x + alpha * d
        g = x' * A - b'
        beta = (g * Ad) / rho
        d = beta * d - g'
    end
    return x
end

function stubby_corbel_H8_big_ms(n = 10, solver = :suitesparse)
    elementtag = "H8"
    println("""
    Stubby corbel example. Element: $(elementtag)
    """)

    fens, fes = H8block(W, L, H, n, 2 * n, 2 * n)
    println("Number of elements: $(count(fes))")
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(csmatout, XYZ, tangents, feid, qpid)
        copyto!(csmatout, csmat)
    end

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    C = connectionmatrix(femm, count(fens))
    # display(spy(C, canvas = DotCanvas))
    # I, J, V = findnz(C)
    # @show bw = maximum(I .- J) + 1
    perm = symrcm(C)
    # C = C[perm, perm]
    # display(spy(C, canvas = DotCanvas))
    # I, J, V = findnz(C)
    # @show bw = maximum(I .- J) + 1 

    # femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u, perm)
    # numberdofs!(u)
    println("nfreedofs(u) = $(nfreedofs(u))")

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F = distribloads(el2femm, geom, u, fi, 2)
    F_f = vector_blocked_f(F, nfreedofs(u))
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K_ff = matrix_blocked_ff(K, nfreedofs(u))
    K = nothing
    println("Stiffness: number of non zeros = $(nnz(K_ff)) [ND]")
    println("Sparsity = $(nnz(K_ff)/size(K_ff, 1)/size(K_ff, 2)) [ND]")
    display(spy(K_ff, canvas = DotCanvas))
    # DataDrop.store_matrix("K$(size(K, 1))", K)
    # DataDrop.store_matrix("F$(size(F, 1))", F)

    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)

    if solver == :suitesparse || solver == :default
        # @show methods(SuiteSparse.CHOLMOD.ldlt, (typeof(K), ))
        # @time K = SuiteSparse.CHOLMOD.ldlt(K)
        @time K = SuiteSparse.CHOLMOD.cholesky(K)
        @show nnz(K)
        # @time K = SparseArrays.ldlt(K)
        # @time K = cholesky(K)
        @time U = K \ (F)
    elseif solver == :cg
        n = size(K_ff, 1)
        mK_ffd = mean(diag(K_ff))
        @time factor = ilu(K_ff, τ = mK_ffd / 100.0) # This may work for compressible materials
        # @time factor = ilu(K, τ = mKd / 1000000.0) # This may work for incompressible materials
        # factor = ilu0(K)
        @show nnz(factor) / nnz(K_ff)
        opM = LinearOperator(Float64, n, n, false, false, (y, v) -> ldiv!(y, factor, v))
        @time (U, stats) = Krylov.cg(K_ff, F_f; M = opM, itmax = Int(round(n / 2)), verbose = 1)
    elseif solver == :cgldl
        n = size(K, 1)
        atol = 1e-10
        rtol = 0.0
        memory = 2000
        @time P = ildl(K, memory = memory)
        # @time U, stats = bicgstab(K, F, N=P, atol=atol, rtol=rtol, verbose=1)
        @time U, stats = cg(K, F, M = P, atol = atol, rtol = rtol, verbose = 1)
        @show stats
    elseif solver == :scaledcg
        n = size(K, 1)
        idKs = Diagonal(1.0 ./ sqrt.(diag(K)))
        sK = idKs * K * idKs
        @show mKd = mean(diag(sK))
        # @time factor = ilu(sK, τ = 0.01) # This may work for compressible materials
        @time factor = ilu(sK, τ = 0.000001) # This may work for incompressible materials
        # @time factor = ilu0(sK)
        @show nnz(factor) / nnz(K)
        opM = LinearOperator(Float64, n, n, false, false, (y, v) -> ldiv!(y, factor, v))
        @time (U, stats) =
            Krylov.cg(sK, idKs * F; M = opM, itmax = Int(round(n / 2)), verbose = 1)
        U = Vector(idKs * U)
    elseif solver == :skyline
        I, J, V = findnz(K)
        @show bw = maximum(abs.(I .- J)) + 1
        M = size(K, 1)
        K = nothing
        GC.gc()
        sky = SkylineSolvers.Ldlt.SkylineMatrix(I, J, V, M)
        I = nothing
        J = nothing
        V = nothing
        GC.gc()
        @show SkylineSolvers.Ldlt.nnz(sky)
        @time SkylineSolvers.Ldlt.factorize!(sky)
        @time U = SkylineSolvers.Ldlt.solve(sky, F)
    elseif solver == :ldlfactorizations
        @time factors = LDLFactorizations.ldlt(K)
        @time U = factors \ F
    elseif solver == :mor0
        Nc = 32
        nbf1max = 4
        mixprop = 0.01
        partitioning = nodepartitioning(fens, Nc)
        mor = CoNCData(fens, partitioning)
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u)
        PhiT = Phi'
        Kr = transfm(K, Phi, PhiT)
        @show size(Kr)
        Krfactor = lu(Kr)
        Ur = Phi * (Krfactor \ (PhiT * F))
        scattersysvec!(u, Ur[:])
        utip = mean(u.values[Tipl, 3], dims = 1)
        println("First Guess of Deflection: $(utip), compared to $(uzex)")
        R = F - K * Ur
        n = size(K, 1)
        Kdinv = 1.0 ./ diag(K)
        morprecond(y, v) = begin
            y .= (1 - mixprop) * (Kdinv .* v)
            y .+= mixprop * (Phi * (Krfactor \ (PhiT * v)))
        end
        morprecondnomix(y, v) = begin
            y .= (Kdinv .* v) + (Phi * (Krfactor \ (PhiT * v)))
        end
        opM = LinearOperator(Float64, n, n, false, false, morprecond)
        U = deepcopy(Ur)
        for iter = 1:50
            @show iter
            (DU, stats) = Krylov.cg(K, F - K * U; M = opM, itmax = 5, verbose = 0)
            U += DU
            @show norm(DU) / norm(U)
            scattersysvec!(u, U[:])
            utip = mean(u.values[Tipl, 3], dims = 1)
            println("Iterated Deflection: $(utip), compared to $(uzex)")
        end
    elseif solver == :mor1
        rtol = 1.0e-9
        Nc = 32
        nbf1max = 4
        mixprop = 1.0
        partitioning = nodepartitioning(fens, Nc)
        mor = CoNCData(fens, partitioning)
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u)
        PhiT = Phi'
        Kr = transfm(K, Phi, PhiT)
        @show size(Kr)
        Krfactor = lu(Kr)
        Kdinv = 1.0 ./ diag(K)
        # invKrd = fill(0.0, size(K, 1))
        # for i in 1:size(K, 1)
        #     invKrd[i] = @views dot(vec(Phi[i, :]), Krfactor \ Vector(Phi[i, :]))
        # end
        # @show norm(invKrd), norm(Kdinv)
        # Kdinv .*= norm(invKrd) / norm(Kdinv)
        # trace1 = scatter(x = 1:length(invKrd),  
        #     y = invKrd,
        #     mode="points",
        #     line_width=1.5,
        #     line_color="RoyalBlue")
        # trace2 = scatter(x = 1:length(invKrd),  
        #     y = Kdinv,
        #     mode="lines",
        #     line_width=1.5,
        #     line_color="red")
        # data = [trace1, trace2]
        # # data = [trace2]
        # layout = Layout(;title="Diagonals")
        # display(plot(data, layout))
        Ur = Phi * (Krfactor \ (PhiT * F))
        scattersysvec!(u, Ur[:])
        utip = mean(u.values[Tipl, 3], dims = 1)
        println("First Guess of Deflection: $(utip), compared to $(uzex)")
        R = F - K * Ur
        n = size(K, 1)
        # morprecondnomix(y, v) = begin
        #     y .= (Kdinv .* v) - (invKrd .* v) + (Phi * (Krfactor \ (PhiT * v)))
        # end
        opM = LinearOperator(
            Float64,
            n,
            n,
            false,
            false,
            (y, v) -> y .= mixprop .* (Kdinv .* v) .+ (Phi * (Krfactor \ (PhiT * v))),
        )
        U = deepcopy(Ur)
        utipprev = utip
        @time for iter = 1:50
            @show iter
            (DU, stats) = Krylov.cg(K, F - K * U; M = opM, itmax = 5, verbose = 0)
            U += DU
            @show norm(DU) / norm(U)
            scattersysvec!(u, U[:])
            utip = mean(u.values[Tipl, 3], dims = 1)
            println("Iterated Deflection: $(utip), compared to $(uzex)")
            if norm(utip - utipprev) / norm(utip) < rtol
                break
            end
            utipprev = utip
        end
    elseif solver == :mor
        Nc = 16
        nbf1max = 4
        mixprop = 0.05
        partitioning = nodepartitioning(fens, Nc)
        mor = CoNCData(fens, partitioning)
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u)
        PhiT = Phi'
        Kr = transfm(K, Phi, PhiT)
        @show size(Kr)
        Ur = Phi * (Kr \ (PhiT * F))
        scattersysvec!(u, Ur[:])
        utip = mean(u.values[Tipl, 3], dims = 1)
        println("First Guess of Deflection: $(utip), compared to $(uzex)")
        Phi = hcat(Ur)
        U = deepcopy(Ur)
        Kd = Diagonal(diag(K))
        n = size(K, 1)

        for iter = 1:50
            @show iter
            Kr = transfm(K, Phi, Phi')
            # Krd = Diagonal(diag(Phi * Kr * Phi'))
            # @show mean(diag(Kd)), mean(diag(Krd))
            # @show norm(diag(Kd) - diag(Krd)), norm(diag(Kd))
            opM = LinearOperator(
                Float64,
                n,
                n,
                false,
                false,
                (y, v) -> morprecond3nomix(y, v, Phi, Kd, Kr, mixprop),
            )
            @time (DU, stats) = Krylov.cg(K, F - K * U; M = opM, itmax = 10, verbose = 0)
            @show norm(DU) / norm(U)
            Phi = hcat(Phi, DU)
            factors = qr(Phi)
            Phi = Matrix(factors.Q)
            U += DU
            scattersysvec!(u, U[:])
            utip = mean(u.values[Tipl, 3], dims = 1)
            println("Iterated Deflection: $(utip), compared to $(uzex)")
        end

    else
        @error "Solver not recognized"
    end
    scattersysvec!(u, U[:])

    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "stubby_corbel_H8_big_ms.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # stubby_corbel_H8_big_ms

function stubby_corbel_H8_big_ms_parallel(N = 10, 
    ntasks = Threads.nthreads(), assembly_only = false)
    elementtag = "MSH8"
    println("""
    Stubby corbel example. Element: $(elementtag)
    """)

    times = Dict{String, Vector{Float64}}()
    
    t1 = time()
    fens, fes = H8block(W, L, H, N, 2 * N, 2 * N)
    times["MeshGeneration"] = [time() - t1]
    println("Number of elements: $(count(fes))")
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(csmatout, XYZ, tangents, feid, qpid)
        copyto!(csmatout, csmat)
    end

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)
    # numberdofs!(u)
    println("nfreedofs(u) = $(nfreedofs(u))")

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F = distribloads(el2femm, geom, u, fi, 2)
    F_f = vector_blocked_f(F, nfreedofs(u))

    function createsubdomain(fessubset)
        FEMMDeforLinearMSH8(MR, IntegDomain(fessubset, GaussRule(3, 2)), material)
    end

    function matrixcomputation!(femm, assembler)
        associategeometry!(femm, geom)
        stiffness(femm, assembler, geom, u)
    end

    t1 = time()
    n2e = FENodeToFEMap(fes.conn, nnodes(u))
    times["FENodeToFEMap"] = [time() - t1]
    println("Make node to element map = $(times["FENodeToFEMap"]) [s]")

    println("Stiffness =============================================================")
    GC.enable(false)

    t0 = time(); 

    t1 = time()
    coloring = element_coloring(fes, n2e)
    times["ElementColors"] = [time() - t1]
    println("    Compute element colors = $(times["ElementColors"]) [s]")

    t1 = time()
    n2n = FENodeToNeighborsMap(n2e, fes)
    times["FENodeToNeighborsMap"] = [time() - t1]
    println("    Make node to neighbor map = $(times["FENodeToNeighborsMap"]) [s]")

    t1 = time()
    K_pattern = sparse_symmetric_csc_pattern(u.dofnums, nalldofs(u), n2n, zero(eltype(u.values)))
    times["SparsityPattern"] = [time() - t1]
    println("    Sparsity pattern = $(times["SparsityPattern"]) [s]")

    t1 = time()
    decomposition = domain_decomposition(fes, coloring, createsubdomain, ntasks)
    times["DomainDecomposition"] = [time() - t1]
    println("    Domain decomposition = $(times["DomainDecomposition"]) [s]")

    t1 = time()
    K = parallel_matrix_assembly!(
        SysmatAssemblerSparsePatt(K_pattern),
        decomposition,
        matrixcomputation!,
        ntasks
    )
    times["AssemblyOfValues"] = [time() - t1]
    println("    Add to matrix = $(times["AssemblyOfValues"]) [s]")

    times["TotalAssembly"] = [time() - t0]
    println("Assembly total = $(times["TotalAssembly"]) [s]")

    GC.enable(true)

    K_ff = matrix_blocked_ff(K, nfreedofs(u))
    F_f = vector_blocked_f(F, nfreedofs(u))
    println("Stiffness: number of non zeros = $(nnz(K_ff)) [ND]")
    println("Sparsity = $(nnz(K_ff)/size(K_ff, 1)/size(K_ff, 2)) [ND]")
    
    if assembly_only
        isdir("$(N)") || mkdir("$(N)")
        n = DataDrop.with_extension(joinpath("$(N)", "stubby_corbel_H8_big_ms_parallel-timing-nth=$(ntasks)"), "json")
        if isfile(n)
            storedtimes = DataDrop.retrieve_json(n)
            for k in keys(storedtimes)
                times[k] = cat(times[k], storedtimes[k], dims = 1)
            end
        end
        DataDrop.store_json(n, times)
        return
    end
    
    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)
    @show norm(K_ff - K_ff') / norm(K_ff)
    @time K_ff_factors = SuiteSparse.CHOLMOD.cholesky(Symmetric(K_ff))
    @show nnz(K_ff_factors)
    # @time K = SparseArrays.ldlt(K)
    # @time K = cholesky(K)
    @time U_f = K_ff_factors \ (F_f)

    scattersysvec!(u, U_f[:])

    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "stubby_corbel_H8_big_ms.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    # @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # stubby_corbel_H8_big_ms

function stubby_corbel_H8_big_ms_serial(N = 10, assembly_only = false)
    elementtag = "MSH8"
    println("""
    Stubby corbel example. Element: $(elementtag). SERIAL
    """)

    times = Dict{String, Vector{Float64}}()
    
    t1 = time()
    fens, fes = H8block(W, L, H, N, 2 * N, 2 * N)
    times["MeshGeneration"] = [time() - t1]
    println("Number of elements: $(count(fes))")
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(csmatout, XYZ, tangents, feid, qpid)
        copyto!(csmatout, csmat)
    end

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)
    # numberdofs!(u)
    println("nfreedofs(u) = $(nfreedofs(u))")

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F = distribloads(el2femm, geom, u, fi, 2)
    F_f = vector_blocked_f(F, nfreedofs(u))

    # femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    t0 = time()
    t1 = time()
    associategeometry!(femm, geom)
    assembler = SysmatAssemblerSparse(1.0)  
    setnomatrixresult(assembler, true)
    stiffness(femm, assembler, geom, u)
    times["ComputeCOO"] = [time() - t1]
    t1 = time()
    setnomatrixresult(assembler, false)
    K = makematrix!(assembler)
    times["BuildCSR"] = [time() - t1]

    times["TotalAssembly"] = [time() - t0]
    println("Assembly total = $(times["TotalAssembly"]) [s]")

    K_ff = matrix_blocked_ff(K, nfreedofs(u))
    F_f = vector_blocked_f(F, nfreedofs(u))
    println("Stiffness: number of non zeros = $(nnz(K_ff)) [ND]")
    println("Sparsity = $(nnz(K_ff)/size(K_ff, 1)/size(K_ff, 2)) [ND]")
    
    if assembly_only
        isdir("$(N)") || mkdir("$(N)")
        n = DataDrop.with_extension(joinpath("$(N)", "stubby_corbel_H8_big_ms_parallel-timing-serial"), "json")
        if isfile(n)
            storedtimes = DataDrop.retrieve_json(n)
            for k in keys(storedtimes)
                times[k] = cat(times[k], storedtimes[k], dims = 1)
            end
        end
        DataDrop.store_json(n, times)
        return
    end
    
    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)
    @show norm(K_ff - K_ff') / norm(K_ff)
    @time K_ff_factors = SuiteSparse.CHOLMOD.cholesky(Symmetric(K_ff))
    @show nnz(K_ff_factors)
    # @time K = SparseArrays.ldlt(K)
    # @time K = cholesky(K)
    @time U_f = K_ff_factors \ (F_f)

    scattersysvec!(u, U_f[:])

    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "stubby_corbel_H8_big_ms.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    # @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # stubby_corbel_H8_big_ms

function allrun()
    println("#####################################################")
    println("# stubby_corbel_H8_by_hand ")
    stubby_corbel_H8_by_hand()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
