module distorted_block_infsup_examples
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using FinEtools.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
using UnicodePlots


function distorted_block_infsup_T10()
	lambdatol = sqrt(1e8*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [2, 4, 6]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		fens, fes = T4toT10(fens, fes)
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

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	@show lambdamin
	a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	display(a)

	# @test norm(lambdamin - [0.262065, 0.1709, 0.126159, 0.100228, 0.0828139]) / norm(lambdamin) <= 1.0e-4

end

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

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	@show lambdamin
	a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	display(a)

	# @test norm(lambdamin - [0.262065, 0.1709, 0.126159, 0.100228, 0.0828139]) / norm(lambdamin) <= 1.0e-4

end

function distorted_block_infsup_H8()
	lambdatol = sqrt(1e8*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [2, 3, 4]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = H8block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt)
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

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)
		Gh = infsup_gh(femm, geom, u);
		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)
		Sh = infsup_sh(femm, geom, u);

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	@show lambdamin
	a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	display(a)

	# @test norm(lambdamin - [0.262065, 0.1709, 0.126159, 0.100228, 0.0828139]) / norm(lambdamin) <= 1.0e-4

end

function distorted_block_infsup_H20()
	lambdatol = sqrt(1e8*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [2, 3, 4]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = H20block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt)
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

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
		Gh = infsup_gh(femm, geom, u);
		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
		Sh = infsup_sh(femm, geom, u);

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	@show lambdamin
	a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	display(a)

	# @test norm(lambdamin - [0.262065, 0.1709, 0.126159, 0.100228, 0.0828139]) / norm(lambdamin) <= 1.0e-4

end

function allrun()
    println("#####################################################")
    println("# distorted_block_infsup_T4 ")
    distorted_block_infsup_T4()
    println("#####################################################")
    println("# distorted_block_infsup_H8 ")
    distorted_block_infsup_H8()
    println("#####################################################")
    println("# distorted_block_infsup_T10 ")
    distorted_block_infsup_T10()
    println("#####################################################")
    println("# distorted_block_infsup_H20 ")
    distorted_block_infsup_H20()
    return true
end # function allrun

end # module distorted_block_infsup_examples
