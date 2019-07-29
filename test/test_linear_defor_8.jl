module mmmtdt1
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
function test()
	MR = DeforModelRed3D
	symmet(a) = a + transpose(a)
	
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6

	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strain3x3tto6v!(av, a)
	@test abs(det(a) - strain6vdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strain3x3tto6v!(av, a)
	@test abs(det(a) - strain6vdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strain3x3tto6v!(av, a)
	@test abs(det(a) - strain6vdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strain3x3tto6v!(av, a)
	@test abs(det(a) - strain6vdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strain3x3tto6v!(av, a)
	@test abs(det(a) - strain6vdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strain3x3tto6v!(av, a)
	@test abs(det(a) - strain6vdet(MR, av)) / abs(det(a)) <= 1.0e-6

	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strain3x3tto6v!(av, a)
	@test abs(tr(a) - strain6vtr(MR, av)) / abs(tr(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strain3x3tto6v!(av, a)
	@test abs(tr(a) - strain6vtr(MR, av)) / abs(tr(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strain3x3tto6v!(av, a)
	@test abs(tr(a) - strain6vtr(MR, av)) / abs(tr(a)) <= 1.0e-6
end
end
using .mmmtdt1
mmmtdt1.test()
