module m4thsymmcheck1
using FinEtools
using FinEtoolsDeforLinear: tens4checksymmetry, tens4symm6x6tot!, tens4symmtto6x6t!
using LinearAlgebra
using Test
function test()
	C = rand(6, 6)
	C = C + C'
	Co = deepcopy(C)
	t = fill(0.0, 3, 3, 3, 3)
	tens4symm6x6tot!(t, C)
	@test tens4checksymmetry(t)
	tens4symmtto6x6t!(C, t)
	@assert norm(C - C') <= eps(1.0)
	@assert norm(C - Co) <= eps(1.0)
end
end
using .m4thsymmcheck1
m4thsymmcheck1.test()