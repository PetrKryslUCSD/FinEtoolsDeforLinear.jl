using SparseArrays
using BandedMatrices
using BenchmarkTools

n=20000;
D=rand(n,);
D1=rand(n-1,);
Dm1=rand(n-1,);
D2=rand(n-2,);
Dm2=rand(n-2,);
d0=Pair(0,D); d1=Pair(1,D1); d2=Pair(2,D2);
dm1=Pair(-1,Dm1); dm2=Pair(-2,Dm2);
FVS=spdiagm(dm2,dm1,d0,d1,d2);
FVB=BandedMatrix(dm2,dm1,d0,d1,d2);
b=rand(n,);

@btime $FVB=BandedMatrix($FVS)
@btime $FVB\$b;
@btime $FVS\$b;
