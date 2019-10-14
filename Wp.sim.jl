

include("Wp.jl")
using Random,Distributions

Random.seed!(42);
Dat = rand(Normal(),8,10);
Dat[1:4,1:5] = Dat[1:4,1:5] .*2 .+ 100;
Dat[1:2,6:8] = Dat[1:2,6:8] .* 2.5 .+ 100;
Dat[5:6,9:10] = Dat[5:6,9:10] .* 2 .+ 100;

res = Wp(Dat,1)
res[1]o




Data =[
			 1 2 3 4 5 5 5 4 0 0 0 ;
			 1 1 2 4 4 5 4 5 0 0 0 ;
			 1 2 3 5 5 0 0 0 0 0 0 ;
			 1 1 2 4 5 0 0 0 0 0 0 ;
			 0 0 0 0 0 0 0 0 0 2 5 ;
			 0 0 0 0 0 0 0 0 1 1 5 ;
			 0 0 0 0 0 0 0 0 0 0 0 ;
			 0 0 0 0 0 0 0 0 0 0 1. ]
Noise = rand(9,11)



include("Wp.jl")
U,Z,Nk,W=Wp(Data,1.5,RequiredK=3)
U,Z,Nk,W=Wp(Data+Noise,1.5,RequiredK=1)

U
Z
Nk
W




1  1  1  1  1  1  1
2  2  2  2  2  2  2
3  3  3  3  2  2  2
4  4  4  4  4  2  2
5  5  5  4  4  2  2
6  6  6  6  6  6  6
7  7  6  6  6  6  6
8  8  8  8  8  8  6
