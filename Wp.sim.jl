include("Wp.jl")

using DelimitedFiles
raw=readdlm("A.std2.t.tsv");

Data=raw[2:end,2:end];
Data=convert(Array{Float64,2}, Data);
typeof(Data),size(Data)

U,Z,Nk,W,WL=Wp(Data[1:10,1:10],2,RequiredK=10,trace=true)

@save "Wp.k1.p2.jld"