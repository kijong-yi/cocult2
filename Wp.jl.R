 c(N,M) %<-% dim(Data)

library(XRJulia)
findJulia(path="/home/users/kjyi/bin/julia")
using("Statistics")


.Wp <- jf("function Wp(Data,p;RequiredK=1.,trace=false)
	# Data = data set (sample * feature) must be already standardized
	# p = Minkowski Exponent p. Good values in [1:0.1:5]
  # RequiredK=1 for full hierarchy
  Data = float(Data)
  # RequiredK=int(RequiredK)
	N,M = size(Data) # N=sample,M=feature
	#U = collect(1:N)' # cluster information, starts from singleton vectors
	U = ones(Int,N,N-RequiredK+1)
	U[:,1] = 1:N
	K = N
	Nk = ones(N,1) # cluster size
	Z = Data[:,:] # centroids. located in certain rows, otherwise marked as Inf.
	W = fill(1/M,K,M) # final Weights will be calculated
	WL = [] # tracing W
	InitialK=K
	AllDistances=fill(Inf,InitialK,InitialK)
	for k1 in 1:InitialK-1, k2 in k1+1:InitialK
		# minkowsky distances of all clusters
		AllDistances[k1,k2] = 0.5 * sum( (abs.(Z[k1,:]-Z[k2,:]).^p) .* (W[1,:].^p) ) 
		AllDistances[k2,k1] = AllDistances[k1,k2]
	end	
	UIndex=1
	while K > RequiredK 
	# 1st steep: look for clusters to merge
		k2Min,k1Min = argmin(AllDistances).I
		# merge k1Min and k2Min
		UIndex += 1
		U[:,UIndex],Nk,Z,K=Merge(k1Min,k2Min,U[:,UIndex-1],Data,Nk,M,Z,p,K)
		W = GetNewW(Data,W,U[:,UIndex],Z,InitialK,M,p)
		if trace
		  push!(WL,W)
		end
		AllDistances=Update_AllDistances(Z,W,p,Nk,AllDistances,k1Min,k2Min,InitialK)
	end
	ZZ=Z[:,:]
	ZZ[ZZ.==Inf].=0
	if trace
	  [U,ZZ,Nk,W,WL]
	else
	  [U,ZZ,Nk,W]
	end
	#U: contains the labels. 
	#Format: Each column is a merge iteration. In U[:,1] each entity will have
	#a different cluster index since each onf them is a different cluster. In
	#U[:,2] two of these clusters will be merged. In U[:,3] Three clusters will
	#be merged, etc.
	#Final size of U: NumberOfEntities x (NumberOfEntities-requiredK+1)
	#Z: This contains the final centroids. Note that the centroids will be
	#located in certain rows only, the others will be marked with inf.
	#Nk: Number of entities per cluster. The index of the numbers will match
	#that of the centroids in Z.
	#W: Final weightes. You should check the weights matching whose row numbers
	#match those of the centroids in Z.
end")

juliaCommand("
function Update_AllDistances(Z,W,p,Nk,AllDistances, k1Min, k2Min, InitialK)
		#remove k2 from Z
		AllDistances[k2Min,:].=Inf
		AllDistances[:,k2Min].=Inf
		#update the distances related to z1
		for k ∈ 1:InitialK
				if k==k1Min||Nk[k,1]==0||k==k2Min
					continue
				end
				avg_W = ((W[k1Min,:] + W[k,:])./2).^p
				AllDistances[k1Min,k] = ((Nk[k1Min] * Nk[k])/(Nk[k1Min] + Nk[k])) * sum((abs.(Z[k1Min,:] - Z[k,:]).^p).*avg_W)
				AllDistances[k,k1Min] = AllDistances[k1Min,k]
		end
		AllDistances
end")

juliaCommand("
function GetNewW(Data,W, Ui, Z, K, M,p)
		D=zeros(K,M)
		for l ∈ 1:K, j ∈ 1:M
			D[l,j] = sum(abs.(Data[Ui.==l,j].-Z[l,j]).^p)
		end
		D .+= 0.0001
		#Calculate the actual Weight for each column
		if p ≠ 1
				for l ∈ 1:K, j ∈ 1:M
					 tmp=D[l,j]
					 W[l,j]= 1/sum((tmp./D[l,:]).^(1/(p-1)))
				end
		else
				for l ∈ 1:K
						MinIndex = argmin(D[l,:])
						W[l,1:M] .= 0 #necessary to zero all others
						W[l,MinIndex] = 1
				end
		end 
		W
end")

juliaCommand("
function Merge(k1Min, k2Min, Up, Data, Nk, M, Z, p, K)
	Nk[k1Min] += Nk[k2Min]
	Up[Up.==k2Min] .= k1Min
	Nk[k2Min] = 0
	Z[k1Min,:] = New_cmt(Data[Up.==k1Min,:],p)
	Z[k2Min,:] .= Inf
	K -= 1
	Up, Nk, Z, K
end")

juliaCommand("
function New_cmt(Data,p) # works only when p > 1
	#Calculates the Minkowski center at a given p.
	#Data MUST BE EntityxFeatures and standardised.
	N,M=size(Data)
	if p==1
		DataCenter=median(Data,dims=1)
		return(DataCenter)
	elseif p==2
		DataCenter=mean(Data,dims=1)
		return(DataCenter)
	elseif N==1
		DataCenter=Data;
		return(DataCenter)
	end
	Gradient = fill(0.001,1,M)
	OnesIndex = ones(Int,N,1)
	DataCenter = sum(Data,dims=1)./N
	DistanceToDataCenter=sum(abs.(Data .- DataCenter).^p, dims=1)
	NewDataCenter=DataCenter+Gradient
	DistanceToNewDataCenter=sum(abs.(Data .- NewDataCenter).^p, dims=1)
	Gradient[1,DistanceToDataCenter[:] .< DistanceToNewDataCenter[:]] .*= -1
	while true
		NewDataCenter = DataCenter + Gradient
		DistanceToNewDataCenter=sum(abs.(Data .- NewDataCenter).^p,dims=1) 
		Gradient[1,DistanceToNewDataCenter[:] .>= DistanceToDataCenter[:]] *= 0.9 
		DataCenter[1,DistanceToNewDataCenter[:].<DistanceToDataCenter[:]]=NewDataCenter[1,DistanceToNewDataCenter[:].<DistanceToDataCenter[:]]
		DistanceToDataCenter[1,DistanceToNewDataCenter[:].<DistanceToDataCenter[:]]=DistanceToNewDataCenter[1,DistanceToNewDataCenter[:].<DistanceToDataCenter[:]]
		DataCenter
		abs.(Gradient)
		if all(abs.(Gradient).<0.0001) break end
	end
	DataCenter
end")

Wp <- function(Data,p=1.5,RequiredK=1L, trace = F){
  data = unname(as.matrix(Data))
  res =.Wp(data, p, RequiredK=as.integer(RequiredK),trace=trace)
  if(length(res)==4){
    names(res) = c("U","Z","Nk","W")
  } else {
    names(res) = c("U","Z","Nk","W","trace")
  }
  rownames(res[[1]]) = rownames(Data)
  rownames(res[[2]]) = rownames(Data)
  colnames(res[[2]]) = colnames(Data)
  rownames(res[[3]]) = rownames(Data)
  rownames(res[[4]]) = rownames(Data)
  colnames(res[[4]]) = colnames(Data)
  res
}
if(F){
  Data = j("[
  1 2 3 4 5 5 5 4 0 0 0 ;
  1 1 2 4 4 5 4 5 0 0 0 ;
  1 2 3 5 5 0 0 0 0 0 0 ;
  1 1 2 4 5 0 0 0 0 0 0 ;
  0 0 0 0 0 0 0 0 0 2 5 ;
  0 0 0 0 0 0 0 0 1 1 5 ;
  0 0 0 0 0 0 0 0 0 0 0 ;
  0 0 0 0 0 0 0 0 0 0 1 ]")
  Data = matrix(c(
    1, 2, 3, 4, 5, 5, 5, 4, 0, 0, 0,
    1, 1, 2, 4, 4, 5, 4, 5, 0, 0, 0,
    1, 2, 3, 5, 5, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 5,
    1, 1, 2, 4, 5, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 5,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
  ), ncol = 11, nrow = 8, byrow = T)
  rownames(Data) = letters[1:8]
  colnames(Data) = LETTERS[1:11]
  Data
  res = Wp(Data,p=1.5,trace=T)
  res2 = Wp_R(Data,p=1.5,RequiredK=T,trace=T)
  res
  set.seed(42);Data = matrix(runif(35000),nrow=35)
  res = Wp(Data,p=1.5,trace=T)
  res2 = Wp_R(Data,p=1.5,RequiredK=T,trace=T)
  system.time({res2 = Wp_R(Data,p=1.5,RequiredK=T,trace=T)})
  
} else {
  cat("
res = Wp_R(Data, p=1~5, RequiredK=1L, trace = F
  res$U
  res$Z
  res$NK
  res$W
  res$WL")
  
}



