using Statistics

function Wp(Data,p;RequiredK=1,trace=false)
	# Data = data set (sample * feature) must be already standardized
	# p = Minkowski Exponent p. Good values in [1:0.1:5]
	# for full hierarchy, RequiredK=1
	Data=float(Data);
	N,M = size(Data) # N=sample,M=feature
	U = ones(Int,N,N-RequiredK+1)
	U[:,1] = 1:N
	K = N
	Nk = ones(N) # cluster size
	Z = Data[:,:] # centroids. located in certain rows, otherwise marked as Inf.
	W = fill(1/M,K,M) # final Weights
	WL = [] # trace W
	InitialK=K
	AllDistances=fill(Inf,InitialK,InitialK)
	for k1 ∈ 1:InitialK-1, k2 ∈ k1+1:InitialK
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
			push!(WL,W[:,:])
		end
		AllDistances=Update_AllDistances(Z,W,p,Nk,AllDistances,k1Min,k2Min,InitialK)
	end
	U,Z,Nk,W,WL
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
end

function Update_AllDistances(Z,W,p,Nk,AllDistances, k1Min, k2Min, InitialK)
	#remove k2 from Z
	AllDistances[k2Min,:].=Inf
	AllDistances[:,k2Min].=Inf
	#update the distances related to z1
	for k ∈ 1:InitialK
		if k==k1Min||Nk[k]==0||k==k2Min
			continue
		end
		avg_W = ((W[k1Min,:] + W[k,:])./2).^p
		@inbounds AllDistances[k1Min,k] = ((Nk[k1Min] * Nk[k])/(Nk[k1Min] + Nk[k])) * sum((abs.(Z[k1Min,:] - Z[k,:]).^p).*avg_W)
		AllDistances[k,k1Min] = AllDistances[k1Min,k]
	end
	AllDistances
end

function GetNewW(Data,W, Ui, Z, K, M,p)
	D=zeros(K,M)
	for l ∈ 1:K, j ∈ 1:M
		@inbounds D[l,j] = sum(abs.(Data[Ui.==l,j].-Z[l,j]).^p)
	end
	D .+= 0.0001
	#Calculate the actual Weight for each column
	if p ≠ 1
		for l ∈ 1:K, j ∈ 1:M
			#@inbounds W[l,j]= Z[l,j]*(1/sum((D[l,j]./D[l,:]).^(1/(p-1))))
			@inbounds W[l,j]= 1/sum((D[l,j]./D[l,:]).^(1/(p-1)))
		end
	else
		for l ∈ 1:K
			MinIndex = argmin(D[l,:])
			W[l,1:M] .= 0 #necessary to zero all others
			W[l,MinIndex] = 1
		end
	end 
	W
end


function Merge(k1Min, k2Min, Up, Data, Nk, M, Z, p, K)
	Nk[k1Min] += Nk[k2Min]
	Up[Up.==k2Min] .= k1Min
	Nk[k2Min] = 0
	Z[k1Min,:] = New_cmt(Data[Up.==k1Min,:],p)
	Z[k2Min,:] .= Inf
	K -= 1
	Up, Nk, Z, K
end

function New_cmt(data,p) # works only when p > 1
	#Calculates the Minkowski center at a given p.
	#Data MUST BE EntityxFeatures and standardised.
	N,M=size(data)
	if p==1
		DataCenter=median(data,dims=1)
		return(DataCenter)
	elseif p==2
		DataCenter=mean(data,dims=1)
		return(DataCenter)
	elseif N==1
		DataCenter=Data;
		return(DataCenter)
	end
	Gradient = fill(0.001,1,M)
	DataCenter = mean(data,dims=1)
	DistanceToDataCenter=sum(abs.(data .- DataCenter).^p, dims=1)
	NewDataCenter=DataCenter+Gradient
	DistanceToNewDataCenter=sum(abs.(data .- NewDataCenter).^p, dims=1)
	@inbounds Gradient[DistanceToDataCenter .< DistanceToNewDataCenter] .*= -1
	while true
		NewDataCenter = DataCenter + Gradient
		DistanceToNewDataCenter = similar(DistanceToDataCenter)
		DistanceToNewDataCenter=sum(abs.(data .- NewDataCenter).^p,dims=1) ###
		@inbounds Gradient[DistanceToNewDataCenter .>= DistanceToDataCenter] *= 0.9 
		@inbounds DataCenter[DistanceToNewDataCenter .< DistanceToDataCenter]=NewDataCenter[DistanceToNewDataCenter .<DistanceToDataCenter]
		@inbounds DistanceToDataCenter[DistanceToNewDataCenter .< DistanceToDataCenter]=DistanceToNewDataCenter[DistanceToNewDataCenter.<DistanceToDataCenter]
		if all(abs.(Gradient).<0.0001) break end
	end
	DataCenter
end

Wp