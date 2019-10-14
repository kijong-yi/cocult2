#Wp.R
library(zeallot)
Wp_R <- function(Data, p=1.5,RequiredK = 1, trace = F){
  c(N,M) %<-% dim(Data)
  U = matrix(1,nrow=N,ncol=N-RequiredK+1)
  U[,1] = 1:N
  K = N
  Nk = rep(1,N)
  Z = Data
  W = matrix(1/M,nrow=K,ncol=M)
  WL = list()
  InitialK=K
  AllDistances=matrix(Inf,nrow=InitialK,ncol=InitialK)
  for (k1 in 1:(InitialK-1)) {
    for (k2 in (k1+1):InitialK) {
      AllDistances[k1,k2] = 0.5 * sum( (abs(Z[k1,]-Z[k2,])^p) * (W[1,]^p) ) 
      AllDistances[k2,k1] = AllDistances[k1,k2]
    }
  }
  UIndex=1
  while (K > RequiredK){
    # 1st steep: look for clusters to merge
    c(k2Min,k1Min) %<-% which(AllDistances == min(AllDistances), arr.ind = T)[1,]
    # merge k1Min and k2Min
    UIndex = UIndex + 1
    c(U[,UIndex], Nk,Z,K) %<-% .Merge(k1Min,k2Min,U[,UIndex-1],Data,Nk,M,Z,p,K)
    W = .GetNewW(Data,W,U[,UIndex],Z,InitialK,M,p)
    if (trace){append(WL,W)}
    AllDistances=.Update_AllDistances(Z,W,p,Nk,AllDistances,k1Min,k2Min,InitialK)
  }
  list(U,Z,Nk,W)
}

.Update_AllDistances <- function(Z,W,p,Nk,AllDistances, k1Min, k2Min, InitialK){
  #remove k2 from Z
  AllDistances[k2Min,]=Inf
  AllDistances[,k2Min]=Inf
  #update the distances related to z1
  for (k in 1:InitialK){
    if (k==k1Min|Nk[k]==0|k==k2Min) {
      next
    }
    avg_W = ((W[k1Min,] + W[k,])/2)^p
    AllDistances[k1Min,k] = ((Nk[k1Min] * Nk[k])/(Nk[k1Min] + Nk[k])) * sum((abs(Z[k1Min,] - Z[k,])^p)*avg_W)
    AllDistances[k,k1Min] = AllDistances[k1Min,k]
  }
  AllDistances
}

.GetNewW <- function(Data,W, Ui, Z, K, M,p){
  D=matrix(0,K,M)
  for (l in 1:K) {
    for (j in 1:M) {
      D[l,j] = sum(abs(Data[Ui==l,j]-Z[l,j])^p)
    }
  }
  D = D + 0.0001
  #Calculate the actual Weight for each column
  if (p != 1){
    for (l in 1:K) {
      for (j in 1:M) {
        tmp=D[l,j]
        W[l,j]= 1/sum((tmp/D[l,])^(1/(p-1)))
      }
    }
  } else {
    for (l in 1:K){
      MinIndex = which.min(D[l,])
      W[l,1:M] = 0 #necessary to zero all others
      W[l,MinIndex] = 1
    }
  }
  W
}

.Merge <- function(k1Min, k2Min, Up, Data, Nk, M, Z, p, K) {
  Nk[k1Min] = Nk[k1Min] + Nk[k2Min]
  Up[Up==k2Min] = k1Min
  Nk[k2Min] = 0
  Z[k1Min,] = .New_cmt(Data[Up==k1Min,],p)
  Z[k2Min,] = Inf
  K = K - 1
  return(list(Up, Nk, Z, K))
}

.New_cmt <- function (Data,p) { # works only when p > 1
  #Calculates the Minkowski center at a given p.
  #Data MUST BE EntityxFeatures and standardised.
  c(N,M) %<-% dim(Data)
  if (p==1) {
    DataCenter=apply(Data,2,median)
    return(DataCenter)
  }
  else if (p==2){
    DataCenter=apply(Data,2,mean)
    return(DataCenter)
  }
  else if (N==1){
    DataCenter=Data
    return(DataCenter)
  }
  Gradient = rep(0.001,M)
  DataCenter = colSums(Data)/N
  DistanceToDataCenter=colSums(abs(sweep(Data, 2, DataCenter, `-`))^p)
  NewDataCenter=DataCenter+Gradient
  DistanceToNewDataCenter=colSums(abs(sweep(Data, 2, DataCenter, `-`))^p)
  Gradient[DistanceToDataCenter < DistanceToNewDataCenter] = -1 * Gradient[DistanceToDataCenter < DistanceToNewDataCenter]
  while (T){
    NewDataCenter = DataCenter + Gradient
    DistanceToNewDataCenter=colSums(abs(sweep(Data, 2, DataCenter, `-`))^p,dims=1) 
    Gradient[DistanceToNewDataCenter >= DistanceToDataCenter] = 0.9 * Gradient[DistanceToNewDataCenter >= DistanceToDataCenter]
    DataCenter[DistanceToNewDataCenter<DistanceToDataCenter]=NewDataCenter[DistanceToNewDataCenter<DistanceToDataCenter]
    DistanceToDataCenter[DistanceToNewDataCenter<DistanceToDataCenter]=DistanceToNewDataCenter[DistanceToNewDataCenter<DistanceToDataCenter]
    if(all(abs(Gradient)<0.0001)){break}
  }
  DataCenter
}
