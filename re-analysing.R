## load ----
library(tidyverse)
library(iNEXT)
library(vegan)
library(doMC)
library(Hmisc)
registerDoMC(24)
counts <- read_tsv("tcr_counts.txt") %>%
  as.data.frame() %>% column_to_rownames("count") %>% as.matrix()
colnames(counts)[1] <- "35.base" 


## basic information ----
par(mfrow = c(1,2))
barplot(colSums(counts),horiz = T, las = 2, main = "number of molecule captured")
barplot(colSums(counts>0),horiz = T, las = 2, main = "number of clone detected")


## rarefaction ----
counts.inext <- iNEXT(counts, q=0, knots=200, se=F, conf=0.95, 
                      nboot=1,doMC = T)
plot(counts.inext,se=F,show.main = F)
plot(counts.inext,se=F,show.main = F, type = 3,show.legend = F)

## simulation ----

### generation of in silico population ----
population <- counts[,1:5] %>% 
  rowSums() %>%
  {.[.>0]} %>% {c((.*500)[.*500>max(.)], rep(.,100))}
base1 <- rrarefy(population, 10000)[,]
base2 <- rrarefy(population, 7000)[,]
base3 <- rrarefy(population, 5000)[,]
base4 <- rrarefy(population, 6000)[,]
mat1 <- cbind(base1,base2,base3,base4)
mat1 <- mat1[rowSums(mat1)>0,]
mat1.inext <- iNEXT(mat1, q=0, knots=100, 
                    se=F, nboot=1,doMC = T)
plot(mat1.inext,se=F,show.main = F,xlim = c(0,140000),ylim = c(0,10000))
plot(counts.inext,se=F,show.main = F,xlim = c(0,140000),ylim = c(0,10000))
rm(base1,base2,base3,base4,mat1,mat1.inext)


### estimate seed size/decay rate/amp rate ----
S <- rrarefy(population, 3*10^5)
NSsimul <- function(start = 3*10^5,
                    seed = 5000,
                    seed.amp.rate = 0.363,
                    decay.rate = 0.85,
                    cycle = 19){
  S <- rrarefy(population, start) %>% t %>% .[.>0,, drop = F]
  S[1] = seed
  for(i in 1:cycle){
    S_ <- S[,ncol(S),drop = F]
    S_[1,1] <- round(S_[1,1]*(seed.amp.rate+1))
    S__ <- t(rrarefy(S_,round(sum(S_)*decay.rate)))
    S <- cbind(S,S_,S__)
  }  
  S
}
S <- NSsimul()
colnames(S)[2:39] <- rep(1:19, each=2)
colnames(S)[1] <- "i"


par(mfcol=c(2,2),mar=c(3,4,0,0),oma=c(0,0,5,0))
barplot(colSums(S), ylab = "sample size")
text(43,130000,paste0("r=",round(S[1,ncol(S)]/sum(S[,ncol(S)]),4)))
barplot(S[1,], col="grey50",add=T)
legend("topright","expanded", fill = "grey50",bty="n")

abline(h = 100000, lty=2)
barplot(colSums(S>0), col = "grey50", ylab = "species richness")
plot(counts.inext,se=F,show.legend = F,show.main = F,xlim = c(0,150000),ylim=c(0,35000), col = "grey")
lines(colSums(S),colSums(S>0),col="blue",lty = 2)
points(colSums(S),colSums(S>0))
abline(v= c(8000,50000,100000),h=c(5000,7500), lty = 2, col = "grey")
plot(counts.inext,se=F,show.legend = F,show.main = F,xlim = c(0,150000),ylim=c(0,10000), col = "grey")
lines(colSums(S),colSums(S>0),type = 'l', col = "blue", lty =2)
points(colSums(S),colSums(S>0), xlim = c(0,150000),ylim=c(0,10000))
abline(v= c(8000,50000,100000),h=c(5000,7500), lty = 2, col = "grey")
mtext("Simulation-based estimation of amplification-decay rate",outer = T, cex = 1.2, font=2,line=2.5)
mtext("initial = 3*10^5, seed = 5000, seed.amp.rate = 0.363, decay.rate = 0.85, cycle = 19",outer = T, cex = 1, font=1,line=1)

## question
par(mfrow=c(2,2),mar = c(4,4,1,1),oma=c(0,0,5,0))
plot(rowSums(counts[,1:5]), counts[,6], xlim = c(0,1000), ylim = c(0,1000), col=  "#00000030")
plot(rowSums(counts[,1:5]), counts[,6], xlim = c(0,100), ylim = c(0,100), col=  "#00000030")
plot(rowSums(counts[,1:5]), counts[,7], xlim = c(0,1000), ylim = c(0,1000), col=  "#00000030")
plot(rowSums(counts[,1:5]), counts[,7], xlim = c(0,100), ylim = c(0,100), col=  "#00000030")
mtext("which clones are expanding?", cex=1.2,font=2,outer=T)


## seed size and amp rate distribution based on top clones ----
par(mfrow=c(2,1),mar = c(4,4,1,1),oma=c(0,0,0,0))
h2 <- counts[,6:34] %>% sweep(2,colSums(.),`/`) %>% hist(breaks = 30, ylim = c(0,100),col="grey",main="")
h <- h2
h$counts <- rev(cumsum(rev(h$counts)))
h$density <- rev(cumsum(rev(h$density)))
plot(h,ylim = c(0,250), col="grey80", main =  "", ylab = "cumulative frequency")
(0:1000/1000)%>% lines(.,1.8/., lty = 2, col = "black")
(0:1000/1000)%>% lines(.,0.805/.^1.5, lty = 2, col = "red")
(0:1000/1000)%>% lines(.,.3600/.^2, lty = 2, col = "darkgreen")
(0:1000/1000)%>% lines(.,.16100/.^2.5, lty = 2, col = "blue")
text(0.3,120,expression(x*y==1.8))
text(0.3,100,expression(x^1.5*y==0.805),col="red")
text(0.3,80,expression(x^2*y==0.36),col="darkgreen")
text(0.3,60,expression(x^2.5*y==0.161),col="blue")

u <- (3*10^5)/rarefy(population, 3*10^5)
seed.distribution.summary <- data.frame("%" = c(40,20,10,7,5,3,1,0.5,0.1), check.names=F) %>%
  mutate(`cumulative n.clone` = h$counts[match(`%`/100,h$breaks)],
         `n.clone in 3*10^5 * 29` = diff(c(0,`cumulative n.clone`)),
         `extrapolated cumsum` = 805/`%`^1.5,
         `extrapolated n.clone` = diff(c(0,`extrapolated cumsum`)),
         `proportion sum %` = `%` * `extrapolated n.clone` / 29,
         `proportion cumsum %` = cumsum(`proportion sum %`),
         ad = ((10^5*`%`/100)/u)^(1/19),
         a = `ad`/0.85,
         `n.clone in population` = `extrapolated n.clone`*length(population)/sum(drarefy(population,300000))/29) %>% round(2)

htmlTable::interactiveTable(seed.distribution.summary) 

## simulation
spawn_specificity <- function(){
  a = rep(1,length(population))
  ai = 1:length(a)
  cv = rep("black", length(population))
  ct = colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(8) %>% rev
  # ct = c("#5E4FA2", "#3F96B7", "#88CFA4", "#D7EF9B", "#FFFFBF", "#FDD380", "#F88D51", "#DC494C", "#9E0142")
  for(i in 1:8){
    size = round(seed.distribution.summary$`n.clone in population`[i])
    sampled <- sample(ai,size)
    a[sampled] = seed.distribution.summary$a[i]
    ai = ai[-sampled]
    cv[sampled] = ct[i]
    names(a) = cv
  }
  a
}
tmp <- spawn_specificity()
table(tmp)

simul <- function(S = rrarefy(population, 300000) %>% t(),
                  sp = spawn_specificity(),
                  cycle=19,
                  decay=0.85,
                  final=rnorm(1,50000,20000)){
  out = S
  for (i in 1:cycle){
    out = round(out*sp)
    out = rrarefy(t(out),round(sum(out)*decay))[,]
  }
  out = rrarefy(t(out),final)[,]
  return(list(count=out,color=names(sp)))
}

S.init <- rrarefy(population, 300000) %>% t()
sp = spawn_specificity()

combine_simul <- function(x,y){list(count=cbind(x$count,y$count),
                                    color = cbind(x$color,y$color))}

S1 <- foreach(i=1:20, .combine = combine_simul) %dopar% {simul()}
# S2 <- foreach(i=1:20, .combine = combine_simul) %dopar% {simul(S=S.init)} # fixed initiation
# S3 <- foreach(i=1:20, .combine = combine_simul) %dopar% {simul(sp=sp)}    # fixed specificity


colnames(S1$count) <- paste0("S1_",1:ncol(S1$count))

B <- foreach(s=c(13000,15000,17000,10000),.combine="cbind") %dopar% {
  t(rrarefy(population, s))
}
colnames(B) <- paste0("B",1:ncol(B))

S <- cbind(S1$count,B) %>% .[rowSums(.)>0,]
S.color <- S1$color[rowSums(cbind(S1$count,B))>0,]
S.color[S.color=="black"] = "#00000030"

S.inext <- iNEXT(S, q=0, knots=200, se=F, conf=0.95, nboot=1,doMC = T)

par(mfrow=c(2,3))
plot(S.inext, se=F, show.main = F)

# example simulation data
plot(S[,"B1"],S[,"S1_1"], col=S.color[,1],xlim=c(0,2800),ylim = c(0,2800), cex=ifelse(S.color[,1] == "#00000030", 1,2), pch=16);abline(0,1, lty = 2, col="grey")
plot(S[,"B1"],S[,"S1_1"], col=S.color[,1],xlim=c(0,150),ylim = c(0,150), cex=ifelse(S.color[,1] == "#00000030", 1,2), pch=16);abline(0,1, lty = 2, col="grey")

# example real count data
plot(rowSums(counts[,1:5]), counts[,"03.1516"], col = "#00000030", xlim = c(0,2800), ylim = c(0,2800));abline(0,1, lty = 2, col="grey")
plot(rowSums(counts[,1:5]), counts[,"03.1516"], col = "#00000030", xlim = c(0,150), ylim = c(0,150));abline(0,1, lty = 2, col="grey")


plot(S.inext$iNextEst[[21]]$m,S.inext$iNextEst[[21]]$qD, type = "l", xlim = c(0,60000), ylim = c(0,15000))
lines(S.inext$iNextEst[[1]]$m,S.inext$iNextEst[[1]]$qD)
abline(h=5000,v=0,lty=2,col="grey")
text(40000,10000,expression(F(x)== over(x*a,a+b)) )
text(c(3000,20000),5200,c("a","b"))


# standardization of counts

colnames(counts)
hist(colSums(counts>0),breaks = 20)
barplot(colSums(counts>0) %>% sort,las=2)
colSums(counts>0) %>% sort

counts.basesum <- rowSums(counts[,1:5])

A <- cbind(counts.basesum, counts)

colnames(A)[1] <- "basesum"
A.inext <- iNEXT(A, q=0, knots=200, se=F, conf=0.95, 
                      nboot=1,doMC = T)

plot(A.inext,se=F,show.main=F)
abline(h=2660)
barplot(colSums(A>0),las=2)
sort(colSums(A>0))[1:5]

A.sizeforq2660 <- lapply(A.inext$iNextEst, function(o){
  splinefun(x=o$qD,y=o$m)(2660)
}) %>% do.call(rbind,.)

A.sizeforq1372 <- lapply(A.inext$iNextEst, function(o){
  splinefun(x=o$qD,y=o$m)(1372)
}) %>% do.call(rbind,.)

barplot(t(A.sizeforq2660), las = 2)
barplot(t(A.sizeforq2660), las = 2,ylim=c(0,20000))
barplot(t(A.sizeforq1372), las = 2,ylim=c(0,20000))


A.sizefactor1 = A.sizeforq2660[1,]/A.sizeforq2660
A.sizefactor2 = A.sizeforq1372[1,]/A.sizeforq1372

barplot(t(A.sizefactor1), las = 2)
barplot(t(A.sizefactor2), las = 2)

A.std1 = A %>% sweep(2,A.sizefactor1, `/`)
A.std2 = A %>% sweep(2,A.sizefactor2, `/`)

A.std.sub1<- foreach(i = 7:35, .combine=cbind) %dopar% {
  tmp <- A.std1[,i]-A.std1[,1]
  tmp.cutoff <- tmp %>% {.[.<0]} %>% quantile(0.01) %>% {-.}
  pmax(tmp-tmp.cutoff,0)
}

A.std.sub2<- foreach(i = 7:35, .combine=cbind) %dopar% {
  tmp <- A.std2[,i]-A.std2[,1]
  tmp.cutoff <- tmp %>% {.[.<0]} %>% quantile(0.01) %>% {-.}
  pmax(tmp-tmp.cutoff,0)
}
par(mfrow=c(2,2))
barplot(colSums(A.std.sub1))
barplot(colSums(A.std.sub2))
barplot(colSums(A.std.sub1>0))
barplot(colSums(A.std.sub2>0))


# 
# plot(rowSums(counts.std[,1:5]), counts.std[,10], col = "#00000020", xlim = c(0,1000),ylim=c(0,1000))
# will take long time to plot!
par(mfrow = c(5,7), mar = c(1,1,1,1),oma=c(0,0,3,0),pty="s")
par(mfrow=c(1,2),pty="s",mar=c(5,4,1,1),oma=c(0,0,3,0))
# for(i in 7:35){
  i = 7
  
  tmp <- A.std[,i]-A.std[,1]
  hist(tmp,breaks=1000,xlim=c(-300,300),ylim=c(0,500),main="")
  tmp.cutoff <- tmp %>% {.[.<0]} %>% quantile(0.01) %>% {-.}
  abline(v=c(-tmp.cutoff,tmp.cutoff),col=c("blue","red"),lty=2)
  ifelse(tmp > tmp.cutoff,"red","#00000020") %>%
    plot(A.std[,1], A.std[,i], xlab="",ylab="", col = ., xlim = c(0,100),ylim=c(0,100),asp=1);
  abline(0,1,col="grey",lty=2);
  abline(-tmp.cutoff,1,col="blue",lty=2)
  abline( tmp.cutoff,1,col="red",lty=2)
  mtext(colnames(A.std)[i],outer=T,font=2,cex=1.2)
  table(tmp > tmp.cutoff)
  table(tmp < -tmp.cutoff)
  i=i+1
# }

A.std.sub<- foreach(i = 7:35, .combine=cbind) %dopar% {
  tmp <- A.std[,i]-A.std[,1]
  tmp.cutoff <- tmp %>% {.[.<0]} %>% quantile(0.01) %>% {-.}
  pmax(tmp-tmp.cutoff,0)
}
colnames(A.std.sub) <- colnames(A.std)[7:35]
A.std.sub

colSums(A.std.sub) %>% barplot(las=2,cex.names = 0.8)
colSums(A.std.sub>0) %>% barplot(las=2,cex.names = 0.8)

A.std.sub.cor <- cor(A.std.sub)
A.std.sub.cor2 <- cor(A.std.sub2)
A.std.cor <- cor(A.std)
A.cor <- cor(A)

adj.corrplot = function(x,...){
  x[x==1]=NA
  x[is.na(x)] = max(x,na.rm=T)
  corrplot::corrplot(x,...)
}



A.std.sub.cor %>% {.[.==1]=NA;.[is.na(.)]=max(.,na.rm=T);.} %>%
  corrplot::corrplot(corr = ., is.corr=F,method="square",order ="hclust")
A.std.sub.cor2 %>% {.[.==1]=NA;.[is.na(.)]=max(.,na.rm=T);.} %>%
  corrplot::corrplot(corr = ., is.corr=F,method="square",order ="hclust")

corrplot::corrplot(A.std.sub.cor,)
corrplot::corrplot(A.std.cor)
corrplot::corrplot(A.cor)

adj.corrplot(A.cor[6:35,6:35],is.corr = F,method = "square",)
adj.corrplot(A.std.cor[6:35,6:35],is.corr = F,method = "square")
adj.corrplot(A.std.sub.cor,is.corr = F,method = "square")


tmp.A <- read_rds("count.div-std.nq-sub.Rds")


tmp.A[1:10,6:10]
A.std.sub[1:10,1:5]

  