## load ----
library(tidyverse)
library(iNEXT)
library(vegan)
library(doMC)
library(Hmisc)
library(ComplexHeatmap)
registerDoMC(24)
counts <- read_tsv("tcr_counts.txt") %>%
  as.data.frame() %>% column_to_rownames("count") %>% as.matrix()
colnames(counts)[1] <- "35.base" 


## basic information ----
par(mfrow = c(1,2))
barplot(colSums(counts),horiz = T, las = 2, main ="number of molecule captured")
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
set.seed(42)
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
plot(counts.inext,se=F,show.legend = F,show.main = F,xlim = c(0,150000),
     ylim=c(0,35000), col = "grey")
lines(colSums(S),colSums(S>0),col="blue",lty = 2)
points(colSums(S),colSums(S>0))
abline(v= c(8000,50000,100000),h=c(5000,7500), lty = 2, col = "grey")
plot(counts.inext,se=F,show.legend = F,show.main = F,xlim = c(0,150000),
     ylim=c(0,10000), col = "grey")
lines(colSums(S),colSums(S>0),type = 'l', col = "blue", lty =2)
points(colSums(S),colSums(S>0), xlim = c(0,150000),ylim=c(0,10000))
abline(v= c(8000,50000,100000),h=c(5000,7500), lty = 2, col = "grey")
mtext("Simulation-based estimation of amplification-decay rate",outer = T, 
      cex = 1.2, font=2,line=2.5)
mtext(paste(
  "initial = 3*10^5, seed = 5000, seed.amp.rate = 0.363,",
  "decay.rate = 0.85, cycle = 19"),
      outer = T, cex = 1, font=1,line=1)

## questional  ----
par(mfrow=c(2,2),mar = c(4,4,1,1),oma=c(0,0,5,0), asp=1,pty="s")
plot(rowSums(counts[,1:5]), counts[,6], xlim = c(0,1000), ylim = c(0,1000), 
     col=  "#00000030")
plot(rowSums(counts[,1:5]), counts[,6], xlim = c(0,100), ylim = c(0,100), 
     col=  "#00000030")
plot(rowSums(counts[,1:5]), counts[,7], xlim = c(0,1000), ylim = c(0,1000), 
     col=  "#00000030")
plot(rowSums(counts[,1:5]), counts[,7], xlim = c(0,100), ylim = c(0,100), 
     col=  "#00000030")
mtext("which clones are expanding?", cex=1.2,font=2,outer=T)


## seed size and amp rate distribution based on top clones ----
par(mfrow=c(2,1),mar = c(4,4,1,1),oma=c(0,0,0,0))
h2 <- counts[,6:34] %>% sweep(2,colSums(.),`/`) %>% 
  hist(breaks = 30, ylim = c(0,100),col="grey",main="")
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
seed.distribution.summary <- data.frame("%" = c(40,20,10,7,5,3,1,0.5,0.1), 
                                        check.names=F) %>%
  mutate(`cumulative n.clone` = h$counts[match(`%`/100,h$breaks)],
         `n.clone in 3*10^5 * 29` = diff(c(0,`cumulative n.clone`)),
         `extrapolated cumsum` = 805/`%`^1.5,
         `extrapolated n.clone` = diff(c(0,`extrapolated cumsum`)),
         `proportion sum %` = `%` * `extrapolated n.clone` / 29,
         `proportion cumsum %` = cumsum(`proportion sum %`),
         ad = ((10^5*`%`/100)/u)^(1/19),
         a = `ad`/0.85,
         `n.clone in population` = 
`extrapolated n.clone`*length(population)/sum(drarefy(population,300000))/29)%>%
  round(2)

htmlTable::interactiveTable(seed.distribution.summary) 

## multi-clonal expansion simulation ----
spawn_specificity <- function(){
  a = rep(1,length(population))
  ai = 1:length(a)
  cv = rep("black", length(population))
  ct = colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(8) %>% rev
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
# S2 <- foreach(i=1:20, .combine = combine_simul) %dopar% {simul(S=S.init)} 
# fixed initiation
# S3 <- foreach(i=1:20, .combine = combine_simul) %dopar% {simul(sp=sp)}    
# fixed specificity


colnames(S1$count) <- paste0("S1_",1:ncol(S1$count))

B <- foreach(s=c(13000,15000,17000,10000),.combine="cbind") %dopar% {
  t(rrarefy(population, s))
}
colnames(B) <- paste0("B",1:ncol(B))

S <- cbind(S1$count,B) %>% .[rowSums(.)>0,]
S.color <- S1$color[rowSums(cbind(S1$count,B))>0,]
S.color[S.color=="black"] = "#00000030"

S.inext <- iNEXT(S, q=0, knots=200, se=F, conf=0.95, nboot=1,doMC = T)


### example simulation & real data ----
par(mfrow=c(2,3))
plot(S.inext, se=F, show.main = F)

plot(S[,"B1"],S[,"S1_1"], col=S.color[,1],xlim=c(0,2800),ylim = c(0,2800), 
     cex=ifelse(S.color[,1] == "#00000030", 1,2), pch=16)
abline(0,1, lty = 2, col="grey")
plot(S[,"B1"],S[,"S1_1"], col=S.color[,1],xlim=c(0,150),ylim = c(0,150), 
     cex=ifelse(S.color[,1] == "#00000030", 1,2), pch=16)
abline(0,1, lty = 2, col="grey")

plot(rowSums(counts[,1:5]), counts[,"03.1516"], col = "#00000030", 
     xlim = c(0,2800), ylim = c(0,2800));abline(0,1, lty = 2, col="grey")
plot(rowSums(counts[,1:5]), counts[,"03.1516"], col = "#00000030", 
     xlim = c(0,150), ylim = c(0,150));abline(0,1, lty = 2, col="grey")

plot(S.inext$iNextEst[[21]]$m,S.inext$iNextEst[[21]]$qD, type = "l", 
     xlim = c(0,60000), ylim = c(0,15000))
lines(S.inext$iNextEst[[1]]$m,S.inext$iNextEst[[1]]$qD)
abline(h=5000,v=0,lty=2,col="grey")
text(40000,10000,expression(X_adj== X*over(a+b,a)) )
text(c(3000,20000),5200,c("a","b"))


## standardization of counts ----

# calculation
A <- cbind(rowSums(counts[,1:5]), counts)
colnames(A)[1] <- "basesum"

A.inext <- iNEXT(A, q=0, knots=200, se=F, conf=0.95, 
                      nboot=1,doMC = T)

A.sizeforq2660 <- lapply(A.inext$iNextEst, function(o){
  splinefun(x=o$qD,y=o$m)(2660)
}) %>% do.call(rbind,.)

A.sizeforq1372 <- lapply(A.inext$iNextEst, function(o){
  splinefun(x=o$qD,y=o$m)(1372)
}) %>% do.call(rbind,.)


A.sizefactor1 = A.sizeforq2660[1,]/A.sizeforq2660
A.sizefactor2 = A.sizeforq1372[1,]/A.sizeforq1372

A.std1 = A %>% sweep(2,A.sizefactor1, `/`)
A.std2 = A %>% sweep(2,A.sizefactor2, `/`)


# plot

par(mfrow = c(2,3))
plot(S.inext$iNextEst[[21]]$m,S.inext$iNextEst[[21]]$qD, type = "l", 
     xlim = c(0,60000), ylim = c(0,15000))
lines(S.inext$iNextEst[[1]]$m,S.inext$iNextEst[[1]]$qD)
abline(h=5000,v=0,lty=2,col="grey")
text(40000,10000,expression(X_adj== X*over(a+b,a)) )
text(c(3000,20000),5200,c("a","b"))

barplot(colSums(A>0) %>% sort,las=2, main = "number of clones", cex.names = 0.7)
abline(h=sort(colSums(A>0))[c(1,3)],col=c("red","blue"))

colSums(A>0) %>% sort %>% .[1:5]

plot(A.inext, se = F, show.main = F, show.legend = F, col = "grey",
     ylim=c(0,20000))
abline(h=sort(colSums(counts>0))[c(1,3)],col=c("red","blue"))
text(120000,sort(colSums(A>0))[c(1,3)],sort(colSums(A>0))[c(1,3)],
     col=c("red","blue"), adj=c(1,-0.5))
sort(colSums(A>0))[1:5]

barplot(t(1/A.sizefactor1), las = 2, main = "size factor for q=2660")
barplot(t(1/A.sizefactor2), las = 2, main = "size factor for q=1372")

# subtraction of cutoff by negative quantile ----
A.std.sub1<- foreach(i = 7:35, .combine=cbind) %dopar% {
  tmp <- A.std1[,i]-A.std1[,1]
  tmp.cutoff <- tmp %>% {.[.<0]} %>% quantile(0.001) %>% {-.}
  pmax(tmp-tmp.cutoff,0)
}

A.std.sub2<- foreach(i = 7:35, .combine=cbind) %dopar% {
  tmp <- A.std2[,i]-A.std2[,1]
  tmp.cutoff <- tmp %>% {.[.<0]} %>% quantile(0.001) %>% {-.}
  pmax(tmp-tmp.cutoff,0)
}
colnames(A.std.sub1) <- colnames(A)[7:35]
colnames(A.std.sub2) <- colnames(A)[7:35]

A.std.cutoff1 <- foreach(i = 7:35, .combine=c) %dopar% {
  tmp <- A.std1[,i]-A.std1[,1]
  tmp.cutoff <- tmp %>% {.[.<0]} %>% quantile(0.001) %>% {-.}
}
A.std.cutoff2 <- foreach(i = 7:35, .combine=c) %dopar% {
  tmp <- A.std2[,i]-A.std2[,1]
  tmp.cutoff <- tmp %>% {.[.<0]} %>% quantile(0.001) %>% {-.}
}
names(A.std.cutoff1) <- colnames(A)[7:35] 
names(A.std.cutoff2) <- colnames(A)[7:35]

# example histogram and scatterplot ----
par(mfcol = c(2,4), mar = c(5,4,3,1),oma=c(0,0,0,0))
for(i in c(7,10,15,20)){
  tmp <- A.std[,i]-A.std[,1]
  tmp.cutoff <- tmp %>% {.[.<0]} %>% quantile(0.001) %>% {-.}
  hist(tmp,freq=F,breaks=c(-Inf,-60:60*5,Inf),xlim=c(-300,300),ylim=c(0,0.001),
       main=colnames(A.std)[i])
  abline(v=c(-tmp.cutoff,tmp.cutoff),col=c("blue","red"),lty=2)
  legend("topright",legend=c(sum(tmp > tmp.cutoff),sum(tmp < -tmp.cutoff)),
         pch=c("+","-"),bty="n")
  plot(A.std[,1], A.std[,i], xlab="",ylab="", 
       col = ifelse(tmp > tmp.cutoff,"red",
                    ifelse(tmp < -tmp.cutoff,"blue","#00000020")), 
       xlim = c(0,200),ylim=c(0,200),asp=1,
       main = colnames(A.std)[i]);
  abline(0,1,col="grey",lty=2);
  abline(-tmp.cutoff,1,col="blue",lty=2)
  abline( tmp.cutoff,1,col="red",lty=2)
}

# plot adj.count.sum, number of positive clone, and cutoffs ----
par(mfcol = c(2,3))
barplot(colSums(A.std.sub1), 
        main = "Sum of adjusted,nq-subtracted count, q=2660")
barplot(colSums(A.std.sub2), 
        main = "Sum of adjusted,nq-subtracted count, q=1372")
barplot(colSums(A.std.sub1>0), main = "clone count, q=2660")
barplot(colSums(A.std.sub2>0), main = "clone count, q=1372")
barplot(A.std.cutoff1, 
        main = "0.001 quantile of negative values of A_adj-B_adj with q=2660")
barplot(A.std.cutoff2, 
        main = "0.001 quantile of negative values of A_adj-B_adj with q=1372")


# save.image("~/Projects/cocult/cocult2/snapshot_20191008.RData")
# load("~/Projects/cocult/cocult2/snapshot_20191008.RData")


# structure of T cell response
A.std2[1:5,1:8]
A.std.sub2[1:5,1:8]

## heatmap basic

cor(A.std.sub2) %>%
  {.[.==1]<-NA;.[is.na(.)]<-max(.,na.rm=T);.} %>%
  Heatmap(col=RColorBrewer::brewer.pal(9,"Blues"), 
          width = unit(10,"cm"),height=unit(10,"cm"),
          name = "Correlation coef",
          top_annotation = HeatmapAnnotation(
            foo1 = 1:29,
            sum = anno_barplot(colSums(A.std.sub2)),
            nclone=anno_barplot(colSums(A.std.sub2>0))),
          right_annotation = rowAnnotation(foo2 = 29:1))

group_info <- tribble(~name,       ~type,   ~pearson_group,
                      "basesum",   "base",  1,
                      "35.base",   "base",  1,
                      "30.base",   "base",  1,
                      "31.base",   "base",  1,
                      "32.base",   "base",  1,
                      "33.base",   "base",  1,
                      "26.rest",   "rest",  0,
                      "27.rest",   "rest",  0,
                      "28.rest",   "rest",  4,
                      "29.rest",   "rest",  0,
                      "03.1516",   "1516",  3,
                      "07.1516",   "1516",  3,
                      "08.1520",   "1520",  0,
                      "01.1520_c", "1520",  0,
                      "04.1599",   "1599",  0,
                      "20.1599",   "1599",  0,
                      "02.1459",   "1459",  0,
                      "15.1459",   "1459",  0,
                      "05.1514",   "other", 2,
                      "06.1449",   "other", 0,
                      "09.1546",   "other", 5,
                      "10.1543",   "other", 0,
                      "11.1590",   "other", 5,
                      "12.1508",   "other", 0,
                      "13.1560",   "other", 4,
                      "14.1496",   "other", 0,
                      "16.1571",   "other", 2,
                      "17.1530",   "other", 0,
                      "18.1510",   "other", 0,
                      "21.1569",   "other", 2,
                      "22.1499",   "other", 2,
                      "24.1466",   "other", 0,
                      "23.1461",   "other", 3,
                      "19.1577",   "other", 3,
                      "25.1583",   "other", 3)
tmp.type.col <- structure(unlist(plyr::revalue(
  group_info$type,list("base"="grey","rest"="blue","1516"="red","1520"="orange",
                       "1599"="purple","1459"="green"), warn_missing = F)), 
  group_info$type)

cor(A.std2[,-1]) %>%
  {.[.==1]<-NA;.[is.na(.)]<-max(.,na.rm=T);.} %>%
  Heatmap(col=RColorBrewer::brewer.pal(9,"Blues"), 
          width = unit(12,"cm"),height=unit(12,"cm"),
          name = "Correlation coef",
          top_annotation = HeatmapAnnotation(
            replicates=group_info$type[-1],
            pearson_group=as.character(group_info$pearson_group[-1]),
            sum_adj = anno_barplot(colSums(A.std2[,-1])),
            nclone = anno_barplot(colSums(A.std2[,-1]>0)),
            col = list(replicates=c("base"="grey60","rest"="grey30",
                                    "1516"="red","1520"="orange",
                                    "1599"="purple","1459"="green",
                                    "other"="white"),
                       pearson_group=c("0"="white",
                                       "1"="blue",
                                       "2"="red",
                                       "3"="orange",
                                       "4"="purple",
                                       "5"="black"))
            )
          )

# high basemeans
tmp <- A.std2[rowSums(A.std.sub2>0)>-1 & rowMeans(A.std2[,3:6]) > 10,]
Heatmap(tmp[,-1], name = "diversity standardized count",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", tmp[,-1][i, j]),
                    x, y, gp = gpar(fontsize = 6,col="white"))},
        width = unit(27,"cm"),height=unit(20,"cm"),
        clustering_distance_columns = "euclidean",
        clustering_distance_rows = "spearman", 
        clustering_method_rows = "ward.D",
        clustering_method_columns="ward.D",
        col=circlize::colorRamp2(
          c(0:11*60),c("black",rev(RColorBrewer::brewer.pal(11,"Spectral")))),
        right_annotation = rowAnnotation(basesum=anno_barplot(tmp[,"basesum"])),
        top_annotation = HeatmapAnnotation(
          replicates=group_info$type[-1],
          pearson_group=as.character(group_info$pearson_group[-1]),
          sum_adj = anno_barplot(colSums(A.std2)[-1]),
          nclone = anno_barplot(colSums(A.std2>0)[-1]),
          col = list(replicates=c("base"="grey60","rest"="grey30",
                                  "1516"="red","1520"="orange", 
                                  "1599"="purple","1459"="green",
                                  "other"="white"),
                     pearson_group=c("0"="white", "1"="blue", "2"="red",
                                     "3"="orange", "4"="purple", "5"="black"))),
        row_title = "33 clones with base1~4 mean > 10")

# common responders
tmp <- A.std2[rowSums(A.std.sub2>0)>=15 ,]
Heatmap(tmp[,-1], name = "diversity standardized count",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", tmp[,-1][i, j]),
                    x, y, gp = gpar(fontsize = 6,col="white"))},
        width = unit(27,"cm"),height=unit(20,"cm"),
        clustering_distance_columns = "euclidean",
        clustering_distance_rows = "euclidean", 
        clustering_method_rows = "ward.D",
        clustering_method_columns="ward.D",
        col=circlize::colorRamp2(
          c(0:11*60),c("black",rev(RColorBrewer::brewer.pal(11,"Spectral")))),
        right_annotation = rowAnnotation(basesum=anno_barplot(tmp[,"basesum"])),
        top_annotation = HeatmapAnnotation(
          replicates=group_info$type[-1],
          pearson_group=as.character(group_info$pearson_group[-1]),
          sum_adj = anno_barplot(colSums(A.std2)[-1]),
          nclone = anno_barplot(colSums(A.std2>0)[-1]),
          col = list(replicates=c("base"="grey60","rest"="grey30",
                                  "1516"="red","1520"="orange", 
                                  "1599"="purple","1459"="green",
                                  "other"="white"),
                     pearson_group=c("0"="white", "1"="blue", "2"="red",
                                     "3"="orange", "4"="purple", "5"="black"))),
        row_title = "35 clones respond(>0) in >=15/29 samples")











## distances ----
pairwise <- function(A,fun){
  n = ncol(A)
  o = matrix(0,ncol=n,nrow=n)
  for (i in 1:n) {
    for (j in 1:i) {
      o[i,j] = fun(A[,i],A[,j])
      o[j,i] = o[i,j]
    }
  }
  colnames(o) = colnames(A)
  rownames(o) = colnames(A)
  as.dist(o)
}

jaccard <- function(A){
  
}
bray <- function(x,y){
  w=ifelse(x==0&y==0,0,1)
  sum(abs(x-y))/sum(abs(x+y)*w)
}

tmp <- A.std2[rowSums(A.std.sub2>0)>15,]
tmp.f <- function(A){pairwise(A,bray)}
tmp.d <- tmp.f(tmp)
as.matrix(tmp.d)[1:5,1:5]


x=runif(10)
y=runif(10)
bray(1:10,runif(10))
pearson <- function(A){
  
}
manhattan <- function(x,y){
  
}

