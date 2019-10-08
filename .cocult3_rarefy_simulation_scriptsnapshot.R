
# prep ------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(iNEXT)
library(vegan)
library(doMC)
library(Hmisc)
registerDoMC(24)

counts <- read_tsv("tcr_counts.txt") %>%
  as.data.frame() %>% column_to_rownames("count") %>% as.matrix()
colnames(counts)[1] <- "35.base" # typo

# rarefaction analysis --------------------------------------------------------------------------------------------

if(F){
  source("~/src/plot.iNEXT.R")
  out <- iNEXT(counts,
               q=c(0,1,2),
               datatype="abundance",
               size=NULL, 
               endpoint=NULL, 
               knots=100, 
               se=TRUE, 
               conf=0.95, 
               nboot=50)
  par(mfrow = c(2,2), mar = c(5,4,1,1))
  plot.iNEXT(out, type = 1, se = F, show.legend = F, show.main = F, Order =1, ylab = "Shannon diversity")
  plot.iNEXT(out, type = 3, se = F, show.legend = F, show.main = F, Order =1, ylab = "Shannon diversity")
  plot.iNEXT(out, type = 2, se = F, show.legend = F, show.main = F, Order =1, ylab = "Shannon diversity")
  
  par(mfrow = c(2,2), mar = c(5,4,1,1))
  plot.iNEXT(out, type = 1, se = F, show.legend = F, show.main = F, Order =2, ylab = "Simpson diversity")
  plot.iNEXT(out, type = 3, se = F, show.legend = F, show.main = F, Order =2, ylab = "Simpson diversity")
  plot.iNEXT(out, type = 2, se = F, show.legend = F, show.main = F, Order =2, ylab = "Simpson diversity")
  
  
  ggiNEXT(out,facet.var = "order") + facet_wrap(~order, scale = "free") + scale_shape_manual(values=c(rep(19,29),rep(17,5)))
  ggiNEXT(out, type=2) + scale_shape_manual(values=c(rep(19,29),rep(17,5)))
  ggiNEXT(out, type=3) + scale_shape_manual(values=c(rep(19,29),rep(17,5)))
  # write_rds(out, "iNEXT.rds")
  out
}


# simulation --------------------------------------------------------------

# in silico population generation -----------------------------------------

if(F){
  tmp <- counts[,1:5] %>% rowSums
  tmp <- tmp[tmp>0]
  population <- c((tmp*500)[tmp*500>max(tmp)],rep(tmp,100))
  system("mkdir -p simulation")
  write_rds(population,"simulation/population.Rds")

  max(tmp)
  max(tmp/sum(tmp))
  sum(tmp==1)
  (population/sum(population)) %>% sort(decreasing = T) %>% head(5)
  (tmp/sum(tmp)) %>% sort(decreasing = T) %>% head(5)
  length(population)
  colSums(counts)[1:5]; mean(colSums(counts)[1:5]); sd(colSums(counts)[1:5])
  colSums(counts)[6:34];
  mean(colSums(counts)[6:34]); sd(colSums(counts)[6:34])
  # need distance between samples - based population size adjustment
} else {
  population <- read_rds("simulation/population.Rds")
}

# subsample bases ---------------------------------------------------------

if(F){
  base1 <- rrarefy(population, rnorm(n = 1,mean = 10000, sd=4000))[,]
  base2 <- rrarefy(population, rnorm(n = 1,mean = 10000, sd=4000))[,]
  base3 <- rrarefy(population, rnorm(n = 1,mean = 10000, sd=4000))[,]
  base4 <- rrarefy(population, rnorm(n = 1,mean = 10000, sd=4000))[,]
  mat1 <- cbind(base1,base2,base3,base4, base_all = base1+base2+base3+base4)
  mat1 <- mat1[rowSums(mat1)>0,]
  rarecurve(t(mat1), step = 10, sample = min(colSums(mat1)))
  rarecurve(t(counts[,1:5]), step = 1, sample = min(colSums(counts[,1:5])))
}
# proliferate/decay simulation --------------------------------------------

if (F){
  recursive_amp_rare <- function(init = t(rrarefy(population,3*10^5)),
                                 amp.r = 2.3,
                                 rare.r = 0.7,
                                 iter = 10,
                                 seed.n = 100) {
    S = cbind(init,init)
    S[1,2] = seed.n
    for(i in 1:iter){
      a=S[,ncol(S), drop = F]
      a[1,1] = a[1,1]+round(a[1,1]*amp.r)
      r=t(rrarefy(t(a), round(sum(a)*rare.r)))
      S = cbind(S,a,r)
    }
    if(ncol(S)>26){
      colnames(S) <- paste0("D",(1:ncol(S)/2))
    } else {
      colnames(S) <- letters[1:ncol(S)]
    }
    S
  }
  
  set.seed(42)
  S.init = t(rrarefy(population,1*10^5))
  S.init = S.init[S.init>0, ,drop = F]
  
  simul <- function(a = 2.3, r = 0.7, s = 10, i = S.init, iter.n = 10,dry = F) {
    cat("simulate rare ...")
    S <- recursive_amp_rare(i, a, r, iter.n, seed.n = s)
    par(mfcol = c(2,2), mar = c(0,0,1,1),oma=c(4,4,4,1))
    x = colSums(S)
    y = colSums(S>0)
    x %>% barplot()
    abline(h=100000)
    y %>% t %>% barplot
    abline(h=8000)
    plot(x,y,xlim=c(0,1.5*10^5), ylim=c(0,34000),)
    lines(x,y,lty=2,col="grey")
    abline(h=c(5000,7500),v=c(10000,50000,100000), col="grey", lty=2)
    plot(x,y,xlim=c(0,1.5*10^5), ylim=c(0,10000))
    lines(x,y,lty=2,col="grey")
    abline(h=c(5000,7500),v=c(10000,50000,100000), col="grey", lty=2)
    if(dry){return(1)}
    ans <- readline("Go? [y/n]: ") 
    if(ans != "y") {stop("user stop")}
    registerDoMC(24)
    cat(" done\nrun iNEXT ...")
    tmp <- iNEXT(S,
                 q=0,
                 se = F,
                 nboot = 1,
                 size = c(0, 5000, 10000, 15000, 20000, 40000, 50000, 75000, 100000, 125000, 150000),
                 datatype="abundance", doMC = T)
    cat(" done\n plotting ...")
    real <-   read_rds("iNEXT.rds")
    par(mfcol = c(2,2), mar = c(0,0,1,1),oma=c(4,4,4,1))
    plot.iNEXT(real, se=F, show.main=F, show.legend=F, xlim=c(0,1.5*10^5), ylim=c(0,34000), xaxt='n', order=0, las=1)
    abline(h=c(5000,7500),v=c(10000,50000,100000), col="grey", lty=2)
    plot.iNEXT(real, xlim=c(0,1.5*10^5), ylim=c(0,10000), se=F, order=0, show.main=F, las=1, show.legend = F)
    abline(h=c(5000,7500),v=c(10000,50000,100000), col="grey", lty=2)
    plot.iNEXT(tmp, se=F, show.main=F, show.legend=F, xlim=c(0,1.5*10^5), ylim=c(0,34000), xaxt='n', yaxt='n', order=0)
    abline(h=c(5000,7500),v=c(10000,50000,100000), col="grey", lty=2)
    plot.iNEXT(tmp, xlim=c(0,1.5*10^5), ylim=c(0,10000), se=F, order=0, show.main=F, yaxt='n')
    abline(h=c(5000,7500),v=c(10000,50000,100000), col="grey", lty=2)
    mtext(paste0("universe = 3*10^5, seed = ",s,", seed.amp.rate = ",a,",\n universe.decay.rate = ",r,", cycle = ",iter.n), outer=T, cex=1.5, font=2)
    cat(" done\n")
    return(tmp)
  }

  a2.3r0.7 <- simul(2.3,0.7)
  a2.4r0.7 <- simul(2.4,0.7)
  a1.7r0.7s100 <- simul(1.7,0.7,100)
  a1.1r0.7s1000 <- simul(1.1,0.7,1000)
  a2.5r0.7 <- simul(2.5,0.7)
  a1.001r0.7s10000 <- simul(1.001,0.7,10000)
  a1.002r0.7s2000 <- simul(1.002,0.7,2000)
  
  S.init2 = t(rrarefy(population,3*10^5))
  S.init2 = S.init2[S.init2>0, ,drop = F]
  
  simul(1.15,0.7,1800, S.init2, 9, T)
  simul(1.15,0.7,1800, S.init2, 9, F)
  
  simul(0.363,0.85,5000, S.init2, 19, T)
  
  S <- recursive_amp_rare(S.init2, 0.363, 0.85, 19, seed.n = 5000)
  par(mfcol = c(2,2), mar = c(0,4,1,1),oma=c(5,0,4,1))
  x = colSums(S)
  y = colSums(S>0)
  x %>% barplot(xaxt='n', ylab = "sample size")
  barplot(5000 * 1.363^c(rep(0:19,each=2)) * 0.85^c(0,0,0,rep(1:18,each=2),19),add = T, col = "grey40", yaxt='n')
  legend("topright",legend=c("expanded"), fill = c("grey40"), bty ='n')
  abline(h=100000)
  y %>% t %>% barplot(col="grey50",ylab="species richness")
  abline(h=8000)
  plot(x,y,xlim=c(0,1.5*10^5), ylim=c(0,34000),xaxt='n', ylab = "species richness")
  lines(x,y,lty=2,col="grey")
  abline(h=c(5000,7500),v=c(10000,50000,100000), col="grey", lty=2)
  plot(x,y,xlim=c(0,1.5*10^5), ylim=c(0,10000), xlab = "sample size", ylab="species richness")
  lines(x,y,lty=2,col="grey")
  abline(h=c(5000,7500),v=c(10000,50000,100000), col="grey", lty=2)
  
  
  simul(0.363,0.85,5000, S.init2, 19, F) # final fit
}

# top clone percentages distribution --------------------------------------

if(F){
  
  d <- counts[,6:34] %>% sweep(2, colSums(.), `/`) %>% unlist %>%
    Hmisc::cut2(cuts = c(40,20,10,7,5,3,1,0.5,0.1)/100) %>% as.numeric %>% table()
  names(d) <- c(0,rev(c(40,20,10,7,5,3,1,0.5,0.1)))
  table1 <- data.frame("observed %" = c(40,20,10,7,5,3,1,0.5,0.1),check.names = F) %>%
    mutate(`n.clone in 3*10^5 * 29` = c(rev(d[5:10]),NA,NA,NA),
           `cumulative n.clone` = cumsum(`n.clone in 3*10^5 * 29`),
           `extrapolated.cumsum` = 805/`observed %`^1.5,
           `extrapolated n.clone` = diff(c(0,805 / `observed %`^1.5)),
           `proportion sum %` = `observed %` * `extrapolated n.clone` / 29,
           `proportion cumsum %` = cumsum(`proportion sum %`),
           `count in 10^5` = `observed %`/100 * 10^5, 
           `a %` = 100*(((10^5*`observed %`/100/4.755642)^(1/19))/0.85 - 1),
           `n.clone in population(23k)` = `extrapolated n.clone`*length(population)/sum(drarefy(population,3*10^5))/29) %>% round(1)
  table1 %>% htmlTable::htmlTableWidget()
  
  
  
  
  par(mfrow = c(3,1), mar = c(5,4,1,1),oma=c(0,0,5,0))
  counts[,6:34] %>% sweep(2, colSums(.), `/`) %>% hist(breaks = 30, ylim = c(0,100), xlim = c(0,0.4), main = "",xlab = "counts[,6:34] %>% sweep(2, colSums(.), `/`) %>% hist(breaks = 30, ylim = c(0,100))", col = "grey80")
  abline(v=c(0.02,0.04, 0.06, 0.08, 0.12, 0.28, 0.52), col = "red")
  text(c(0.03,0.05,0.07,0.1,0.2,0.4),rep(100,6),c("3%","5%","7%","10%","20%","40%"))
  counts[,6:34] %>% sweep(2, colSums(.), `/`) %>% unlist %>% Hmisc::cut2(cuts = c(0.02,0.04, 0.06, 0.08, 0.12, 0.28, 0.52)) %>% table %>%
    .[2:7] %>% {text(c(0.03,0.05,0.07,0.1,0.2,0.4),rep(95,6),.);.} %>%
    rev() %>% cumsum %>% {text(rev(c(0.03,0.05,0.07,0.1,0.2,0.4)),rep(90,6),.)}
  
  plot(table1$`observed %`[1:6],table1$`cumulative n.clone`[1:6], xlim = c(0,40), ylim = c(0,250))
  # points(table1$`observed %`[7],table1$`extrapolated n.clone`[7],pch=2)
  (1:1500/30) %>% lines(9*20    /.    , col = "black", lty = 2)
  (1:1500/30) %>% lines(9*20^1.5/.^1.5, col = "red", lty = 2)
  (1:1500/30) %>% lines(9*20^2  /.^2  , col = "darkgreen", lty = 2)
  (1:1500/30) %>% lines(9*20^2.5/.^2.5, col = "blue", lty = 2)
  
  text(30,110,expression(x * y == 180), col="black",cex=1.3,adj=c(0,0))
  text(30,90,expression(x^1.5 * y == 805), col="red",cex=1.3,adj=c(0,0))
  text(30,70,expression(x^2 * y == 3600), col="darkgreen",cex=1.3,adj=c(0,0),font=2)
  text(30,50,expression(x^2.5 * y == 16100), col="blue",cex=1.3, adj = c(0,0))
  
  
  
  h<-counts[,6:34] %>% sweep(2, colSums(.), `/`) %>% {.*100} %>% hist(breaks = 30, plot=F)
  h$counts <- rev(cumsum(rev(h$counts)))
  h$density <- rev(cumsum(rev(h$density)))
  plot(h, freq=TRUE, ylim = c(0,250), xlim = c(0,40), main = "", col="grey80",border="black")
  (1:1500/30) %>% lines(9*20    /.    , col = "black", lty = 2)
  (1:1500/30) %>% lines(9*20^1.5/.^1.5, col = "red", lty = 2)
  (1:1500/30) %>% lines(9*20^2  /.^2  , col = "darkgreen", lty = 2)
  (1:1500/30) %>% lines(9*20^2.5/.^2.5, col = "blue", lty = 2)
  
  text(30,110,expression(x * y == 180), col="black",cex=1.3,adj=c(0,0))
  text(30,90,expression(x^1.5 * y == 805), col="red",cex=1.3,adj=c(0,0))
  text(30,70,expression(x^2 * y == 3600), col="darkgreen",cex=1.3,adj=c(0,0),font=2)
  text(30,50,expression(x^2.5 * y == 16100), col="blue",cex=1.3, adj = c(0,0))

  mtext("Top clone percentage distribution in samples",outer=T, cex=1.3,font=2,line=1)
  
}


# seeding on the population & sampling & cocult ---------------------------

# from table1
spawn_specificity <- function(pop = population){
  mat <- tribble(~fold, ~size, ~color,
                 1.893,     3, "#9E0142",
                 1.825,     6, "#E25249",
                 1.760,    17, "#FBA45C",
                 1.727,    18, "#FEE899",
                 1.697,    29, "#EDF7A3",
                 1.652,    84, "#A1D9A4",
                 1.559,   656, "#48A0B2",
                 1.503,  1486, "#5E4FA2")
    fold = rep(1, length(population))
    color = rep("black", length(population))
    I = 1:length(pop)
    for (i in nrow(mat):1){
      idx = sample(I, size = mat$size[i])
      fold[idx] <- mat$fold[i]
      color[idx] <- mat$color[i]
      I = I[-idx]
    }
    attr(fold, "color") = color
    fold
}

growth.specificity1 <- spawn_specificity()

inspect <- function(l){
  n <- l %>% lapply(sum) %>% unlist
  S.obs <- l %>% lapply(function(x){sum(x>0)}) %>% unlist
  qD <- l %>% lapply(function(x){rarefy(x,17000, se = F, MARGIN = 1)}) %>% unlist
  SC <- l %>% lapply(function(x){1-rareslope(x,17000)}) %>% unlist
  hist(n); mtext(mean(n))
  hist(S.obs); mtext(mean(S.obs))
  hist(qD); mtext(mean(qD))
  hist(SC); mtext(mean(SC))
}

cocult <- function(P, r=0.85, sp = growth.specificity1, iter=19){
  for(i in 1:19){
    P <- P*sp
    P <- rrarefy(t(P),round(r*colSums(P))) %>% t
  }
  P
}


rnorm(100,mean = 3*10^5, sd = 10000) %>% hist

P <- foreach(size=rnorm(100,mean = 3*10^5, sd = 10000),.combine=rbind) %dopar% {
  t(rrarefy(population, size))}

P.result <- cocult(P)



par(mfcol = c(4,4))
P              %>% inspect
P %>% amp(0.5) %>% inspect
P %>% amp(1)   %>% inspect
P %>% amp(1.5) %>% inspect

par(mfcol = c(4,3))
P %>% amp(0.5) %>% rare(0.1) %>% inspect
P %>% amp(1)   %>% rare(0.1) %>% inspect
P %>% amp(1.5) %>% rare(0.1) %>% inspect

P %>% amp(0.5) %>% rare(0.1) %>% amp(0.5) %>% inspect
P %>% amp(0.5) %>% rare(0.1) %>% amp(1)   %>% inspect
P %>% amp(0.5) %>% rare(0.1) %>% amp(1.5) %>% inspect

P %>% amp(1)   %>% rare(0.1) %>% amp(0.5) %>% inspect
P %>% amp(1)   %>% rare(0.1) %>% amp(1)   %>% inspect
P %>% amp(1)   %>% rare(0.1) %>% amp(1.5) %>% inspect

P %>% amp(1.5) %>% rare(0.1) %>% amp(0.5) %>% inspect
P %>% amp(1.5) %>% rare(0.1) %>% amp(1)   %>% inspect
P %>% amp(1.5) %>% rare(0.1) %>% amp(1.5) %>% inspect



# one-time amplify model

stimulate <- function (pop = population,
                       start.size = 3*10^5,
                       end.size = round(rnorm(n = 1,mean = 40000, sd=10000)),
                       growth = spawn_specificity()
){
  s.pop <- vegan::rrarefy(pop, start.size)
  s.fold <- rep(1,length(s.pop))
  s.color <- rep("black",length(s.pop))
  
  for(i in growth){
    if(is.null(i$idx)) {
      idx <- sample(which(s.pop>0), i$nclones)
    } else {
      idx <- i$idx
    }
    s.fold[idx] <- i$fold
    s.color[idx] <- i$color
    s.pop[idx] <- round(s.pop[idx] * i$fold)
  }
  s.pop <- vegan::rrarefy(s.pop, end.size)[,]
  return(list(s.pop = s.pop, s.fold = s.fold, s.color = s.color))
}

# RARAR model



stimulate2 <- function (pop = population,
                        start.size = 3*10^5,
                        middle.size = 3*10^4,
                        end.size.prop = 0.5,
                        growth = spawn_specificity()){
  s.pop <- vegan::rrarefy(pop, start.size)
  # s.fold <- rep(1,length(population))
  # s.color <- rep("black",length(population))
  for(i in growth){
    if(is.null(i$idx)) {
      idx <- sample(which(s.pop>0), i$nclones)
    } else {
      idx <- i$idx
    }
    # s.fold[idx] <- i$fold
    # s.color[idx] <- i$color
    s.pop[idx] <- round(s.pop[idx] * i$fold ^ .5)
  }
  s.pop <- vegan::rrarefy(s.pop, middle.size + (sum(s.pop) - start.size))[,]  
  for(i in growth){
    if(is.null(i$idx)) {
      idx <- sample(which(s.pop>0), i$nclones)
    } else {
      idx <- i$idx
    }
    s.pop[idx] <- round(s.pop[idx] * i$fold ^ .5)
  }
  # s.pop <- vegan::rrarefy(s.pop, sum(s.pop) * end.size.prop)[,]
  return(list(s.pop = s.pop
              # , s.fold = s.fold, s.color = s.color
              ))
}

# simulation
if (F){
  registerDoMC(5)
  sim <- foreach(x = 1:50) %dopar% {
    stimulate() %>% .$s.pop %>% .[.>0] %>% 
    {
      list(n = sum(.), 
           S.obs = length(.),
           qD = rarefy(.,17000, se = FALSE, MARGIN = 1),
           SC = 1-rareslope(.,17000))
    }
  }
  sim.result <- matrix(unlist(sim), ncol = 4, byrow = TRUE)
  colnames(sim.result) <- c("n", "S.obs", "qD", "SC")
  sim.result <- as.data.frame(sim.result)
  
  sim2 <- foreach(x = 1:20) %dopar% {
    stimulate2() %>% .$s.pop %>% .[.>0] %>% 
    {
      list(n = sum(.), 
           S.obs = length(.),
           qD = rarefy(.,17000, se = FALSE, MARGIN = 1),
           SC = 1-rareslope(.,17000))
    }
  }
  sim2.result <- matrix(unlist(sim2), ncol = 4, byrow = TRUE)
  colnames(sim2.result) <- c("n", "S.obs", "qD", "SC")
  sim2.result <- as.data.frame(sim2.result)
  

  
  par(mfrow = c(4,3))
  
  hist(colSums(counts[,6:34]))
  hist(sim.result$n)
  hist(sim2.result$n)
  
  hist(colSums(counts[,6:34] > 0))
  hist(sim.result$S.obs)
  hist(sim2.result$S.obs)
  
  rarefy(counts[,6:34], 17000, se = FALSE, MARGIN = 2) %>% hist
  hist(sim.result$qD)
  hist(sim2.result$qD)
  
  counts[,6:34] %>% t %>% rareslope(17000) %>% {1-.} %>% hist(main = "Histogram of sample coverage of counts[,6:34]")
  hist(sim.result$SC)
  hist(sim2.result$SC)
}






if(F){  
  counts[,6:34] %>% rarefy(17000, se = FALSE, MARGIN = 2) %>% summary
  counts[,6:34] %>% t %>% rareslope(17000) %>% summary
  counts[,6:34] %>% rarefy(17000, se = FALSE, MARGIN = 2) %>% hist
  counts[,6:34] %>% t %>% rareslope(17000) %>% hist
  
  
  growth.specificity1 <- spawn_specificity()
  
  sim1 <- stimulate(population, growth = growth.specificity1)
  sim2 <- stimulate(population, growth = growth.specificity1)
  sim3 <- stimulate(population, growth = growth.specificity1)
  sim4 <- stimulate(population, growth = spawn_specificity())
  sim5 <- stimulate(population, growth = spawn_specificity())
  sim6 <- stimulate(population, growth = spawn_specificity())
  
  
  mat2 <- cbind(sim1=sim1$s.pop,
                sim2=sim2$s.pop,
                sim3=sim3$s.pop,
                sim4=sim4$s.pop,
                sim5=sim5$s.pop,
                sim6=sim6$s.pop,
                base1,
                base2,
                base3,
                base4)
  mat2 <- mat2[rowSums(mat2)>0,]
  dim(mat2)
  dim(counts)
  counts[,1:10] %>% rowSums() %>% {. > 0} %>% table # need adjust of parameters
  
  rarecurve(x = t(mat2), sample = 17000, step = 1)
  
  
  par(mfrow = c(1,2))
  rarecurve(t(mat2), step = 5, sample = min(colSums(mat1)))
  rarecurve(t(counts[,1:11]), step = 5, sample = min(colSums(counts[,1:11])))
  
  
  # plot(sim1$s.pop, base1, col = sim1$s.color); abline(v = 0, h = 0)
  # plot(as.data.frame(mat1))
  
  inext1 <- iNEXT(mat1,
                  q=c(0,1),
                  datatype="abundance",
                  size=NULL, 
                  endpoint=NULL, 
                  knots=30, 
                  se=F
                  # conf=0.95, 
                  # nboot=50
  )
  plot(inext1) 
}



