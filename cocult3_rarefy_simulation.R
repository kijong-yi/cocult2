
# prep ------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(iNEXT)
library(vegan)
library(doMC)
library(Hmisc)
registerDoMC(24)

counts <- read_tsv("tcr_counts.txt") %>%
  as.data.frame() %>% column_to_rownames("count") %>% as.matrix()
colnames(counts)[1] <- "35.base" 

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


# simulation ------------------------------------------------------------------------------------------------------

# population
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
} else {
  population <- read_rds("simulation/population.Rds")
}



# subsample bases
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


# amplificatoin
# amplification strength based on cocult top-clone percentages
if(F){
  par(mfrow = c(3,1), mar = c(5,4,1,1))
  counts[,6:34] %>% sweep(2, colSums(.), `/`) %>% hist(breaks = 30, ylim = c(0,100), xlim = c(0,0.4), main = "",xlab = "counts[,6:34] %>% sweep(2, colSums(.), `/`) %>% hist(breaks = 30, ylim = c(0,100))", col = "grey80")
  counts[,6:34] %>% sweep(2, colSums(.), `/`) %>% unlist %>% Hmisc::cut2(cuts = c(0.02,0.04, 0.06, 0.08, 0.12, 0.28, 0.52)) %>% table
  abline(v=c(0.02,0.04, 0.06, 0.08, 0.12, 0.28, 0.52), col = "red")
  text(c(0.03,0.05,0.07,0.1,0.2,0.4),rep(100,6),c("3%","5%","7%","10%","20%","40%"))
  counts[,6:34] %>% sweep(2, colSums(.), `/`) %>% unlist %>% Hmisc::cut2(cuts = c(0.02,0.04, 0.06, 0.08, 0.12, 0.28, 0.52)) %>% table %>%
    .[2:7] %>% {text(c(0.03,0.05,0.07,0.1,0.2,0.4),rep(95,6),.);.} %>%
    rev() %>% cumsum %>% {text(rev(c(0.03,0.05,0.07,0.1,0.2,0.4)),rep(90,6),.)}
  plot(rev(c(3,5,7,10,20,40)), c(105,32,29,23,21,1) %>% rev %>% cumsum, xlim = c(0,45))
  lines(1:1500/30, 5000*3/(1:1500), col = "grey", lty = 2)
  plot(c(40,25,10,7,5,3,1,0.5,0.1),
       c(1,21,23,29,32,105,289,500,4000) %>% cumsum, xlim =c(0,45), ylim = c(0,6000), pch = c(rep(1,6),c(2,2,2)))
  lines(1:1500/30, 5000*3/(1:1500), col = "grey", lty = 2)
  
  # table
  drarefy(x = population, sample = 3*10^5) %>% mean
  drarefy(x = population, sample = 3*10^5) %>% mean %>% {c(1,21,23,29,32,105,289,500,4000)/29/.}
  data.frame("%" = c(40,20,10,7,5,3,1,0.5,0.1),check.names = F) %>%
    mutate(`count` = `%` * 400, 
           `mean fold` = `count` / (3*10^5/rarefy(population, 3*10^5)),
           `n.clone in 3*10^5` = c(1,21,23,29,32,105,289,500,4000),
           `n.clone in population` = drarefy(x = population, sample = 3*10^5) %>% 
             mean %>% {c(1,21,23,29,32,105,289,500,4000)/29/.} %>% round) %>% knitr::kable()
}


#

spawn_specificity <- function(pop = population){
  mat <- tribble(~fold, ~size, ~color,
    3364.437652,    1, "#5E4FA2",
    1682.218826,   21, "#3F96B7",
     841.109413,   23, "#88CFA4",
     588.776589,   29, "#D7EF9B",
     420.554706,   32, "#FFFFBF",
     252.332824,  106, "#FDD380",
      84.110941,  292, "#F88D51",
      42.055471,  505, "#DC494C",
       8.411094, 4038, "#9E0142")
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

# strength tuning
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

amp <- function (l,amp=1,sp=growth.specificity1) {
  l*sp^amp
}

rare <- function(l, r) {
  rrarefy(t(l),round(r*colSums(l))) %>% t
}

# 24 base, 3*10^5

set.seed(42)
S <- rrarefy(population,3*10^5) %>% t()
S <- cbind(S,S,S,S,S,S,S,S,S,S,S,S)
colnames(S) = letters[1:12]
S.bak = S

# amp 
S <- S.bak

colSums(S)


S.init = t(rrarefy(population,1*10^5))
S.init = S.init[S.init>0, ,drop = F]

recursive_amp_rare <- function(init = t(rrarefy(population,3*10^5)),
                               amp.r = 2.3,
                               rare.r = 0.7,
                               iter = 10) {
  S = init
  S[1,1] = 10
  for(i in 1:iter){
    a=S[,ncol(S), drop = F]
    a[1,1] = a[1,1]+round(a[1,1]*amp.r)
    r=t(rrarefy(t(a), round(sum(a)*rare.r)))
    S = cbind(S,a,r)
  }
  colnames(S) <- letters[1:ncol(S)]
  S
}
S <- recursive_amp_rare(S.init, amp.r = 2.3, rare.r = 0.7)
colSums(S) %>% barplot()

dim(S)

registerDoMC(24)
tmp <- iNEXT(S,
             q=0,
             se = T,
             nboot = 40,
             size = c(0,5000,10000,15000,20000,40000,50000,75000,100000,125000,150000),
             datatype="abundance", doMC = T)
par(mfrow = c(2,1), mar = c(1,4,1,0),oma=c(4,0,3,1))
plot(tmp, xlim = c(0,1.5*10^5), ylim = c(0,35000), se = T, show.main = F, xaxt='n', show.legend = F)
plot(tmp, xlim = c(0,1.5*10^5), ylim = c(0,10000), se = T, show.main = F)
mtext("Universe = 1*10^5, Seed = 10, Seed.amp.rate = 2.3,\n Universe.reduction.rate = 0.7, n.iter = 10", outer = TRUE, cex = 1.5, font = 2)

count.inext <- read_rds("iNEXT.rds")
plot(count.inext, xlim = c(0,1.5*10^5), ylim = c(0,35000), se = T, show.main = F, xaxt='n', show.legend = F)
plot(count.inext, xlim = c(0,1.5*10^5), ylim = c(0,10000), se = T, show.main = F)
par(mfrow = c(1,1))

# -----------------------------------
library(tidyverse)
library(iNEXT)
library(vegan)
library(doMC)
library(Hmisc)

load("~/Projects/cocult/cocult2/tmpenv.RData")

tmp <- iNEXT(S[rowSums(S[,1:16])>0,1:16],
             q=0,
             se = F,
             size = c(5000,10000,15000,20000,50000,100000,150000,200000,300000,400000),
             datatype="abundance", doMC = T)

par(mfrow = c(2,2))
plot(tmp, xlim = c(0,1.5*10^5), ylim = c(0,35000), type = 1)
plot(tmp, type = 2)
plot(tmp, xlim = c(0,1.5*10^5), ylim = c(0,10000))


plot(count.inext, xlim = c(0,1.5*10^5), ylim = c(0,35000))
plot(count.inext, xlim = c(0,1.5*10^5), ylim = c(0,10000))

rarefy(population, 3*10^5*29)


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



