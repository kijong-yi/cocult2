
# modules -----------------------------------------------------------------------

library(tidyverse)
library(ComplexHeatmap)
library(corrplot)
library(Rtsne)


# collect count data -------------------------------------------------------------

read_migec_count <- function(dir) {
  A = read_tsv(paste(dir,"cdr_A_filter","TCRA.filtered.cdrblast.txt", sep = "/"),
               col_types = "ddcccccdddddddd")
  B = read_tsv(paste(dir,"cdr_B_filter","TCRB.filtered.cdrblast.txt", sep = "/"),
               col_types = "ddcccccdddddddd")
  a <- A$Count; names(a) <- paste0("(A)",A$`CDR3 amino acid sequence`)
  b <- B$Count; names(b) <- paste0("(B)",B$`CDR3 amino acid sequence`)
  c(a,b)
}

col_bind_with_fill_zero <- function(l) {
  u <- unique(unlist(lapply(l, names)))
  ll <- lapply(l, function(x){o <- x[u]; o[is.na(o)] <- 0; o})
  lm <- do.call(cbind, ll)
  colnames(lm) <- names(l)
  rownames(lm) <- u
  lm
}

counts <- system("ls migec/* -d", intern = T) %>%
  {o = lapply(.,read_migec_count); names(o) = basename(.);o} %>%
  col_bind_with_fill_zero()

colnames(counts)[8]
colnames(counts)[8] <- "08.1520" # fix typo

# reordering
counts <- counts[,c("183rest-3hr1", "30.base", "31.base", "32.base", "33.base",
                    "26.rest", "27.rest", "28.rest", "29.rest",
                    "03.1516", "07.1516", "08.1520", "01.1520_c",
                    "04.1599", "20.1599", "02.1459", "15.1459",
                    "05.1514", "06.1449", "09.1546", "10.1543", "11.1590", "12.1508", 
                    "13.1560", "14.1496", "16.1571","17.1530", "18.1510",
                    "21.1569", "22.1499", "24.1466","23.1461", "19.1577", "25.1583")]

counts %>% as.data.frame %>% rownames_to_column("count") %>% write_tsv("tcr_counts.txt")

# percentage data

read_migec <- function(dir) {
  A = read_tsv(paste(dir,"cdr_A_filter","TCRA.filtered.cdrblast.txt", sep = "/"),
               col_types = "ddcccccdddddddd")
  B = read_tsv(paste(dir,"cdr_B_filter","TCRB.filtered.cdrblast.txt", sep = "/"),
               col_types = "ddcccccdddddddd")
  a <- A$Percentage; names(a) <- paste0("(A)",A$`CDR3 amino acid sequence`)
  b <- B$Percentage; names(b) <- paste0("(B)",B$`CDR3 amino acid sequence`)
  c(a,b)
}

percentage <- system("ls migec/* -d", intern = T) %>%
{o = lapply(.,read_migec); names(o) = basename(.);o} %>%
  col_bind_with_fill_zero()

colnames(percentage)[8]
colnames(percentage)[8] <- "08.1520" # fix typo
percentage <- percentage[,colnames(counts)]
percentage %>% as.data.frame %>% rownames_to_column("percentage") %>% write_tsv("tcr_percentage.txt")
# colSums of percentage is slightly less than 2, maybe because of small number rounding
rm(read_migec, read_migec_count, col_bind_with_fill_zero)


# additional statistics -------------------------------------------------------------

tcr <- data.frame(
  totalsums = counts %>% rowSums,
  basesums = counts[,1:5] %>% rowSums,
  basemeans = percentage[,1:5] %>% rowMeans(),
  recurrence = rowSums(ifelse(counts == 0, 0, 1)),
  unique = rowSums(ifelse(counts == 0, 0, 1)) > 0
)


samples <- data.frame(
  is_base = c(T,T,T,T,rep(F,30)),
  is_rest = c(rep(F,4),rep(T,4),rep(F,26)),
  replicate = c("base","base","base","base","base",
                "rest","rest","rest","rest",
                "1516", "1516", "1520", "1520", "1599", "1599", "1459", "1459",
                rep("other",17)),
  total_captured_UMI=colSums(counts),
  total_number_of_clones=colSums(ifelse(counts == 0, 0, 1)),
  percentage_of_base_sum_over_10_UMI = 0,                                               
  sum_of_dominant_clones_over_1_percent= 0,
  row.names = colnames(counts), check.names = F
)


# tcr, recurrence
par(mfrow=  c(1,2))
suppressWarnings(plot(tcr$basesums,tcr$recurrence))
suppressWarnings(plot(tcr$basesums,tcr$recurrence, xlim = c(0,100)))
mtext("xlim=0~100")
abline(v=10,col="red")
# TCR with basesum>10 apears > 4 of samples (564 clones among 66579)
table(tcr$basesums > 10)
sum(tcr$basesums[tcr$basesums > 10])/sum(tcr$basesums)
sum(tcr$totalsums[tcr$basesums > 10])/sum(tcr$totalsums)

# 
par(mar = c(5,5,1,1))

percentage %>%
  {as.data.frame(.[,c(12,13)])*100} %>%
  {plot(.,
     pch=c(21,24)[grepl("^\\(A",rownames(percentage))+1],
     asp=1,
     # xlim = c(0,58),ylim = c(0,58),
     xlim = c(0,5),ylim = c(0,5),
     xlab = "1520 rep 1 (%)",
     ylab = "1520 rep 2 (%)",
     col=ifelse(apply(.,1,min) == 0,"red","grey50"), 
     las = 1,
     bg=c("#DBDBDB60", "#D4D4D460", "#66666660", "#5E5E5E60", "#454545",
          "#363636", "#242424", "#141414", "#0D0D0D")[Hmisc::cut2(tcr$basemeans, 
                        c(0,0.00001,0.0001,0.001,0.002,0.004,0.008,0.016))])}
abline(v = 0, h = 0)
# text(3,56,"(0)", cex=.7);
# text(3,23,"(0)", cex=.7);
# text(19.7,0,"(0)",cex=.7);
# text(15.7,4,"(0.006%)",cex=.7,adj=c(.5))


# imaginary cartoon ----------------------------------------------------------------

par(mar = c(1,1,1,1), mfrow = c(1,2))
set.seed(42)
pool = cbind(rep(1:100,100), rep(1:100, each = 100))
cex1 = rep(0.5, 100*100)
cex1[sample(1:10000, 100, replace = F)] = 1
col1 = ifelse((pool[,1]-40)^2+(pool[,2]-50)^2 < 20^2, "black","grey")
cex1[sample(1:10000, 40, replace = F)] = 2
cex1[cex1 == 2 & col1 == "grey"] = 0.5
plot(pool, col = col1, cex = cex1, pch = 21, bg = "white", bty = "n", xaxt= "n", yaxt = "n")
plotrix::draw.circle(60,50,25,border = "black", lty = 2)

col2 = ifelse((pool[,1]-60)^2+(pool[,2]-50)^2 < 25^2, "black","grey")
cex2 = cex1
bg2 = rep("white",10000)
cex2[sample(1:10000, 500, replace = F)] = 1
cex2[sample(1:10000, 100, replace = F)] = 5
cex2[sample(1:10000, 20, replace = F)] = 10
cex2[sample(1:10000, 5000, replace = F)] = 0
cex2[cex2 > 4 & col2 == "grey"] = 0.5
col2[cex2 >0.5 & col2 == "black" & cex2 != cex1]="red"

plot(pool[order(cex2),], col = col2[order(cex2)], cex = cex2[order(cex2)], pch = 21, bg = "white",
     bty = "n",xaxt= "n", yaxt = "n")
plotrix::draw.circle(40,50,20,border = "black", lty = 2)
rm(bg2, cex1, cex2, col1, col2, pool)
dev.off()
# 
X = counts[,10][counts[,10] > 0 & tcr$basesums > 0]
X0 = tcr$basesums[counts[,10] > 0 & tcr$basesums > 0]
par(mfrow = c(4,1), mar = c(0,5,1,1), oma = c(4,0,0,0))
hist(X - X0, breaks = 1000,                  main = "Histogram of (03.1516 - basesums)", xaxt = "n")
hist(X - X0, breaks = 1000, ylim = c(0,100), main = "Histogram of (03.1516 - basesums)", xaxt = "n")
plot(X - X0, X0, xaxt = "n", col = ifelse(X0==1,"red","black"))
plot(X - X0, X, xaxt = "n", col = ifelse(X==1,"red","black"))

table(X-X0)[c("-5","-4","-3","-2","-1","0","1","2","3","4","5")]
X[which(X-X0 == 0)] %>% table
X[which(X-X0 == 1)] %>% table
X[which(X-X0 == 2)] %>% table
X[which(X-X0 == -1)] %>% table

mylog <- function(x,base=2){
  ifelse(x==0,0,log(x, base=base)+1)
}


library(vegan)
data(BCI)
S <- specnumber(BCI) # observed number of species
(raremax <- min(rowSums(BCI)))
Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)



