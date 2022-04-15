#### Linkage mapping ####
# Brien, Orteu et al. 2022

library(qtl)
library(dplyr)
library(qtl2convert)
library(qtl2)
library(ggplot2)

# load map
c<-read.cross(format="csv", file="YY_map.csv", genotypes = c("AB", "BB"), alleles = c("A", "B"),  estimate.map = FALSE, convertXdata = T)
c<-read.cross(format="csv", file="WW_map.csv", genotypes = c("AB", "BB"), alleles = c("A", "B"),  estimate.map = FALSE, convertXdata = T)

# data checking
summary.map(c)
plot.map(c)
# plot genotype frequencies
g <- pull.geno(c)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,2), las=1) 
for(i in 1:2) 
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AB", "BB")[i], ylim=c(0,1))
par(mfrow=c(1,1))

########################################################

#### Genome scans ####

c<-jittermap(c)
# calculate genotype probabilities
c<- calc.genoprob(c, step=2.0, off.end=0.0, error.prob=1.0e-4, map.function="haldane", stepwidth="fixed")

# as covariates
cross<-as.numeric(pull.pheno(c, "Family"))

# genome scan
c_scan1<- scanone(c, pheno.col=4, model="binary", method="hk", addcovar = cross)
plot(c_scan1)

perm.c <-scanone(c, pheno.col=4, n.perm=5000,  method="hk")
summary(perm.c) 
add.threshold(c_scan1, perms=perm.c, alpha=0.05, lty=2)

#to plot in rqtl2
all_qtl1 <- scan_qtl_to_qtl2(c_scan1)
plot(all_qtl1$scan1, all_qtl1$map, bgcolor="white", altbgcolor="aliceblue", hlines="NA", col="dodgerblue4", col.axis="white")
add_threshold(all_qtl1$map, thresholdA = 3.02, lty=2)

#### significant markers ####
summary(c_scan1, perms=perm.c, lodcolumn=1, alpha=0.05, pvalues=TRUE) 

# intervals
bayesint(c_scan1, chr=" 9" , prob=0.95, lodcolumn = 1)
find.marker(c, chr=" 9", pos=28.0)

#effect sizes
sim <- sim.geno(c, n.draws=128, step=2, err=0.001)
qtl<- makeqtl(sim, chr=c(" 9"), pos=c(28.9)) 
out.fq <- fitqtl(sim, pheno.col=4, qtl=qtl, formula=y~Q1) 
summary(out.fq)



