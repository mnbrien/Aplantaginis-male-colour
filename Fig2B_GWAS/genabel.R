#### GWAS ####
# Brien, Orteu et al. 2022

install.packages("https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz", repos=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz", repos=NULL, type="source")

## start with your VCF and convert to BED using PLINK. Then comvert BED to TPED.

library(GenABEL)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(qqman)


# Import raw data and phenotypes into gwaa object
wy <- load.gwaa.data(genofile = "yy_mapped.raw", phenofile = "yy_mapped_phenotypes.txt")
wy <- load.gwaa.data(genofile = "ww_mapped.raw", phenofile = "ww_mapped_phenotypes.txt")


# QC
# Round 1
QC1 <- check.marker(wy, perid.call = 0, callrate = 0.5, p.level = 0)
wyQC <- wy[QC1$idok, QC1$snpok]
# Second round of QC
QC2 <- check.marker(wyQC, perid.call = 0.95, callrate = 0.95, p.level = 0.01) 
# Filter for second round of QC
wyQC <- wyQC[QC2$idok, QC2$snpok]


## GWAS without correction for population structure ###
fo.QT <- qtscore(colourWY ~ 1, data = wyQC, trait = "binomial")
summary(fo.QT, top=30)
plot(fo.QT, col = c("black", "grey"), ystart = 1)


### Estimation of genomic kinship matrix

# Get SNP names
all_snp_names <- wyQC@gtdata@snpnames
# Exclude unassigned scaffolds
included_snps <- all_snp_names[grep("YY_tarseq_206_arrow:", all_snp_names, invert = F)]

# Get genomic kinship matrix using autosomes
gkins_f2 <- ibs(wyQC[,autosomal(wyQC)], weight = "freq")
mds <- cmdscale(as.dist(0.5 - gkins_f2))
fo.QTibs <- qtscore(colourWY ~ mds[, 1] + mds[, 2],wyQC, trait.type = "binomial")
results<-summary(fo.QTibs, top=162)
sum(fo.QTibs[, "P1df"] <= 0.0001) # how many snps over threshold
plot(fo.QT, col = c("black", "black"), ystart = 1, cex=0.4)
abline(h=bf5HQ, col = "red", lty = 2)

#### Estimate Bonferroni corrected thresholds

bf5HQ <- -log(0.05/nsnps(wyQC), base = 10)
bf1HQ <- -log(0.01/nsnps(wyQC), base = 10)


### plot single scaffold only
results<-summary(fo.QTibs, top=416280)
scaf206<-subset(results, Chromosome=="WW_tarseq_419_arrow")
scaf206<- scaf206 %>% tibble::rownames_to_column()
names(scaf206)<-c("SNP", "CHR", "BP","Strand", "A1", "A2", "N", "effB" ,"se_effB","chi2.1df" , "P1df", "effAB" ,"effBB" ,     "chi2.2df", "P" ,"Pc1df")
scaf206$CHR <- gsub("WW_tarseq_419_arrow", "206", scaf206$CHR)
scaf206$CHR <-as.numeric(scaf206$CHR)
manhattan(scaf206, logp=T, suggestiveline = F, genomewideline = F)
abline(a=7.62, b=0, col="red", lty=2)

### plot using ggplot
gwas<-read.csv("white.yellow.gwas_yymapped_only206.csv")

p3 <- gwas%>% 
    ggplot(aes(Position, -log10(P1df)))+
    geom_point(colour="sienna2", size=0.9)+
    theme_bw(14)+
    theme(panel.grid = element_blank(), panel.background = element_blank(),
           panel.border = element_blank())+
    ylab("-Log10(p value)") +
    xlab("Position (bp) along scaffold 206") +
    geom_hline(yintercept=7.62, linetype="dashed", colour="gray28")+
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE))+
    scale_y_continuous(expand=c(0,0), limits=c(0,8.5)) 
