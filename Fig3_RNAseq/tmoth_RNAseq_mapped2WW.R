##################
#Prepare workspace 
##################

 library(dplyr)
 library(tidyverse)
 #library(rtracklayer)
 library(ggplot2)
 library(limma)
 library(edgeR)
 library(rstudioapi)
 library(ggrepel)






#set wd to where the script is saved
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##################
#Load all the data
##################

#Read in the count data 
counts <- read.csv("featureCounts_mapped2WW", sep = "\t", skip = 1, stringsAsFactors = F)
names(counts) <- gsub(".*(CAM.*)\\_Aligned.sortedByCoord.out.bam","\\1", names(counts))

#information about the samples
info <- read.csv("../data_info/seq_samples_RNAseq.csv", header = T)

info <- info %>% 
  dplyr::select(CAMID, Genotype, Stage, Family, Pupation.day, Dissection.day)

info$Genotype <- gsub("yy", "YY", info$Genotype)
info$Stage <-  gsub("5d", "d5", info$Stage)
info$Stage <-  gsub("-", "", info$Stage)



genotype <- as.factor(info$Genotype)
stage <- as.factor(info$Stage)
family <- as.factor(info$Family)
Pupation.day <- as.factor(info$Pupation.day)
Dissection.day <- as.factor(info$Dissection.day)

genStage <- as.factor(paste(info$Genotype, info$Stage, sep=""))
genStage <- as.factor(gsub("[-_]", "",genStage))

info$genstate <- as.factor(paste(info$Genotype, info$Stage, sep=""))
info$genstate <- factor(x = info$genstate, levels = c("WW72h","YY72h", "WWd5","YYd5", "WWMel","YYMel", "WWPremel", "YYPremel"))


########################
###Create the dge object and normalise
######################## 


# design matrix 
design <- model.matrix(~0+genStage)
colnames(design) <- levels(genStage)

# 1. create DGE object 
info$CAMID <- gsub("CAM0", "CAM", info$CAMID) #make camids match
info <- info[info$CAMID %in% colnames(counts)[7:length(counts)],] #select info only for sequenced samples
info <- info[order(info$CAMID),] #put info in same order as samples in count matrix
counts_matrix <- counts[,7:length(counts)]
d0 <- DGEList(counts=counts_matrix, genes=counts$Geneid)
#filter
keep.exprs <- filterByExpr(d0, design)
dge <- d0[keep.exprs,,keep.lib.sizes=FALSE]


# 3.Calculate normalization factors
d <- calcNormFactors(dge, method="TMM")
d

# 4. Voom
dev.off()
par(mfrow=c(2,1))
y <- voom(d, design, plot = T) #plot=T if want to see voom plot:mean-variance trend

corfit <- duplicateCorrelation(y,design,block=family)

corfit$consensus

fit <- lmFit(y,design,block=family,correlation=corfit$consensus)

cont.matrix <- makeContrasts(h72=WW72h-YY72h, d5=WWd5-YYd5,
                             Premel=WWPremel-YYPremel,
                             Mel=WWMel-YYMel, 
                             levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2, main="Final model: Mean-variance trend")

tmp <- eBayes(fit2)
idx=tmp[tmp$genes$genes=="jg1310",]


results.limma <- decideTests(tmp, p.value = 0.05)
summary(results.limma)

topTable(tmp, coef = 3, sort.by = "P", p.value = 0.05, number = Inf)
dev.off()
vennDiagram(results.limma)

################
# plot MDS 
################
dev.off()
a <- plotMDS(y, col=as.numeric(genStage), labels =genStage)
MDS$samples <- names(a$y)
MDS$logFCdim1 <-  a$x
MDS$logFCdim2 <-  a$y

info

MDS <- tibble(CAMID=info$CAMID,logFCdim1=a$x, logFCdim2=a$y )
MDS %>% 
  full_join(info) %>% 
  ggplot(aes(logFCdim1, logFCdim2, colour=Genotype))+
  geom_text(aes(label=Stage, fontface = "bold"), size=10)+
  theme_bw(base_size = 24)+
  theme(panel.grid.major = element_blank())+
  theme(strip.background = element_rect(fill = "#f0efef"))+
  theme(panel.grid.minor = element_blank())+
  ylab("Leading LogFC dim 2")+
  xlab("Leading LogFC dim 1")+
  scale_colour_manual(values = cbPalette)+
  theme(legend.position = "none")


################
# Volcano plots
################

# define a function to plot volcano plots
my_volcaqnoPlot <- function(x){
  ggplot(x,aes(logFC, -log10(adj.P.Val)))+
    geom_point(alpha=0.5, size=3)+
    #geom_point(data = filter(slopes, fmr.adj.P.Val<0.05, FMR<0), aes(colour="Signif_fmr"))+
    geom_point(data = filter(x, adj.P.Val<0.05), aes(colour="Signif"), size=3)+
    theme_bw(base_size = 24)+
    theme(panel.grid.major = element_blank())+
    theme(strip.background = element_rect(fill = "#f0efef"))+
    theme(panel.grid.minor = element_blank())+
    theme(legend.position = "none")+
    theme(axis.text = element_text(size=22))+
    theme(axis.title = element_text(size=24))+
    scale_color_manual(values = c("magenta"))
}


# pre-mel volcano plot
limmapremel <- topTable(tmp, coef = 3, sort.by = "P", n=Inf)
limmapremel
g <- my_volcaqnoPlot(limmapremel)+
  #geom_text_repel(data = filter(limmapremel, genes=="jg1310"), 
                  #aes(label = "yellow e"), size=6)+
  geom_text_repel(data = filter(limmapremel, genes=="jg1308"), 
                  aes(label = "valkea", fontface = "italic"), size=12)
g


# 72h 
limma72 <- topTable(tmp, coef = 1, sort.by = "P", n=Inf)
my_volcaqnoPlot(limma72)+geom_text_repel(data = filter(limma72, adj.P.Val<0.05), 
                                         aes(label = genes))

# 5 days 
limma5 <- topTable(tmp, coef = 2, sort.by = "P", n=Inf)
my_volcaqnoPlot(limma5)+geom_text_repel(data = filter(limma5, adj.P.Val<0.05), 
                                        aes(label = genes))

# Melanised 
limmamel <- topTable(tmp, coef = 4, sort.by = "P", n=Inf)
my_volcaqnoPlot(limmamel)+geom_text_repel(data = filter(limmamel, adj.P.Val<0.05), 
                                          aes(label = genes))


################
# List genes DE per stage
################

filter_stage <- function(x, y,z){
  x %>% 
    filter(adj.P.Val<z) %>% 
    mutate(Stage=y)
}

signif001 <- filter_stage(limma72, "72h",0.01) %>% 
  full_join(filter_stage(limma5, "5days",0.01)) %>% 
  full_join(filter_stage(limmapremel, "premel",0.01)) %>% 
  full_join(filter_stage(limmamel, "mel",0.01))

write.table(signif001,"mapww_significant_genes_allStages_001.csv", quote = FALSE, row.names = FALSE, sep = "\t")

signif005 <- filter_stage(limma72, "72h",0.05) %>% 
  full_join(filter_stage(limma5, "5days",0.05)) %>% 
  full_join(filter_stage(limmapremel, "premel",0.05)) %>% 
  full_join(filter_stage(limmamel, "mel",0.05))

write.table(signif005,"mapww_significant_genes_allStages_005.csv", quote = FALSE, row.names = FALSE, sep = "\t")


################
# Plot read counts of genes of interest
################

# build a data frame with gene log2 counts with info
lcpm <- as.data.frame(cpm(d0, log=TRUE)) %>% 
  mutate(gene=d0$genes$genes) %>% 
  gather(key="CAMID", value="log2cpm", -gene) %>% 
  full_join(info)

# plot 
my_plotCounts <- function(x,y){x %>% 
    filter(gene==y) %>% 
    ggplot(aes(Stage, log2cpm, fill=Genotype, shape=gene))+
    geom_boxplot(width=0.5, position = position_dodge(width = 0.75))+
    geom_point(size=2, position = position_jitterdodge(jitter.width = 0.00001))+
    #  facet_grid(cols=vars(Stage))+
    theme_bw(base_size = 18)+
    theme(panel.grid.major = element_blank())+
    theme(strip.background = element_rect(fill = "#f0efef"))+
    theme(panel.grid.minor = element_blank())+
    ylab("Log2 CPM")+
    scale_fill_manual(values = cbPalette)}

lcpm$Stage <- factor(lcpm$Stage, levels=c("72h", "d5","Premel", "Mel"))

# Plot candidate genes 
# plot yellow-e counts 
my_plotCounts(lcpm, "jg1308")
# plot valkea counts 
my_plotCounts(lcpm, "jg1310")

# genes de in 3 stages 
intersect(filter(limma72, adj.P.Val<0.05)$genes, filter(limma5, adj.P.Val<0.05)$genes)
# gene DE in multiple stages
my_plotCounts(lcpm, "jg15945")

