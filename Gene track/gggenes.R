#### Gene track with gggenes ####
# Brien, Orteu et al. 2022

library(ggplot2)
library(gggenes) 

## WW gene track
genes<- read.csv("trioYannontation.dupreg.gtf", sep = "\t", header = FALSE)

whitegenes <- read.csv("genes_white_region.csv",  header = T)

yellgenes <- read.csv("genes_yellow_region.csv",  header = T)


# white genes scaffold 419
plot1 <- ggplot(whitegenes, aes(xmin = start, xmax = end, y = CHR, fill = gene, label = gene, forward=orientation)) +
    geom_gene_arrow(arrow_body_height = grid::unit(x = 6, units = "mm"), 
                    arrowhead_height = grid::unit(x = 9, units = "mm"),
                    arrowhead_width = grid::unit(x = 5, units = "mm")) +
    ggrepel::geom_text_repel(data = whitegenes %>% mutate(start = (start + end)/2), 
                             aes(x = start, y = CHR, label = gene), 
                             inherit.aes = F, nudge_y = -0.1, size=5, segment.color = NA)+  
    theme_bw(base_size = 18)+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank())+
    theme_genes()+
    theme(legend.position = "none", axis.text.y  = element_blank(), 
          axis.title = element_blank()) 

# yellow genes scaffold 206
plot2 <- ggplot(yellgenes, aes(xmin = start, xmax = end, y = CHR, fill = gene, label = gene, forward=orientation)) +
    geom_gene_arrow(arrow_body_height = grid::unit(x = 6, units = "mm"), 
                    arrowhead_height = grid::unit(x = 9, units = "mm"),
                    arrowhead_width = grid::unit(x = 5, units = "mm")) +
    ggrepel::geom_text_repel(data = yellgenes %>% mutate(start = (start + end)/2), 
                             aes(x = start, y = CHR, label = gene), 
                             inherit.aes = F, nudge_y = -0.1, size=5, segment.color = NA)+  
    theme_bw(base_size = 18)+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank())+
    theme_genes()+
    theme(legend.position = "none", axis.text.y  = element_blank(), 
          axis.title = element_blank()) 
