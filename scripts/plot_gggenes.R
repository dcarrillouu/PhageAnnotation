#

plot_ORF_callers_annot <- function(phanot, prod11, prodTAG, prodTGA, outfigure) {
  library(gggenes)
  library(ggplot2)
  library(dplyr)
  
  data1 <- read.csv(phanot, header = T, sep = "\t", check.names = F)  
  data2 <- read.csv(prod11, header = T, sep = "\t", check.names = F)  
  data3 <- read.csv(prodTAG, header = T, sep = "\t", check.names = F)  
  data4 <- read.csv(prodTGA, header = T, sep = "\t", check.names = F)  
  
  data <- rbind(data1, data2, data3, data4)
  data$direction <- ifelse(data$strand == "+", 1, -1)
  
  data %>% filter(caller == "tRNAscan-SE") %>% arrange(caller)
  
  a <-ggplot(data, aes(xmin = start, xmax = end, y = caller, forward = direction, fill = VOG_annot)) +
    geom_gene_arrow(show.legend = F) +
    facet_wrap(~ caller, scales = "free", ncol = 1) +
    #scale_fill_brewer(palette = "Set3") +
    theme_genes()
  
  
  ggsave(outfigure, a)
  
  
}

plot_ORF_callers_annot(snakemake@input[[1]],
                       snakemake@input[[2]],
                       snakemake@input[[3]],
                       snakemake@input[[4]],
                       snakemake@output[[1]])

