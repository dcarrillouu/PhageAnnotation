#

plot_ORF_callers_annot <- function(phanot, prod11, prodTAG, prodTGA, outfigure) {
  library(gggenes)
  library(ggplot2)
  library(dplyr)
  
  data1 <- read.csv(phanot, header = T, sep = "\t", na.strings=c("","NA")) 
  data2 <- read.csv(prod11, header = T, sep = "\t", na.strings=c("","NA")) 
  data3 <- read.csv(prodTAG, header = T, sep = "\t", na.strings=c("","NA")) 
  data4 <- read.csv(prodTGA, header = T, sep = "\t", na.strings=c("","NA")) 
  
  data <- rbind(data1, data2, data3, data4)
  data$direction <- ifelse(data$strand == "+", 1, -1)
  
  vogs_uniq <- sort(unique(data$VOG_annot))
  
  if (length(vogs_uniq) == 0) {
    col_vector <- "white"
  } else {
    col_vector <- unname(distinctColorPalette(length(vogs_uniq)))
  }
  
  if ("unknown" %in% vogs_uniq) {
    col_vector <- replace(col_vector, match("unknown", vogs_uniq), "gray" )
  }
  
  a <-ggplot(data, aes(xmin = start, xmax = end, y = caller, forward = direction, fill = VOG_annot)) +
    geom_gene_arrow() +
    facet_wrap(~ caller, scales = "free", ncol = 1) +
    scale_fill_manual(values=col_vector, na.value = "white") +
    #scale_fill_brewer() +
    theme_genes() +
    theme(legend.text = element_text(color = "gray13", size = 8),
          legend.position = "bottom",
          legend.direction = "horizontal") 
  
  
  ggsave(outfigure, a, height= 10, width= 30, units = "cm")
  
  
}

plot_ORF_callers_annot(snakemake@input[[1]],
                       snakemake@input[[2]],
                       snakemake@input[[3]],
                       snakemake@input[[4]],
                       snakemake@output[[1]])

