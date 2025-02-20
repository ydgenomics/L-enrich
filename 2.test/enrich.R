library(ggplot2)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Cthalictroides.eg.db)
length(keys(org.Cthalictroides.eg.db))
library(org.Pcirratum.eg.db)
length(keys(org.Pcirratum.eg.db))
library(tidyverse)
library(optparse)
option_list <- list(
  make_option(c("--gene_csv"),
    type = "character", default = NULL,
    help = "input the csv of leiden_0.5"
  ),
  make_option(c("--minp"),
    type = "numeric", default = NULL,
    help = "filter marker gene limited by min pvalue_adj"
  ),
  make_option(c("--species"),
    type = "character", default = NULL,
    help = "differ different species use different database"
  ),
  make_option(c("--filepath"),
    type = "character", default = NULL,
    help = "store all files"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

filepath <- paste0(opt$filepath, "_enrichGO")
dir.create(filepath)
setwd(filepath)


markers <- read.csv(opt$gene_csv, header = TRUE, stringsAsFactors = FALSE)
if (opt$species == "Cer") {markers$gene_id <- sub(".v2.1$", ".1.p", markers$names)}
if (opt$species == "Pog") {markers$gene_id <- gsub("_", "-", markers$names)}
head(markers)
for(i in unique(markers$group)){
    marker_subset <- filter(markers, group == i)
    length(marker_subset$gene_id)
    gene_list <- marker_subset %>% filter(pvals_adj < opt$minp)
    gene_list <- gene_list$gene_id
    length(gene_list)
    # run enrichGO
    if (opt$species == "Cer") {
    data <- enrichGO(gene = gene_list,OrgDb = org.Cthalictroides.eg.db,keyType = 'GID',ont = 'ALL',qvalueCutoff = 0.05,pvalueCutoff = 0.05)}
    if (opt$species == "Pog") {
    data <- enrichGO(gene = gene_list,OrgDb = org.Pcirratum.eg.db,keyType = 'GID',ont = 'ALL',qvalueCutoff = 0.05,pvalueCutoff = 0.05)}
    data <- as_tibble(data)
    if (nrow(data) > 0) {
        print("Data is not empty. Proceeding with the code.")
        length(data$ID)
        data$Name <- paste0(data$ID,"_",data$Description)
        #write.table(data, file = paste0(i,"_enrichGO_data.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        data_subset <- data %>% arrange(desc(Count)) %>% head(60)
        length(data_subset$ID)
        data_subset <- data_subset %>% mutate(GeneRatio = as.numeric(gsub("/.*", "", GeneRatio)) / as.numeric(gsub(".*/", "", GeneRatio)))
        data_subset$Name <- ifelse(nchar(data_subset$Name) > 100, substr(data_subset$Name, 1, 100), data_subset$Name)
        # viusal1
        if (length(data_subset$ID) > 0) {
        pdf(paste0(i,"_plot1.pdf"))
        plot1 <- ggplot(data_subset, aes(y = GeneRatio, x = reorder(Name, GeneRatio))) + 
            geom_bar(stat = "identity", aes(fill = p.adjust), width = 0.8) +  
            scale_fill_gradient(low = "red", high = "blue") +  
            facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +  
            coord_flip() + xlab("Name") + ylab("GeneRatio") + labs(title = paste0("Group ", i, " GO Terms Enrich")) + 
            theme(
                axis.text.x = element_text(size = 10), 
                axis.text.y = element_text(size = 5), 
                axis.title.x = element_text(size = 12),  
                axis.title.y = element_text(size = 12)) +
            geom_text(aes(label = Count), vjust = 0, size = 1.5) +
            scale_size_continuous(range = c(0.1, 3)) 
        print(plot1)
        dev.off()}
    } else {
        print("Data is empty. Skipping the code.")
    }
}
