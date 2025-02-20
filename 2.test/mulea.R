library(ggplot2)
library(tidyverse)
library(dplyr)
library(mulea)
library(tidyverse)
library(optparse)
option_list <- list(
  make_option(c("--gene_csv"),
    type = "character", default = NULL,
    help = "input the csv of leiden_0.5"
  ),
  make_option(c("--preframe_csv"),
    type = "character", default = NULL,
    help = "input the pre-frame for ora"
  ),
  make_option(c("--species_gmt"),
    type = "character", default = NULL,
    help = "input the gmt for ora"
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


filepath <- paste0(opt$filepath, "_mulea")
dir.create(filepath)
setwd(filepath)

markers <- read.csv(opt$gene_csv, header = TRUE, stringsAsFactors = FALSE)
pre_frame <- read.csv(opt$preframe_csv, header = TRUE, stringsAsFactors = FALSE)
tf_ontology_filtered <- read_gmt(opt$species_gmt)
if (opt$species == "Cer") {markers$gene_id <- sub(".v2.1$", ".1.p", markers$names)}
if (opt$species == "Pog") {markers$gene_id <- gsub("_", "-", markers$names)}
head(markers)
for(i in unique(markers$group)){
    marker_subset <- filter(markers, group == i)
    length(marker_subset$gene_id)
    gene_list <- marker_subset %>% filter(pvals_adj < opt$minp)
    gene_list <- gene_list$gene_id
    # run mulea
    target_set <- gene_list
    background_set  <- marker_subset$gene_id
    ora_model <- ora(gmt = tf_ontology_filtered, 
                 element_names = target_set, 
                 background_element_names = background_set, 
                 p_value_adjustment_method = "eFDR", 
                 number_of_permutations = 10000,
                 nthreads = 2, 
                 random_seed = 1) 
    ora_results <- run_test(ora_model)
    n_rows <- ora_results %>% filter(eFDR < 0.05) %>% nrow()
    if (n_rows > 0) {
        ora_reshaped_results <- reshape_results(model = ora_model, model_results = ora_results, p_value_type_colname = "eFDR")
        length(ora_reshaped_results$ontology_id)
        ora_reshaped_results <- left_join(ora_reshaped_results, pre_frame %>% dplyr::select(GO, Description, Ontology), by = c("ontology_id" = "GO"))
        length(ora_reshaped_results$ontology_id)
        ora_reshaped_results$Name <- paste0(ora_reshaped_results$ontology_id,"_",ora_reshaped_results$Description)
        ora_reshaped_results$Name <- ifelse(nchar(ora_reshaped_results$Name) > 100, substr(ora_reshaped_results$Name, 1, 100), ora_reshaped_results$Name)
        ora_reshaped_results <- ora_reshaped_results %>% mutate(Ontology = case_when(Ontology == "molecular_function" ~ "MF", 
                                                                                     Ontology == "biological_process" ~ "BP", 
                                                                                     Ontology == "cellular_component" ~ "CC",TRUE ~ Ontology))
        #visual mulea
        pdf(paste0(i,"_eFDR.pdf"))
        p1 <- plot_lollipop(reshaped_results = ora_reshaped_results,ontology_id_colname = "Name",p_value_max_threshold = 0.05,p_value_type_colname = "eFDR")
        p4 <- plot_heatmap(reshaped_results = ora_reshaped_results,ontology_id_colname = "ontology_id",p_value_max_threshold = 0.05,p_value_type_colname = "eFDR")
        print(p1)
        print(p4)
        dev.off()
        ora_reshaped_results <- group_by(ora_reshaped_results,ontology_id) %>% summarise(eFDR = eFDR[1], Description = Description[1], Ontology = Ontology[1], Name = Name[1]) %>% ungroup()
        if (!"Count" %in% names(ora_reshaped_results)) {ora_reshaped_results$Count <- 0}
        if (!"GeneRatio" %in% names(ora_reshaped_results)) {ora_reshaped_results$GeneRatio <- 0}
        for (j in ora_reshaped_results$ontology_id) {
            matched_rows <- ora_results[ora_results$ontology_id == j, ]
            if (nrow(matched_rows) > 0) {
                ora_reshaped_results$Count[ora_reshaped_results$ontology_id == j] <- sum(matched_rows$nr_common_with_tested_elements)
                ora_reshaped_results$GeneRatio[ora_reshaped_results$ontology_id == j] <- sum(matched_rows$nr_common_with_tested_elements) / length(gene_list)
            }
        }
        head(arrange(ora_reshaped_results, eFDR))
        ora_reshaped_results <- ora_reshaped_results %>% filter(eFDR < 0.05)
        ora_reshaped_results <- ora_reshaped_results %>% arrange(desc(Count)) %>% head(60)
        # viusal2
        pdf(paste0(i,"_plot2.pdf"))
        plot2 <- ggplot(ora_reshaped_results, aes(y = GeneRatio, x = reorder(Name, GeneRatio))) + 
            geom_bar(stat = "identity", aes(fill = eFDR), width = 0.8) +  
            scale_fill_gradient(low = "red", high = "blue") +  
            facet_grid(Ontology ~ ., scales = "free", space = "free") +  
            coord_flip() + xlab("Name") + ylab("GeneRatio") + labs(title = paste0("Group ", i, " GO Terms Enrich")) + 
            theme(
                axis.text.x = element_text(size = 10), 
                axis.text.y = element_text(size = 5), 
                axis.title.x = element_text(size = 12),  
                axis.title.y = element_text(size = 12)) +
            geom_text(aes(label = Count), vjust = 0, size = 1.5) +
            scale_size_continuous(range = c(0.1, 3)) 
        print(plot2)
        dev.off()}
    #write.table(ora_reshaped_results, file = paste0(i,"_ora_reshaped_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    #ora_reshaped <- read.table("ora_reshaped_results.txt", header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE) 
}
