library(optparse)
option_list <- list(
    make_option(c("-i", "--input_csv"), type = "character", help = "Input CSV file"),
    make_option(c("-l", "--key_log2fc"), type = "character", help = "Log2 fold change column name"),
    make_option(c("-p", "--key_p_val_adj"), type = "character", help = "Adjusted p-value column name"),
    make_option(c("-o", "--output_csv"), type = "character", help = "Output CSV file")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# parameters
input_csv <- opt$input_csv
key_log2fc <- opt$key_log2fc
key_p_val_adj <- opt$key_p_val_adj
output_csv <- opt$output_csv
# read data
res <- read.csv(input_csv, row.names = 1)
res <- rownames_to_column(res, var = "Gene") %>% filter(!is.na(Gene) & !duplicated(Gene))
gene <- pull(res, Gene)
geneList <- res[[key_log2fc]]
names(geneList) <- gene
head(geneList)
kegglist <- read.csv('/data/users/lili10/lili10_642e569efa3b4d56a57481c396194c66/online/input/Clusterprofiler/ITAG4.1_Gene_KEGGID_KEGGNAME_241225.csv')
gmt_go <- "/data/users/lili10/lili10_642e569efa3b4d56a57481c396194c66/online/input/Enrich/ITAG4.1_description_GO_Gene_241220.gmt"
gmt_kegg <- '/data/users/lili10/lili10_642e569efa3b4d56a57481c396194c66/online/input/Enrich/ITAG4.1_KEGGNAME_KEGGID_Gene_241225.gmt'
kegg2gene <- data.frame(kegglist$KEGGID,kegglist$Gene) 
kegg2name <- data.frame(kegglist$KEGGID,kegglist$KEGGNAME)
file_prefix <- str_remove(basename(input_csv), ".csv")
# run kk enrich analysis
kk <- enricher(gene, 
                TERM2GENE = kegg2gene,
                TERM2NAME = kegg2name,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)
kk_df <- as.data.frame(kk)
head(kk_df)
write.csv(kk_df, output_csv, row.names = FALSE)
# visualization
base_theme <- theme(
    plot.title = element_text(size = 16, hjust = 0.5, vjust = 1)
)
if (nrow(kk)==0){
    print('kk nrow == 0')
    blank_plot <- ggplot() +
        theme_void() +
        labs(title = "No Data Available")
    p9 <- blank_plot
    p10 <- blank_plot
    p11 <- blank_plot
    p12 <- blank_plot
    p13 <- blank_plot
    p14 <- blank_plot
    p15a <- blank_plot
} else {
    p9 <- dotplot(kk, showCategory = 30) + 
        labs(title = str_c("KEGG_DotPlot_", file_prefix)) + 
        base_theme
    
    p10 <- cnetplot(
        kk,
        foldChange = geneList,
        showCategory = 5,
        node_label = "category",
        layout = 'kk',
        color_edge = "category",
        color_category = "#279e68",
        max.overlaps = Inf
    ) +
    labs(title = str_c("KEGG_CnetPlot_", file_prefix)) +
    base_theme + 
    scale_color_gradient2(low='darkblue', mid = "#ff7f0e", high='darkred')
    
    p11 <- cnetplot(
        kk,
        foldChange = geneList,
        showCategory = 5,
        node_label = "all",
        layout = 'fr',
        color_category = "#279e68",
        color_edge = "category",
        max.overlaps = Inf,
        cex_label_gene = 0.6  
    ) +
    labs(title = str_c("KEGG_CnetPlot_circular_", file_prefix)) +
    base_theme + 
    scale_color_gradient2(low='darkblue', mid = "#ff7f0e", high='darkred')
    
    p12 <- upsetplot(kk) +
        labs(title = str_c("KEGG_UpsetPlot_", file_prefix)) +
        base_theme
    
    p13_1 <- heatplot(kk, showCategory=5)
    p13_2 <- heatplot(kk1, foldChange=geneList1, showCategory=5)
    p13 <- cowplot::plot_grid(p13_1, p13_2, ncol=1, labels=LETTERS[1:2]) + ggtitle(str_c("KEGG_heatPlot_", file_prefix)) + base_theme
    # 计算KEGG term之间的相似性
    kk_p <- pairwise_termsim(kk)
    p14_1 <- treeplot(kk_p)
    p14_2 <- treeplot(kk_p, cluster.params = list(method = "average"))
    p14 <- aplot::plot_list(p14_1, p14_2, tag_levels='A')
    
    p15 <- emapplot(kk_p, layout = 'fr') + ggtitle(str_c("KEGG_Enrichment_Map_", file_prefix)) + base_theme
}