library(optparse)
option_list <- list(
    make_option(c("-i", "--input_csv"), type = "character", help = "Input CSV file"),
    make_option(c("-l", "--log2fc"), type = "character", help = "Log2 fold change column name"),
    make_option(c("-p", "--p_val_adj"), type = "character", help = "Adjusted p-value column name"),
    make_option(c("-o", "--output_csv"), type = "character", help = "Output CSV file")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# parameters
input_csv <- opt$input_csv
key_log2fc <- opt$log2fc
key_p_val_adj <- opt$p_val_adj
output_csv <- opt$output_csv
# read data
res <- read.csv(input_csv, row.names = 1)
res <- rownames_to_column(res, var = "Gene") %>% filter(!is.na(Gene) & !duplicated(Gene))
gene <- pull(res, Gene)
geneList <- res[[key_log2fc]]
names(geneList) <- gene
head(geneList)
file_prefix <- str_remove(basename(input_csv), ".csv")
# load tomato database
Tomato <-loadDb(file="/data/users/lili10/online/script/Enrich/clusterProfiler/ITAG4.1_241207.OrgDb")
# run go enrich analysis
ego <- enrichGO(
    gene = gene,
    OrgDb = Tomato,
    keyType = 'GID',
    ont = 'ALL',
    qvalueCutoff = 0.01,
    pvalueCutoff = 0.01
)
ego_df <- as.data.frame(ego)
head(ego_df)
write.csv(ego_df, output_csv, row.names = FALSE)
# visualization
base_theme <- theme(
    plot.title = element_text(size = 16, hjust = 0.5, vjust = 1)
)
p1 <- dotplot(ego, showCategory = 30, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scale = "free", space = "free") +
    scale_colour_gradient(low = "red", high = "blue") +  # 颜色标度仅在此处添加
    labs(title = str_c("GO_DotPlot_", file_prefix)) +
    base_theme

p2 <- cnetplot(
        ego,
        foldChange = geneList,
        showCategory = 5,
        node_label = "category",
        layout = 'kk',
        color_edge = "category",
        color_category = "#279e68",
        max.overlaps = Inf
    ) +
    labs(title = str_c("GO_CnetPlot_", file_prefix)) +
    base_theme + 
    scale_color_gradient2(low='darkblue', mid = "#ff7f0e", high='darkred')

p3 <- cnetplot(
        ego,
        foldChange = geneList,
        showCategory = 5,
        layout = 'fr',
        node_label = "all",
        color_edge = "category",
        color_category = "#279e68",
        max.overlaps = Inf, 
        cex_label_gene = 0.4 
    ) +
    labs(title = str_c("GO_CnetPlot_circular_", file_prefix)) +
    base_theme + 
    scale_color_gradient2(low='darkblue', mid = "#ff7f0e", high='darkred')

p4 <- upsetplot(ego) +
    labs(title = str_c("GO_UpsetPlot_", file_prefix)) +
    base_theme

p5_1 <- heatplot(ego, showCategory=5)
p5_2 <- heatplot(ego, foldChange=geneList, showCategory=5)
p5 <- cowplot::plot_grid(p5_1, p5_2, ncol=1, labels=LETTERS[1:2]) + ggtitle(str_c("GO_heatPlot_", file_prefix)) + base_theme

# 计算GO term之间的相似性
ego_p <- pairwise_termsim(ego)

p6_1 <- treeplot(ego_p)
p6_2 <- treeplot(ego_p, cluster.params = list(method = "average"))
p6 <- aplot::plot_list(p6_1, p6_2, tag_levels='A')

p7 <- emapplot(ego_p, layout = 'fr') + ggtitle(str_c("GO_Enrichment_Map_", file_prefix)) + base_theme

# 计算生物学过程的富集结果
ego_BP <- enrichGO(
    gene = gene,
    OrgDb = Tomato,
    keyType = 'GID',
    ont = 'BP',
    qvalueCutoff = 0.05,
    pvalueCutoff = 0.05
)

p8 <- goplot(ego_BP, showCategory = 10) + ggtitle(str_c("GO_topGO_", file_prefix)) + base_theme 