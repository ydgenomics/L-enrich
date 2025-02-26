gmt_go
gmt_kegg
input_csv

pathway_lines_go <- strsplit(readLines(gmt_go), "\t")
go_pathways <- lapply(pathway_lines_go, function(x) setNames(list(tail(x, -2)), head(x, 1)))
pathway_lines_kegg <- strsplit(readLines(gmt_kegg), "\t")
kegg_pathways <- lapply(pathway_lines_kegg, function(x) setNames(list(tail(x, -2)), head(x, 1)))
data <- read.csv(input_csv, row.names = 1) 
deg <- rownames_to_column(data1, var = "Gene") %>% filter(!is.na(Gene) & !duplicated(Gene))
head(deg)
go_ranks <- setNames(deg[["log2fc_1"]], deg$Gene)
go_ranks <- sort(go_ranks, decreasing = FALSE)

fgsea <- fgsea(pathways = go_pathways,
                  stats = go_ranks, 
                  eps = 0.0,
                  minSize = 15,
                  maxSize = 500)
head(fgsea)

p16 <- plotEnrichment(go_pathways[[head(fgsea[order(pval), ], 1)$pathway]],
                         go_ranks) + labs(title = head(fgsea[order(pval), ], 1)$pathway)

top_pathways_up <- fgsea[ES > 0][head(order(pval), n = 20), pathway]
top_pathways_down <- fgsea[ES < 0][head(order(pval), n = 20), pathway]
top_athways <- c(top_pathways_up, rev(top_pathways_down))
p17 <- plotGseaTable(go_pathways[top_pathways], go_ranks, fgsea, gseaParam = 0.5)

collapsed_pathways <- collapsePathways(fgsea[order(pval)][padj < 0.01], go_pathways, go_ranks)

main_pathways <- fgsea[pathway %in% collapsed_pathways$mainPathways][order(-NES), pathway]

p3a <- plotGseaTable(TomatoGOPathways[mainPathways1], TomatoGORanks1, fgseaRes1, gseaParam = 0.5)
p3a