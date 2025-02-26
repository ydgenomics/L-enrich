library(AnnotationDbi)
library(stringr) 
library(clusterProfiler)
library(tidyverse)
library(ggnewscale) 
library(tools)
library(ggpubr) 
library(enrichplot) 
library(ggupset) 
library(stats)
library(data.table)
library(lifecycle)
library(dplyr)
library(ggplot2)
library(topGO) 
library(Rgraphviz) 
library(RColorBrewer)
library(wordcloud) 
library(pathview) 
library(UpSetR) 
library(fgsea)

library(optparse)
option_list <- list(
    make_option(c("-i", "--input_csv"), type = "character", help = "Input CSV file"),
    make_option(c("-l", "--log2fc"), type = "character", help = "Log2 fold change column name"),
    make_option(c("-p", "--p_val_adj"), type = "character", help = "Adjusted p-value column name"),
    make_option(c("-o", "--output_csv"), type = "character", help = "Output CSV file")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_csv <- opt$input_csv
key_log2fc <- opt$log2fc
key_p_val_adj <- opt$p_val_adj
output_csv <- opt$output_csv

## 读取数据
res <- read.csv(input_csv, row.names = 1)
head(res)
## 基于筛选指标筛选差异基因
log2fc <- key_log2fc
p_val_adj <- key_p_val_adj
# Convert column names to symbols for use in dplyr
log2fc <- sym(log2fc)
p_val_adj <- sym(p_val_adj)
# Filter genes based on log2fc and p_val_adj
deg <- rownames_to_column(res, var = "Gene") %>%
    filter(abs(!!log2fc) > 1 & !!p_val_adj < 0.05) %>%
    filter(!is.na(Gene) & !duplicated(Gene)) %>%
    { 
        # 判断行数并排序截断
        if(nrow(.) > 3000) {
            arrange(., desc(abs(!!log2fc_sym))) %>% 
            slice_head(n = 3000)
        } else {
            .
        }
    }
deg$Gene <- sub("^gene:", "", deg$Gene)
head(deg1)
## 保存筛选结果
write.csv(deg, file = output_csv, row.names = FALSE)
