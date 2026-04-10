

#R

#做下分析core_net
rm(list = ls())
getwd()
setwd("/Users/chaochen/Desktop/Acanthus_r")
library(dplyr)
library(openxlsx)

core_net<-read.csv("core_net.csv")
core_net

core_net$ID
#R

library(dplyr)
library(openxlsx)
file_path <-"/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-b/05.Quantitation/Unigene.count.annot.xlsx"

all_sample_data_RNA0<- read.xlsx(file_path,1)
dim(all_sample_data_RNA0)
head(all_sample_data_RNA0,5)
diff_genes_HL<-read.csv("DESeq2_diff_genes_H_vs_L.csv")
diff_genes_ML<-read.csv("DESeq2_diff_genes_M_vs_L.csv")

length(diff_genes_ML$X)

diff_geneML<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%diff_genes_ML$X,]

dim(diff_geneML)
diff_geneHL<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%diff_genes_HL$X,]
dim(diff_geneHL)

library(clusterProfiler)
head(diff_geneHL,1)

library(dplyr)
library(tidyr)
library(clusterProfiler)
packageVersion("clusterProfiler")

colnames(all_sample_data_RNA0)[colnames(all_sample_data_RNA0) == "ID"] <- "GeneID"

colnames(all_sample_data_RNA0)

# 拆分GO列：分号切分成多条，再取第一个字段（GO term）
library(dplyr)
library(tidyr)
library(stringr)
library(GO.db)

# 预构建GO到level映射表
library(parallel)
# all_go <- keys(GO.db)
# get_level <- function(goid) {
#   if (!goid %in% all_go) return(NA_integer_)
#   ont <- Ontology(GOTERM[[goid]])
#   ancestors <- switch(ont,
#                       BP = as.list(GOBPANCESTOR)[[goid]],
#                       MF = as.list(GOMFANCESTOR)[[goid]],
#                       CC = as.list(GOCCANCESTOR)[[goid]],
#                       NULL)
#   if (is.null(ancestors)) 1L else length(ancestors) + 1L
# }

# num_cores <- detectCores() - 1  # 留一个核心给系统
# levels <- mclapply(all_go, get_level, mc.cores = num_cores)

# # 把结果合成数据框
# go_level_map <- data.frame(GOterm = all_go, level = unlist(levels))

# write.csv(go_level_map,"go_level_map.csv")
go_level_map<-read.csv("go_level_map.csv")
# 对你的数据处理
# go_long <- all_sample_data_RNA0 %>%
#   separate_rows(GO, sep = ";") %>%
#   mutate(GO = trimws(GO)) %>%
#   filter(GO != "--") %>%
#   mutate(
#     GOterm = str_extract(GO, "GO:\\d{7}"),
#     description = str_trim(str_replace(GO, "GO:\\d{7},?", ""))
#   ) %>%
#   filter(!is.na(GOterm) & GOterm != "") %>%
#   distinct() %>%
#   left_join(go_level_map, by="GOterm")

# library(dplyr)
# library(tidyr)
# library(stringr)
# go_longs<-go_long %>% dplyr::select(GeneID, GOterm, description,level) %>%
#   filter(!is.na(GOterm) & GOterm != "") %>%
#   distinct()


# names(go_longs)

# write.csv(go_longs,"go_longs.csv")

go_longs<-read.csv("go_longs.csv")
head(go_longs)
go_level_2 <- go_longs %>%
  filter(level %in% c(6,7))

dim(go_level_2)

# 你的差异基因列表

gene_list<-core_net$ID  ################这里换

#TERM2GENE = go_long %>% select(GOterm, ID) %>% as.data.frame()
# 做GO富集分析，TERM2GENE格式为 (GOterm, GeneID)

ego <- enricher(gene = gene_list,
                TERM2GENE = go_level_2[, c("description", "GeneID")],
                pAdjustMethod = "BH",
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05)
head(ego)

# TERM2GENE <- go_long %>% select(GOterm, GeneID) %>% as.data.frame()
# intersect_genes <- intersect(gene_list, unique(TERM2GENE$GeneID))
ego@result
# 可视化
# library(enrichplot)
# names(ego@result) 
# format_scientific <- function(x, digits = 2) {
#   formatC(x, format = "e", digits = digits - 1)
# }

# ego@result$p.adjust <- format_scientific(ego@result$p.adjust, digits = 2)
# head(ego@result$p.adjust)

# ego@result$p.adjust

dotplot(ego, showCategory=20, font.size=12, title="GO Enrichment Dotplot of H-L") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))  # 设一个很大的width防止换行


#orderBy="p.adjust"
head(ego@result,20)
names(ego@result)
H_LDescription<-ego@result$Description[1:20]

ego
####乱图
library(clusterProfiler)
library(enrichplot)

# ego是enrichResult对象
cnetplot(ego, showCategory=20, foldChange=NULL, circular=TRUE, colorEdge=TRUE) +
  ggplot2::labs(title="GO Enrichment Circle Network Plot")
 
ego@result

# 从ego对象提取GO term与基因关联表
df <- ego@result %>%
  dplyr::select(ID, Description, geneID) %>%
  rename(term = ID, geneID = geneID) %>%
  mutate(geneID = strsplit(as.character(geneID), split = "/")) %>%
  unnest(geneID)

# 去除NA并转字符
df <- df %>% filter(!is.na(term), !is.na(geneID)) %>%
  mutate(term = as.character(term), geneID = as.character(geneID))

# 构建邻接矩阵
mat <- table(df$term, df$geneID)

# 绘制弦图
chordDiagram(mat, annotationTrack = "grid", preAllocateTracks = 1)

# 添加标题（用base R函数）
title("GO Enrichment Chord Diagram")