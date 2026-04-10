#R
rm(list = ls())
getwd()
setwd("/Users/chaochen/Desktop/Acanthus_r")
library(dplyr)
library(openxlsx)
file_path <-"/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-b/05.Quantitation/Unigene.count.annot.xlsx"

all_sample_data_RNA0<- read.xlsx(file_path,1)
head(all_sample_data_RNA0,5)
diff_genes_HL<-read.csv("DESeq2_diff_genes_H_vs_L.csv")
diff_genes_ML<-read.csv("DESeq2_diff_genes_M_vs_L.csv")

dim(diff_genes_HL)
dim(diff_genes_ML)

length(diff_genes_ML$X)

diff_geneML<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%diff_genes_ML$X,]

dim(diff_geneML)
diff_geneHL<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%diff_genes_HL$X,]

dim(diff_geneHL)
head(diff_geneHL)
library(clusterProfiler)
head(diff_geneHL,1)

library(dplyr)
library(tidyr)

packageVersion("clusterProfiler")

colnames(all_sample_data_RNA0)[colnames(all_sample_data_RNA0) == "ID"] <- "GeneID"
colnames(all_sample_data_RNA0)

#    /Users/chaochen/Desktop/Acanthus_r/ath_ko.txt

diff_geneHL <- diff_geneHL[!grepl("--", diff_geneHL$KEGG), ]
diff_geneML <- diff_geneML[!grepl("--", diff_geneML$KEGG), ]

dim(diff_geneHL)
dim(diff_geneML)

# 假设 all_sample_data_RNA0 是你的数据框，KEGG列是KO编号，先提取KO编号
ko_ids_raw <- diff_geneML$KEGG
# 从KO编号列清理得到纯KO编号（去除描述，只留下类似K10527）
ko_ids <- gsub(" .*", "", ko_ids_raw)
ko_ids
ko_ids_unique <- unique(ko_ids)
length(ko_ids)

# 读入KO到拟南芥基因ID映射，例如你从KEGG数据库或其他渠道获得的映射表，格式示例：
# ko2gene <- data.frame(KO = c("K10527", "Kxxxx", ...), 
#                       GeneID = c("AT1Gxxxx", "AT2Gxxxx", ...), stringsAsFactors=FALSE)
# 假设你有这个映射表，示例手动创建一个（请替换成真实数据）
ko2gene <- read.table("/Users/chaochen/Desktop/Acanthus_r/ath_ko.txt", header = FALSE, na.strings = c("NA"))
names(ko2gene)<-c("KO","GeneID")
head(ko2gene)
ko2gene$KO <- gsub("^ko:", "", ko2gene$KO)
ko2gene$GeneID <- gsub("^ath:", "", ko2gene$GeneID)
# 根据你的KO编号匹配得到对应的拟南芥基因ID
mapped_genes <- ko2gene %>% filter(KO %in% ko_ids) %>% pull(GeneID) %>% unique()
head(mapped_genes)
# 进行enrichKEGG富集分析
kegg_result <- enrichKEGG(
  gene = mapped_genes,
  organism = "ath",      # 拟南芥的KEGG物种代码
  keyType = "kegg",      # 表示输入基因ID符合KEGG基因ID格式（AT开头）
  pvalueCutoff = 0.05
)
# 查看富集结果
head(kegg_result)
#dotplot(kegg_result, showCategory=20)
library(ggplot2)
library(stringr)
barplot(kegg_result, showCategory=20, font.size=12, title="KEGG Enrichment Barplot of M-L") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))  # 设一个很大的width防止换行



####################hL
# 假设 all_sample_data_RNA0 是你的数据框，KEGG列是KO编号，先提取KO编号
ko_ids_raw <- diff_geneHL$KEGG
# 从KO编号列清理得到纯KO编号（去除描述，只留下类似K10527）
ko_ids <- gsub(" .*", "", ko_ids_raw)
ko_ids
ko_ids_unique <- unique(ko_ids)
length(ko_ids)

# 读入KO到拟南芥基因ID映射，例如你从KEGG数据库或其他渠道获得的映射表，格式示例：
# ko2gene <- data.frame(KO = c("K10527", "Kxxxx", ...), 
#                       GeneID = c("AT1Gxxxx", "AT2Gxxxx", ...), stringsAsFactors=FALSE)
# 假设你有这个映射表，示例手动创建一个（请替换成真实数据）
ko2gene <- read.table("/Users/chaochen/Desktop/Acanthus_r/ath_ko.txt", header = FALSE, na.strings = c("NA"))
names(ko2gene)<-c("KO","GeneID")
head(ko2gene)
ko2gene$KO <- gsub("^ko:", "", ko2gene$KO)
ko2gene$GeneID <- gsub("^ath:", "", ko2gene$GeneID)
# 根据你的KO编号匹配得到对应的拟南芥基因ID
mapped_genes <- ko2gene %>% filter(KO %in% ko_ids) %>% pull(GeneID) %>% unique()
head(mapped_genes)
# 进行enrichKEGG富集分析
kegg_resulthl <- enrichKEGG(
  gene = mapped_genes,
  organism = "ath",      # 拟南芥的KEGG物种代码
  keyType = "kegg",      # 表示输入基因ID符合KEGG基因ID格式（AT开头）
  pvalueCutoff = 0.05
)
# 查看富集结果
head(kegg_resulthl)

# 绘制柱状图和点图示例

barplot(kegg_resulthl, showCategory=20, font.size=12, title="KEGG Enrichment Barplot of H-L") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))  # 设一个很大的width防止换行


kegg_result$Description[1:20]
kegg_resulthl$Description[1:20]

intersection_genes <- intersect(kegg_result$Description[1:20], kegg_resulthl$Description[1:20])
print(intersection_genes)

desc1RNA <- kegg_result$Description[1:20]
desc2RNA <- kegg_resulthl$Description[1:20]

# 找出desc1中有但desc2中没有的元素，即特有成分
unique_to_desc1 <- setdiff(desc1RNA, desc2RNA)
unique_to_desc1
unique_to_desc2 <- setdiff(desc2RNA, desc1RNA)
unique_to_desc2

################代谢和基因的交集
resultHLkegg<-read.csv("resultHLkegg.csv")
df_top10HL <- resultHLkegg %>%
  dplyr::arrange(desc(X.log10.p.)) %>%
  dplyr::slice(1:20)
df_top10HL$X

resultMLkegg<-read.csv("resultMLkegg.csv")
df_top10ML <- resultMLkegg %>%
  dplyr::arrange(desc(X.log10.p.)) %>%
  dplyr::slice(1:20)
df_top10ML$X

intersection_genes <- intersect(desc1RNA, resultMLkegg)

print(intersect(desc1RNA, df_top10ML$X))
print(intersect(desc2RNA, df_top10HL$X))

