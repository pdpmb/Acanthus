
#R
rm(list = ls())
getwd()
setwd("/Users/chaochen/Desktop/Acanthus_r")
library(dplyr)
library(openxlsx)
file_path <-"/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-b/05.Quantitation/Unigene.count.annot.xlsx"
all_sample_data_RNA0<- read.xlsx(file_path,1)

head(all_sample_data_RNA0,5)

all_sample_data_RNA <- all_sample_data_RNA0[, c(1, 2:10)]  # 取第1、3、5列
head(all_sample_data_RNA)
all_sample_data_RNA[100:105,]
library(DESeq2)
packageVersion("DESeq2")
# 读入原始计数数据，行为基因，列为样本

row.names(all_sample_data_RNA) <- all_sample_data_RNA$ID
countData<- all_sample_data_RNA[, -1]

head(countData)

# 构建样本分组信息，确保顺序和count矩阵列顺序一致，比如这里是L组和M组3个样本各3个
colData <- data.frame(
  row.names = colnames(countData),
  group = factor(c(rep("L",3), rep("M",3),rep("H",3)))  # 修改对应的样本分组，如L组与M组对比
)
colData
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~group)
#过滤掉低表达基因（可选，但建议）
dds <- dds[rowSums(counts(dds)) > 1, ]

dds <- DESeq(dds)
############ 这里比较M组 vs L组
res_ML <- results(dds, contrast=c("group", "M", "L"))  
summary(res_ML )
# 筛选显著差异基因，常用阈值padj < 0.05且|log2FoldChange| > 1
diff_genes_ML  <- subset(res_ML , padj < 0.01 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(diff_genes_ML), file="DESeq2_diff_genes_M_vs_L.csv")
# 统计上调基因（log2FoldChange > 0）
num_upregulated <- sum(diff_genes_ML$log2FoldChange > 0, na.rm = TRUE)

# 统计下调基因（log2FoldChange < 0）
num_downregulated <- sum(diff_genes_ML$log2FoldChange < 0, na.rm = TRUE)
# 输出结果
cat("ML上调基因数量:", num_upregulated, "\n")
cat("ML下调基因数量:", num_downregulated, "\n")



res_HL <- results(dds, contrast=c("group", "H", "L"))
summary(res_HL)
# Filter for significantly different genes, commonly using threshold p-adjusted < 0.05 and |log2FoldChange| > 1
diff_genes_HL <- subset(res_HL, padj < 0.01 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(diff_genes_HL), file="DESeq2_diff_genes_H_vs_L.csv")
num_upregulated <- sum(diff_genes_HL$log2FoldChange > 0, na.rm = TRUE)

# 统计下调基因（log2FoldChange < 0）
num_downregulated <- sum(diff_genes_HL$log2FoldChange < 0, na.rm = TRUE)
# 输出结果
cat("ML上调基因数量:", num_upregulated, "\n")
cat("ML下调基因数量:", num_downregulated, "\n")


library(DESeq2)
library(pheatmap)

# 假设dds是你用DESeq2分析的对象

# 1. 进行rlog转换（对计数数据进行归一化和变换）
rld <- rlog(dds, blind=FALSE)

# 2. 提取差异表达分析结果，筛选显著差异基因（示例以padj<0.05）
res <- results(dds)
sig_genes <- rownames(subset(res, padj < 0.05 & !is.na(padj)))

# 3. 提取显著基因的表达矩阵
expr_mat <- assay(rld)[sig_genes, ]

# 4. 对表达矩阵按行标准化（基因）
expr_mat_scaled <- t(scale(t(expr_mat)))

# 5. 构建样本分组的注释信息（根据你的分组设计）
annotation_col <- data.frame(
  Group = colData(dds)$group
)
rownames(annotation_col) <- colnames(expr_mat_scaled)

# 6. 绘制热图
pheatmap(expr_mat_scaled,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         annotation_col=annotation_col,
         show_rownames=FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100))



## ---提取显著差异基因的随机。
set.seed(123)
expr_mat_scaled

set.seed(12545)
random_rows <- expr_mat_scaled[sample(nrow(expr_mat_scaled), 3), ]
random_rows 
kankan2random_rows<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%rownames(random_rows),]
kankan2random_rows

kankan2random_rows[,c("ID","NR")]
# #### 随机的结果
# 1731     Cluster-11016.0 PREDICTED: calcium-dependent protein kinase 29 [Tarenaya hassleriana]
# 169838   Cluster-81632.3                  hypothetical protein SASPL_105552 [Salvia splendens]
# 176826  Cluster-83147.111                  hypothetical protein SASPL_146904 [Salvia splendens]
