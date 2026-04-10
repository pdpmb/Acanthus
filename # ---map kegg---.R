# ---map kegg---
# 创建数据框
library(ggplot2)
library(dplyr)
library(tidyr)
###读取RNA seq2 均一化的
normalized_counts<-  read.csv( "normalized_counts.csv")


head(normalized_counts)
# 原始数据
  ID = c("Cluster-16166.0", "Cluster-18867.0", "Cluster-22427.0", "Cluster-27723.10",
         "Cluster-50390.0", "Cluster-51648.0", "Cluster-55426.0", "Cluster-57456.0",
         "Cluster-58239.0", "Cluster-78165.6", "Cluster-82875.4", "Cluster-84503.6",
         "Cluster-86453.1")


df<-normalized_counts[normalized_counts$X%in%ID,]



library(ggplot2)
library(dplyr)
library(tidyr)

# Assume data frame 'df' with columns: X (gene), L1,L2,L3, M1,M2,M3, H1,H2,H3

# 1. Calculate means for each gene in groups L, M, H
library(dplyr)
library(tidyr)

# 把数据转成长格式
df_long <- df %>%
  pivot_longer(cols = -X, names_to = "sample", values_to = "expr") %>%
  mutate(group = case_when(
    grepl("^L", sample) ~ "L",
    grepl("^M", sample) ~ "M",
    grepl("^H", sample) ~ "H"
  ))

# 按基因拆分数据，生成列表
gene_list <- df_long %>%
  group_split(X)

# 创建一个空dataframe存放结果
anova_results <- data.frame(X = character(), p_value = numeric())

# 循环对每个基因做ANOVA
for (gene_data in gene_list) {
  gene_name <- unique(gene_data$X)
  aov_res <- aov(expr ~ group, data = gene_data)
  p <- summary(aov_res)[[1]][["Pr(>F)"]][1]
  anova_results <- rbind(anova_results, data.frame(X = gene_name, p_value = p))
}

# 计算每组均值用于绘图
plot_data <- df_long %>%
  group_by(X, group) %>%
  summarize(mean_expr = mean(expr, na.rm = TRUE), .groups = 'drop') %>%
  left_join(anova_results, by = "X")

# 绘图部分略，使用plot_data，按需要画条形图，并在标题加入p_value即可

library(ggplot2)
library(dplyr)

# 为了示范，保证group顺序是L, M, H
plot_data$group <- factor(plot_data$group, levels = c("L", "M", "H"))

# 对每个基因画一个柱状图，展示三组均值，标题加p值
unique_genes <- unique(plot_data$X)

for (gene in unique_genes) {
  gene_data <- filter(plot_data, X == gene)
  p_val <- unique(gene_data$p_value)
  p_val_text <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", signif(p_val, 3)))
  
  p <- ggplot(gene_data, aes(x = group, y = mean_expr, fill = group)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("L" = "#66c2a5", "M" = "#fc8d62", "H" = "#8da0cb")) +
    labs(title = paste0(gene, " 组间表达差异ANOVA ", p_val_text),
         x = "组别",
         y = "平均表达量") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
  print(p)  # 输出图形
}

##-------------一个图-----------
library(ggplot2)
library(dplyr)

# 保证group顺序
plot_data$group <- factor(plot_data$group, levels = c("L", "M", "H"))

# p值格式化，准备显示
pvals <- plot_data %>%
  distinct(X, p_value) %>%
  mutate(p_label = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 3))))

# 合并p值显示
plot_data <- left_join(plot_data, pvals, by = "X")

# 利用facet_wrap绘制所有基因分面图
ggplot(plot_data, aes(x = group, y = mean_expr, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), color = "black") +
  scale_fill_manual(values = c("L" = "#66c2a5", "M" = "#fc8d62", "H" = "#8da0cb")) +
  #facet_wrap(~ X, scales = "free_y", ncol = 3) +  # 每行3个图，根据需要调整ncol

  facet_wrap(~ X, scales = "free", ncol = 3) 
  labs(x = "组别", y = "平均表达量") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12)) +
  # 添加小标题为p值，使用strip标签下方添加辅助函数
  geom_text(data = pvals, aes(x = 2, y = Inf, label = p_label), inherit.aes = FALSE, vjust = 1.5, hjust = 0.5)


library(openxlsx)
file_path <-"/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-b/05.Quantitation/Unigene.count.annot.xlsx"

all_sample_data_RNA0<- read.xlsx(file_path,1)

kegg_ID<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%ID,]
kegg_ID
write.csv(kegg_ID,"kegg_ID.csv")



##  ----qPCR-----
my_qpcrlist <- list("Cluster-27723.10",  #蛋白
                "Cluster-63782.1",
                "Cluster-67187.22",
                "Cluster-73702.5",
                "Cluster-75870.16")

gene_list <- list(
  "Cluster-75620.679",
  "Cluster-67187.22",
  "Cluster-73702.5",
  "Cluster-63782.1",
  "Cluster-75870.16",
  "Cluster-27723.10",
  "Cluster-25138.24",
  "Cluster-25578.39",
  "Cluster-73620.0"
)

my_qpcr<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%my_qpcrlist,]
my_qpcr

my_qcount<-normalized_counts[normalized_counts$X%in%my_qpcrlist,]
my_qcount
write.csv(my_qpcr,"my_qpcr.csv")
normalized_counts

###寻找稳定基因
# 将行名设置为基因名
rownames(normalized_counts) <- normalized_counts$X

# 移除基因名列，只保留表达矩阵
expr_matrix <- as.matrix(normalized_counts[ , -1])

threshold <- 10
# 计算每行大于阈值的样本数是否超过总样本数的一半
genes_over_50pct <- apply(expr_matrix, 1, function(x) {
  sum(x > threshold) > (length(x) / 2)
})

# 筛选满足条件的基因（行）
filtered_expr_matrix <- expr_matrix[genes_over_50pct, ]

# 输出筛选后的矩阵
print(filtered_expr_matrix)

# 计算每个基因在样本中的表达方差（方差越小说明越稳定）
gene_variances <- apply(filtered_expr_matrix, 1, var)

# 根据方差排序，方差小的靠前
stable_genes <- sort(gene_variances)

# 输出最稳定的基因和对应方差值

stable_genes[1:10]
stable_genes[1]
top10_names <- names(stable_genes)[1:1000]
# 输出
print(top10_names)
normalized_counts[top10_names,]

top10_namesStablegenes<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%top10_names,]
top10_namesStablegenes
write.csv(top10_namesStablegenes,"top10_namesStablegenes.csv")



kankan<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%gene_list,]

kankan[,c("ID","NR")]
###取出序列
# 安装Biostrings包（如果没装过）
if (!requireNamespace("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")

library(Biostrings)

# 读取FASTA文件
fasta_file <- "/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-b/02.Assembly/Unigene.fa"
fasta_sequences <- readDNAStringSet(fasta_file)

# 获取FASTA序列的名称
seq_names <- names(fasta_sequences)

# 假设你的基因名列表如下


# 在FASTA序列中匹配基因名，提取对应序列
matched_seqs <- fasta_sequences[seq_names %in% gene_list]

# 查看提取出的序列
print(matched_seqs)

# 如果需要保存为新的FASTA文件
writeXStringSet(matched_seqs, filepath = "extracted_sequences.fa")

####设计引物
# 读取FASTA文件中的序列和名称
# 读取FASTA文件中的序列和名称
library(Biostrings)
fasta_file <- "extracted_sequences.fa"   # 请替换为你的FASTA文件路径
dna_seqs <- readDNAStringSet(fasta_file)

# 创建Primer3输入文件的函数
create_primer3_input <- function(seqs, input_file) {
  fileConn <- file(input_file, open = "w")
  for (i in seq_along(seqs)) {
    seq_name <- names(seqs)[i]
    sequence <- as.character(seqs[[i]])
    cat(
      paste0("SEQUENCE_ID=", seq_name, "\n"),
      "PRIMER_TASK=generic\n",
      "PRIMER_PICK_LEFT_PRIMER=1\n",
      "PRIMER_PICK_RIGHT_PRIMER=1\n",
      "PRIMER_OPT_SIZE=20\n",
      "PRIMER_MIN_SIZE=18\n",
      "PRIMER_MAX_SIZE=25\n",
      "PRIMER_OPT_TM=60.0\n",
      "PRIMER_MIN_TM=57.0\n",
      "PRIMER_MAX_TM=63.0\n",
      "PRIMER_MAX_POLY_X=3\n",
      "PRIMER_PRODUCT_SIZE_RANGE=80-150\n",
      paste0("SEQUENCE_TEMPLATE=", sequence, "\n"),
      "=\n",
      file = fileConn,
      sep = "",
      append = TRUE
    )
  }
  close(fileConn)
}

# 生成Primer3输入文件
primer3_input_file <- "primer3_input.txt"
create_primer3_input(dna_seqs, primer3_input_file)

# 系统调用primer3_core程序（需先安装并配置环境变量）
primer3_exe <- "/opt/homebrew/bin/primer3_core"  # 或完整路径

primer3_output_file <- "primer3_output.txt"
system_command <- paste(primer3_exe, primer3_input_file, ">", primer3_output_file)
system(system_command)

# 读取输出文件内容（原始文本）
primer3_results <- readLines(primer3_output_file)
cat(primer3_results, sep = "\n")

# 后续你可以编写代码解析primer3_results，提取引物信息

##  ------随机选三个基因
set.seed(123)
random_rows <- filtered_expr_matrix[sample(nrow(filtered_expr_matrix), 3), ]
random_rows 

kankan2random_rows<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%rownames(random_rows),]
kankan2random_rows

kankan2random_rows[,c("ID","NR")]
#                         L1        L2        L3        M1        M2        M3         H1         H2         H3
# Cluster-25138.24  18.88101  20.10640  47.63572   0.00000   0.00000   0.00000 268.119034 162.786099 249.126954
# Cluster-25578.39 523.51882 599.17073 519.71790 580.35654 597.17748 558.21016 342.449064 355.287120 274.400702
# Cluster-73620.0   30.03797  23.45747  57.40715  42.68663  51.82605  22.44061   5.309288   5.167813   2.407024

#                      ID                                                        NR
# 23726  Cluster-25138.24  LOB domain-containing protein 40-like [Salvia hispanica]
# 24655  Cluster-25578.39      hypothetical protein Pfo_013142 [Paulownia fortunei]
# 135723  Cluster-73620.0 hypothetical protein GOBAR_DD24731 [Gossypium barbadense]