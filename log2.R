
#R
rm(list = ls())
getwd()
setwd("/Users/chaochen/Desktop/Acanthus_r")
library(openxlsx)
library(readr)
library(ropls)
library(ggplot2)
library(dplyr)
file_path <-"/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-c/1.Data_Assess/all_group/all_sample_data.xlsx"
all_sample_data_metabolite0<- read.xlsx(file_path,1)

head(all_sample_data_metabolite0)

data_mat<- all_sample_data_metabolite0[, c(1,16:24)] 
data_mat[1:5, 1:10]

my_list <- list(
  "Hesperidin",
  "Bayin",
  "2-(2,4-dimethoxyphenyl)-3,5,6,7,8-pentamethoxychromen-4-one*",
  "3,3',4',5,6,7,8-heptamethoxyflavone*",
  "6'-Hydroxy-3,4,2',3',4',5'-Hexamethoxychalcone",
  "5,7,8,4'-Tetramethoxyflavone*"
)
print(my_list)


my_listmet<-all_sample_data_metabolite0[all_sample_data_metabolite0$Compounds%in%my_list,]
write.csv(my_listmet,"my_listmet.csv")

# 假设data_mat是你的数据框，行是代谢物，列是样本

row.names(data_mat) <- data_mat$Index
data_mat<- data_mat[, -1]
daixiewu<-data_mat["mws0036",]

write.csv(daixiewu,"daixiewu.csv")

data_mat[1:5, 1:9]

# 分组信息，按列顺序对应样本
group <- c(rep("L", 3), rep("M", 3), rep("H", 3))


calc_volcano <- function(data, group, grp1, grp2) {
  cols1 <- which(group == grp1)
  cols2 <- which(group == grp2)
  
  mean1 <- rowMeans(data[, cols1])
  mean2 <- rowMeans(data[, cols2])
  
  log2FC <- log2(mean2 + 1e-8) - log2(mean1 + 1e-8)
  
  p_values <- apply(data, 1, function(x) {
    x1 <- x[cols1]
    x2 <- x[cols2]
    # 判断是否有变异
    if (length(unique(x1)) == 1 && length(unique(x2)) == 1 && unique(x1) == unique(x2)) {
      return(NA)
    } else {
      tryCatch(
        t.test(x1, x2)$p.value,
        error = function(e) NA
      )
    }
  })
  
  result <- data.frame(
    Metabolite = rownames(data),
    log2FC = log2FC,
    p_value = p_values
  )
  
  result$negLog10P <- -log10(result$p_value)
  result$Significant <- "Not Significant"
  result$Significant[result$log2FC > 1 & result$p_value < 0.05] <- "Up"
  result$Significant[result$log2FC < -1 & result$p_value < 0.05] <- "Down"
  
  return(result)
}

plot_volcano <- function(df, title) {
  ggplot(df, aes(x = log2FC, y = negLog10P, color = Significant)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Significant" = "gray")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    labs(title = title, x = "log2(Fold Change)", y = "-log10(P value)")
}

# 计算M-L组差异
res_ML <- calc_volcano(data_mat, group, "L", "M")
print(res_ML %>% filter(Significant != "Not Significant"))
plot_volcano(res_ML, "Volcano plot: M vs L")

# 计算H-L组差异
res_HL <- calc_volcano(data_mat, group, "L", "H")
print(res_HL %>% filter(Significant != "Not Significant"))
plot_volcano(res_HL, "Volcano plot: H vs L")

res_ML_filt<-res_ML %>% filter(Significant != "Not Significant")
res_HL_filt<-res_HL %>% filter(Significant != "Not Significant")

all_sample_data_metabolite0_filteredML <- all_sample_data_metabolite0[all_sample_data_metabolite0$Index %in% res_ML_filt$Metabolite, ]
all_sample_data_metabolite0_filteredHL <- all_sample_data_metabolite0[all_sample_data_metabolite0$Index %in% res_HL_filt$Metabolite, ]
dim(all_sample_data_metabolite0_filteredML)
dim(all_sample_data_metabolite0_filteredHL)
head(all_sample_data_metabolite0_filteredHL)
# 按照Compounds和Metabolite列做left join合并
library(dplyr)
merged_df_ML <- left_join(all_sample_data_metabolite0_filteredML, res_ML_filt, by = c("Index" = "Metabolite"))
merged_df_HL <- left_join(all_sample_data_metabolite0_filteredHL, res_HL_filt, by = c("Index" = "Metabolite"))






head(merged_df_ML)
4

merged_df_ML$CAS

mean(merged_df_ML$CAS == "-", na.rm = TRUE)
mean(merged_df_HL$CAS == "-", na.rm = TRUE)

library(webchem)
library(KEGGREST)

write.xlsx(merged_df_ML, file = "merged_df_ML.xlsx", asTable = TRUE, overwrite = TRUE)
write.xlsx(merged_df_HL, file = "merged_df_HL.xlsx", asTable = TRUE, overwrite = TRUE)

head(merged_df_ML)


merged_df_ML$Index
merged_df_HL$Index

vip_greater_1<-read.csv("vip_greater_1.csv")
colnames(vip_greater_1)<-c("Index","Vip")
dim(vip_greater_1)
head(vip_greater_1)
dim(merged_df_ML)
dim(merged_df_HL)

merged_df_ML_logvip <- left_join(merged_df_ML, vip_greater_1, by = c("Index" = "Index"))
head(merged_df_ML_logvip)
write.csv(merged_df_ML_logvip,"merged_df_ML_logvip.csv")
merged_df_HL_logvip <- left_join(merged_df_HL, vip_greater_1, by = c("Index" = "Index"))
head(merged_df_HL_logvip)
write.csv(merged_df_HL_logvip,"merged_df_HL_logvip.csv")




merged_df_ML1 <- merged_df_ML[merged_df_ML$Index %in% vip_greater_1$Index, ]
merged_df_HL1 <- merged_df_HL[merged_df_HL$Index %in% vip_greater_1$Index, ]

dim(merged_df_ML1)
dim(merged_df_HL1)


# 安装RefMet包
#devtools::install_github("metabolomicsworkbench/RefMet")
# 加载RefMet包
library(RefMet)
# 代谢物名称向量（示例）
metabolite_names <- merged_df_ML1$Compounds
# 批量转换标准名
resultML <- refmet_map_df(metabolite_names)
# 查看转换结果
print(resultML)
write.xlsx(resultML, file = "resultML.xlsx", asTable = TRUE, overwrite = TRUE)

metabolite_nameshl <- merged_df_HL1$Compounds
# 批量转换标准名
resultHL <- refmet_map_df(metabolite_nameshl)
write.xlsx(resultHL, file = "resultHL.xlsx", asTable = TRUE, overwrite = TRUE)
####后面在https://www.metaboanalyst.ca 上面做。


head(resultML)
dim(resultML)
head(resultHL)
dim(resultHL)

# # 安装必要包（如果没有安装）
# # BiocManager::install("clusterProfiler")
# # BiocManager::install("DOSE")
# library(clusterProfiler)
# # 示例RefMet ID
# refmet_ids <- resultML$KEGG_ID
# refmet_ids <- refmet_ids[refmet_ids != "-"]

# refmet_ids

# # 富集分析
# cpd_enrich_result <- enrichKEGG(
#     gene      = refmet_ids,      # 你的代谢物KEGG ID
#     organism  = "cpd",         # 化合物数据库
#     minGSSize = 1,
#     pvalueCutoff = 1           # 可自行调整
# )
# # 查看结果
# head(as.data.frame(cpd_enrich_result))
# # 可视化

# clusterProfiler::dotplot(cpd_enrich_result, showCategory=10)

resultMLkegg<-read.csv("resultMLkegg.csv")
resultMLkegg

df_top10ML <- resultMLkegg %>%
  dplyr::arrange(desc(X.log10.p.)) %>%
  dplyr::slice(1:20)

ggplot(resultMLkegg, aes(x = Hits, y = reorder(X, Hits))) +
  geom_point(aes(size = Hits, color = Raw.p)) +  # 点大小对应Hits，颜色对应pvalue
  scale_color_gradient(low = "red", high = "blue") +  # pvalue颜色渐变，显著红色
  labs(title = "KEGG Pathway Enrichment Bubble Plot (M-L)",
       x = "Hits (Count)",
       y = "Pathway",
       color = "p-value",
       size = "Hits") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

df_top10ML$X

resultHLkegg<-read.csv("resultHLkegg.csv")

df_top10HL <- resultHLkegg %>%
  dplyr::arrange(desc(X.log10.p.)) %>%
  dplyr::slice(1:20)

ggplot(df_top10HL, aes(x = Hits, y = reorder(X, Hits))) +
  geom_point(aes(size = Hits, color = Raw.p)) +  # 点大小对应Hits，颜色对应pvalue
  scale_color_gradient(low = "red", high = "blue") +  # pvalue颜色渐变，显著红色
  labs(title = "KEGG Pathway Enrichment Bubble Plot(H-L)",
       x = "Hits (Count)",
       y = "Pathway",
       color = "p-value",
       size = "Hits") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

df_top10HL$X
df_top10ML$X


library(ggvenn)
set1 = df_top10HL$X
set2 = df_top10ML$X

data_list <- list(
  set1 = df_top10HL$X,
  set2 = df_top10ML$X
)

ggvenn(data_list, show_elements = FALSE)

# 输出交集和特有元素
common <- intersect(set1, set2)
unique_set1 <- setdiff(set1, set2)
unique_set2 <- setdiff(set2, set1)

list(
  Intersection = common,
  Unique_to_Set1 = unique_set1,
  Unique_to_Set2 = unique_set2
)

######大类统计下

merged_df_ML1$Class.II
merged_df_ML1$Class.I

# 统计Class.II列的频数
counts <- table(merged_df_ML1$Class.I)
# 绘制饼图
pie(counts, labels = names(counts), main = "M-L classification pie chart", col = rainbow(length(counts)))
# 绘饼图，不显示标签
# 统计Class.II列的频数

counts <- table(merged_df_HL1$Class.I)
# 绘制饼图
pie(counts, labels = names(counts), main = "H-L classification pie chart", col = rainbow(length(counts)))
# 绘饼图，不显示标签
# 统计Class.II列的频数

library(plotrix)
par(mfrow = c(1, 1))  # 一行两列
# 假设数据来自merged_df_ML1$Class.II
counts <- table(merged_df_ML1$Class.I)
# 去除为0的类别（防止多余标签）
counts <- counts[counts > 0]
# 创建标签，添加百分比，确保没有NA和长度匹配
labels <- paste0(names(counts), " ", round(100*counts/sum(counts), 1), "%")
# 绘制3D饼图
pie3D(counts, labels = labels, explode = 0.1, labelcex=0.7,col = terrain.colors(length(counts)),
      main = "Classification of M-L metabolites")
counts_sorted <- sort(counts)
counts_sorted 



counts <- table(merged_df_HL1$Class.I)
# 去除为0的类别（防止多余标签）
counts <- counts[counts > 0]
# 创建标签，添加百分比，确保没有NA和长度匹配
labels <- paste0(names(counts), " ", round(100*counts/sum(counts), 1), "%")
# 绘制3D饼图
pie3D(counts, labels = labels, explode = 0.1, labelcex=0.7,col = terrain.colors(length(counts)),
      main = "Classification of H-L metabolites")


counts
counts_sorted <- sort(counts)
counts_sorted 
