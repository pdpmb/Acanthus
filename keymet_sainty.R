#R
rm(list = ls())
getwd()
setwd("/Users/chaochen/Desktop/Acanthus_r")
library(dplyr)
library(openxlsx)

library(tidyr)

######################做两个自己的分析
merged_met_ML_logvip<-read.csv("merged_df_ML_logvip.csv")
merged_met_HL_logvip<-read.csv("merged_df_HL_logvip.csv")

names(merged_met_ML_logvip)

filtered_met_ML <- merged_met_ML_logvip %>%
  filter(!is.na(Vip)) %>%         # 筛选VIP不为NA的行
  dplyr::select(2, 17:25)                # 选择第2列和第17到25列

head(filtered_met_ML)
filtered_met_HL <- merged_met_HL_logvip %>%
  filter(!is.na(Vip)) %>%         # 筛选VIP不为NA的行
  dplyr::select(2, 17:25)                # 选择第2列和第17到25列
head(filtered_met_HL)


#####################################

rownames(filtered_met_HL)<-filtered_met_HL[,1]
filtered_met_HL<-filtered_met_HL[,-1]
filtered_met_HL<-t(filtered_met_HL)
filtered_met_HL[1:5,1:5]

# Autoscaling归一化（每列减去均值，除以标准差）

metabolite_data <- scale(filtered_met_HL)

# 添加盐浓度信息作为注释
annotation <- data.frame(
  Salt_Concentration = c(rep(0.18,3), rep(0.44,3), rep(0.58,3))
)
rownames(annotation) <- rownames(metabolite_data)

salinity<-annotation$Salt_Concentration
##################{
# 计算代谢物与盐度的相关性
cor_results <- apply(metabolite_data, 2, function(x){
  cor.test(x, salinity, method = "spearman",exact = FALSE) # 斯皮尔曼相关
})

# 提取相关系数和p值
cor_coefficients <- sapply(cor_results, function(x) x$estimate)
p_values <- sapply(cor_results, function(x) x$p.value)

# 组合结果为数据框
correlation_df <- data.frame(
  metabolite = names(cor_coefficients),
  spearman_rho = cor_coefficients,
  p_value = p_values
)
# 根据p值筛选显著相关代谢物（例如p < 0.05）
significant_metabolites <- correlation_df[correlation_df$p_value < 0.05, ]
# 按相关系数绝对值降序排序
significant_metabolites <- significant_metabolites[order(-abs(significant_metabolites$spearman_rho)), ]

# 输出相关性最高的10个代谢物
top20_metabolites <- head(significant_metabolites, 20)
print(top20_metabolites)
high_corr_metabolites <- significant_metabolites[
  abs(significant_metabolites$spearman_rho) > 0.95, 
]
# 输出结果
print(high_corr_metabolites)
dim(high_corr_metabolites)
metabolite_data[1:5,1:5]

high_corr_metabolites$metabolite <- sub("\\.rho$", "", high_corr_metabolites$metabolite)

data_scaled<-metabolite_data[, colnames(metabolite_data) %in% high_corr_metabolites$metabolite]

length(metabolite_data)
dim (metabolite_data)

colnames(metabolite_data)
################}

colnames(metabolite_data)
high_corr_metabolites$metabolite
write.csv(high_corr_metabolites,"high_corr_metabolites_HL.csv")

# 1. 按顺序将 data_scaled 列名找到对应的新名称
new_colnames <- colnames(data_scaled)  # 先拿当前列名向量

# 2. 找到这些列名对应在 merged_met_ML_logvip$Index 中的位置
match_idx <- match(new_colnames, merged_met_HL_logvip$Index) 

# 3. 用对应位置的 Compounds 替换，NA表示没有匹配上的保持原名
#new_colnames[!is.na(match_idx)] <- merged_met_ML_logvip$Compounds[match_idx[!is.na(match_idx)]]
data_scaled_filtered <- data_scaled[, !is.na(match_idx)]
dim(data_scaled_filtered)
# 替换为对应的Compounds列名
new_colnames <- merged_met_HL_logvip$Compounds[match_idx[!is.na(match_idx)]]
new_colnames

# 4. 赋值回 data_scaled
colnames(data_scaled_filtered) <- new_colnames

# 查看结果
print(colnames(data_scaled_filtered))


# 绘制热图
library(ComplexHeatmap)

ht<-Heatmap(data_scaled_filtered,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        column_names_rot = 45,
        column_names_side = "bottom",
        column_names_gp = grid::gpar(fontsize = 8),  # 字体小一点
        heatmap_legend_param = list(title = "Autoscaled metabolites"))
merged_met_HL_logvip
head(merged_met_HL_logvip)[1:6]
draw(ht)
# Add title on top
grid::grid.text("H-L group key differential metabolites", y = unit(1, "npc") - unit(2, "mm"), gp = grid::gpar(fontsize = 14, fontface = "bold"))


########################ML########################
############
#####################################

rownames(filtered_met_ML)<-filtered_met_ML[,1]
filtered_met_ML<-filtered_met_ML[,-1]
filtered_met_ML<-t(filtered_met_ML)
filtered_met_ML[1:5,1:5]

# Autoscaling归一化（每列减去均值，除以标准差）

metabolite_data <- scale(filtered_met_ML)

# 添加盐浓度信息作为注释
annotation <- data.frame(
  Salt_Concentration = c(rep(0.18,3), rep(0.44,3), rep(0.58,3))
)
rownames(annotation) <- rownames(metabolite_data)

salinity<-annotation$Salt_Concentration
##################{
# 计算代谢物与盐度的相关性
cor_results <- apply(metabolite_data, 2, function(x){
  cor.test(x, salinity, method = "spearman",exact = FALSE) # 斯皮尔曼相关
})

# 提取相关系数和p值
cor_coefficients <- sapply(cor_results, function(x) x$estimate)
p_values <- sapply(cor_results, function(x) x$p.value)

# 组合结果为数据框
correlation_df <- data.frame(
  metabolite = names(cor_coefficients),
  spearman_rho = cor_coefficients,
  p_value = p_values
)
# 根据p值筛选显著相关代谢物（例如p < 0.05）
significant_metabolites <- correlation_df[correlation_df$p_value < 0.05, ]
# 按相关系数绝对值降序排序
significant_metabolites <- significant_metabolites[order(-abs(significant_metabolites$spearman_rho)), ]

# 输出相关性最高的10个代谢物
top20_metabolites <- head(significant_metabolites, 20)
print(top20_metabolites)
high_corr_metabolites <- significant_metabolites[
  abs(significant_metabolites$spearman_rho) > 0.95, 
]
# 输出结果
print(high_corr_metabolites)
dim(high_corr_metabolites)
metabolite_data[1:5,1:5]

high_corr_metabolites$metabolite <- sub("\\.rho$", "", high_corr_metabolites$metabolite)

data_scaled<-metabolite_data[, colnames(metabolite_data) %in% high_corr_metabolites$metabolite]

length(metabolite_data)
dim (metabolite_data)

colnames(metabolite_data)
################}

colnames(metabolite_data)
high_corr_metabolites$metabolite
write.csv(high_corr_metabolites,"high_corr_metabolites_ML.csv")

# 1. 按顺序将 data_scaled 列名找到对应的新名称
new_colnames <- colnames(data_scaled)  # 先拿当前列名向量

# 2. 找到这些列名对应在 merged_met_ML_logvip$Index 中的位置
match_idx <- match(new_colnames, merged_met_ML_logvip$Index) 

# 3. 用对应位置的 Compounds 替换，NA表示没有匹配上的保持原名
#new_colnames[!is.na(match_idx)] <- merged_met_ML_logvip$Compounds[match_idx[!is.na(match_idx)]]
data_scaled_filtered <- data_scaled[, !is.na(match_idx)]
dim(data_scaled_filtered)
# 替换为对应的Compounds列名
new_colnames <- merged_met_ML_logvip$Compounds[match_idx[!is.na(match_idx)]]
new_colnames

# 4. 赋值回 data_scaled
colnames(data_scaled_filtered) <- new_colnames

# 查看结果
print(colnames(data_scaled_filtered))

# 绘制热图
library(ComplexHeatmap)

ht<-Heatmap(data_scaled_filtered,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        column_names_rot = 45,
        column_names_side = "bottom",
        column_names_gp = grid::gpar(fontsize = 8),
        heatmap_legend_param = list(title = "Autoscaled metabolites"))

draw(ht)
# Add title on top
grid::grid.text("M-L group key differential metabolites", y = unit(1, "npc") - unit(2, "mm"), gp = grid::gpar(fontsize = 14, fontface = "bold"))


