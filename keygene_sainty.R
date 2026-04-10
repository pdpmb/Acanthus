
## ------------keyDEgene---
###########################
###########################
###########################
rm(list = ls())
getwd()
setwd("/Users/chaochen/Desktop/Acanthus_r")
library(dplyr)
library(openxlsx)

library(tidyr)
library(dplyr)
library(openxlsx)
file_path <-"/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-b/05.Quantitation/Unigene.count.annot.xlsx"

all_sample_data_RNA0<- read.xlsx(file_path,1)
head(all_sample_data_RNA0,5)
diff_genes_HL<-read.csv("DESeq2_diff_genes_H_vs_L.csv")
diff_genes_all<-read.csv("DESeq2_diff_genes_M_vs_L.csv")
diff_genes_HL
diff_genes_all

dim(diff_genes_HL)
dim(diff_genes_all)

diff_genes_all<-merge(diff_genes_HL, diff_genes_all, by = "X", all = TRUE)
dim(diff_genes_all)

dim(diff_genes_HL)
dim(diff_genes_all)

length(diff_genes_all$X)

diff_geneall<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%diff_genes_all$X,]
dim(diff_geneall)
diff_geneall
names(diff_geneall)
filtered_rna_all <- diff_geneall %>%
  #filter(!is.na(Vip)) %>%         # 筛选VIP不为NA的行
  dplyr::select(1:10)                # 选择第2列和第17到25列
filtered_rna_all

rownames(filtered_rna_all)<-filtered_rna_all[,1]
filtered_rna_all<-filtered_rna_all[,-1]

sampleInfo <- data.frame(
  row.names = colnames(filtered_rna_all),
  condition = factor(c(rep("L",3), rep("M",3),rep("H",3)))  # 修改对应的样本分组，如L组与M组对比
)
sampleInfo

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = filtered_rna_all,
                              colData = sampleInfo,
                              design = ~ condition)  # condition为分组变量名

dds <- DESeq(dds)

# 提取归一化后的count值
normalized_counts <- counts(dds, normalized=TRUE)

# 保存归一化count结果（可用于关联分析）
write.csv(normalized_counts, file = "normalized_counts.csv")

head(normalized_counts)

filtered_rna_all<-t(normalized_counts)
filtered_rna_all[1:5,1:5]

rna_data <- scale(filtered_rna_all)
head(rna_data)


# 添加盐浓度信息作为注释
annotation <- data.frame(
  Salt_Concentration = c(rep(0.18,3), rep(0.44,3), rep(0.58,3))
)
rownames(annotation) <- rownames(rna_data)

salinity<-annotation$Salt_Concentration
salinity
salinity_scaled <- scale(salinity)
salinity_scaled
##################{
# 计算代谢物与盐度的相关性
cor_results <- apply(rna_data, 2, function(x){
  cor.test(x, salinity, method = "spearman",exact = FALSE) # 斯皮尔曼相关
})

# 提取相关系数和p值
cor_coefficients <- sapply(cor_results, function(x) x$estimate)
p_values <- sapply(cor_results, function(x) x$p.value)

# 组合结果为数据框
correlation_df <- data.frame(
  rna = names(cor_coefficients),
  spearman_rho = cor_coefficients,
  p_value = p_values
)

# 根据p值筛选显著相关代谢物（例如p < 0.05）
significant_rna <- correlation_df[correlation_df$p_value < 0.05, ]


# 按相关系数绝对值降序排序
significant_rna <- significant_rna[order(-abs(significant_rna$spearman_rho)), ]

significant_rna

# 输出相关性最高的10个代谢物
top20_rna <- head(significant_rna, 40)

#####我这里要调整 调为0.8
print(top20_rna)
high_corr_rna <- significant_rna[
  abs(significant_rna$spearman_rho) > 0.8, 
]
# 输出结果
print(high_corr_rna)
dim(high_corr_rna)



high_corr_rna$rna <- sub("\\.rho$", "", high_corr_rna$rna)

data_scaled<-rna_data[, colnames(rna_data) %in% high_corr_rna$rna]

length(rna_data)
dim (data_scaled)



t(data_scaled)
################}
colnames(rna_data)
high_corr_rna$rna


# 绘制热图
library(ComplexHeatmap)


set.seed(125)  # 保证每次运行随机选择结果一致
selected_genes <- sample(colnames(data_scaled), 25)
data_scaled
# 提取这100个基因的表达数据
data_subset <- data_scaled[,selected_genes]


ht<-Heatmap(data_subset,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        column_names_rot = 45,
        column_names_side = "bottom",
        column_names_gp = grid::gpar(fontsize = 8),  # 字体小一点
        heatmap_legend_param = list(title = "Autoscaled rna"))
draw(ht)
# Add title on top
grid::grid.text("Key DE genes", y = unit(1, "npc") - unit(2, "mm"), gp = grid::gpar(fontsize = 14, fontface = "bold"))

print(high_corr_rna)

high_corr_rna$rna
length(high_corr_rna$rna)

library(openxlsx)
file_path <-"/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-b/05.Quantitation/Unigene.count.annot.xlsx"

all_sample_data_RNA0<- read.xlsx(file_path,1)

names(all_sample_data_RNA0)

all_key_rna<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%high_corr_rna$rna,]

dim(all_key_rna)
write.csv(all_key_rna,"all_key_rna.csv")

all_key_rna["Cluster-80458.16",]

all_key_rna_second<-subset(all_key_rna, grepl("Biosynthesis of secondary metabolites", KEGG.pathway, ignore.case = TRUE))
names(all_key_rna_second)

all_key_rna_second[, c("ID", "KEGG.pathway")]

write.csv(all_key_rna_second[, c("ID", "KEGG.pathway")],"all_key_rna_second.csv")


write.xlsx(all_key_rna_second[, c("ID", "KEGG.pathway")], file = "TableS3.all_key_rna_second.xlsx", sheetName = "Sheet1", names = FALSE)

dim(all_key_rna_second)
##  -------------rna and metalite r?------------------
###########################
###########################
###########################
all_key_rna_second$ID
Rna_matrix0<-t(data_scaled)
Rna_matrix<-Rna_matrix0[rownames(Rna_matrix0)%in%all_key_rna_second$ID,]

Rna_matrix
dim(Rna_matrix)
dim(all_key_rna_second)
data_scaled[1:5,1:5]



dim(data_scaled)
high_corr_metabolites_ML<-read.csv("high_corr_metabolites_ML.csv")
high_corr_metabolites_HL<-read.csv("high_corr_metabolites_HL.csv")
selected_data_HL <- high_corr_metabolites_HL[high_corr_metabolites_HL$spearman_rho > 0, ]

selected_data_HL
selected_data_ML <- high_corr_metabolites_ML[high_corr_metabolites_ML$spearman_rho > 0, ]
selected_data_ML 
high_corr_metabolites_HL$spearman_rho
high_corr_metabolites_ML$metabolite
high_corr_metabolites_HL$metabolite
merged_metabolites <- unique(c(selected_data_HL$metabolite, selected_data_ML$metabolite))
merged_metabolites
library(openxlsx)
file_path <-"/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-c/1.Data_Assess/all_group/all_sample_data.xlsx"
all_sample_data_metabolite0<- read.xlsx(file_path,1)
head(all_sample_data_metabolite0)
data_mat<- all_sample_data_metabolite0[, c(1,16:24)] 
data_mat[1:5, 1:10]





all_high_corr_metabolites<-data_mat[data_mat$Index%in%merged_metabolites,]
all_high_corr_metabolites
dim(all_high_corr_metabolites)
rownames(all_high_corr_metabolites) <- all_high_corr_metabolites$Index   # 先将Index列设置为行名
all_high_corr_metabolites$Index <- NULL           # 删除Index列

index2compound0<-all_sample_data_metabolite0[,c(1,2)]
index2compound0
head (index2compound0)
# 取旧行名
old_rownames <- rownames(all_high_corr_metabolites)

new_rownames<-index2compound0[index2compound0$Index%in%old_rownames,]
new_rownames$Compounds
# 建立Index到compound的映射（命名向量）
# 用映射替换旧行名，如果没有对应的新名字保持原名
# 赋新行名
rownames(all_high_corr_metabolites) <-new_rownames$Compounds






met_normalized <- t(apply(all_high_corr_metabolites, 1, scale))
colnames(met_normalized)<-colnames(all_high_corr_metabolites) 
met_normalized
Rna_matrix


# 载入必要包
library(pheatmap)

# 假设met_normalized和Rna_matrix行名分别是代谢物和基因名，列名是组名且对应
common_samples <- intersect(colnames(met_normalized), colnames(Rna_matrix))

met_sub <- met_normalized[, common_samples]
rna_sub <- Rna_matrix[, common_samples]

# 初始化相关系数矩阵
cor_mat <- matrix(NA, nrow=nrow(met_sub), ncol=nrow(rna_sub),
                  dimnames = list(rownames(met_sub), rownames(rna_sub)))

# 计算相关系数（皮尔逊）
for(i in 1:nrow(met_sub)) {
  for(j in 1:nrow(rna_sub)) {
    cor_mat[i, j] <- cor(as.numeric(met_sub[i, ]), as.numeric(rna_sub[j, ]), method="pearson")
  }
}
cor_mat

# 筛选高相关对（阈值可调整）
threshold <- 0.80
high_cor_indices <- which(abs(cor_mat) > threshold, arr.ind = TRUE)
dim(high_cor_indices)
#high_cor_indices <- which(cor_mat > threshold, arr.ind = TRUE)


# 构建筛选的相关性子矩阵
high_cor_mat <- cor_mat[unique(high_cor_indices[,1]), unique(high_cor_indices[,2])]

dim(cor_mat)
dim(high_cor_mat)

high_cor_mat
colnames(high_cor_mat)

# 可视化高相关子矩阵（代谢物 vs 基因相关系数热图）


pheatmap(high_cor_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Metabolite-Gene Correlation Heatmap",
         fontsize_row = 12,
         fontsize_col = 8,
         angle_col = 45  # 旋转45度
)



library(igraph)

# 假设cor_mat是代谢物（行）与基因（列）相关系数矩阵
# 相关系数阈值

# 筛选满足阈值的相关关系，获得边表（代谢物、基因、相关性）


edges <- which(abs(cor_mat) > threshold, arr.ind = TRUE)    # 我试试正相关的
#edges <- which(cor_mat> threshold, arr.ind = TRUE)    # 我试试正相关的
edge_list <- data.frame(
  from = rownames(cor_mat)[edges[,1]],
  to = colnames(cor_mat)[edges[,2]],
  weight = cor_mat[edges]
)


write.csv(edge_list,"edge_list.csv")

# 创建节点表，代谢物和基因分别作为不同类型节点
metabolites <- unique(edge_list$from)
genes <- unique(edge_list$to)
nodes <- data.frame(name = c(metabolites, genes),
                    type = c(rep("Metabolite", length(metabolites)), rep("Gene", length(genes))))


# 创建图对象
g <- graph_from_data_frame(d=edge_list, vertices=nodes, directed=FALSE)

# 设定节点颜色和形状
V(g)$color <- ifelse(V(g)$type == "Metabolite", "orange", "skyblue")
V(g)$shape <- ifelse(V(g)$type == "Metabolite", "circle", "square")

# 画图
plot(g,
     vertex.label.cex = 0.7,
     vertex.size = 10,
     edge.width = abs(E(g)$weight)*5,   # 根据相关权重调整边宽
     main = "Metabolite-Gene Correlation Network")

####
library(igraph)
metabolites
# 假设 metabolite_nodes 是你关注的两个代谢物节点名称

metabolite_nodes <- metabolites[c(2,4,5,6,7,8)]
metabolites[c(3,4,5,6,7,8)]
metabolites

all_sample_data_metabolite0

metabolites_k<-all_sample_data_metabolite0[all_sample_data_metabolite0$Compounds%in%metabolites[c(3,4,5,6,7,8)],]

metabolites_k$kegg_map


metabolite_nodes
#metabolite_nodes <- "5,7,8,4'-Tetramethoxyflavone*"
# 找到节点在igraph中的索引
metab_indexes <- which(V(g)$name %in% metabolite_nodes)
metab_indexes

neighbors_list <- adjacent_vertices(g, metab_indexes, mode = "all")
connected_nodes <- unique(unlist(neighbors_list))
connected_nodes
neighbor_names <- V(g)$name[connected_nodes]
# 打印所有邻居名字
print(neighbor_names)

library(ggraph)
library(igraph)
library(tidygraph)
V(g)$name
### 转到cyto
library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g, title = "My Network Metablites", collection = "MyCollection")
# 假设g是igraph对象，节点名字含代谢物标识


V(g)$type <- ifelse(grepl("Cluster", V(g)$name, ignore.case=TRUE), "Gene", "Metabolite")
V(g)$type <- factor(V(g)$type)
V(g)$type
V(g)$shape <- ifelse(grepl("metabolite", V(g)$name, ignore.case=TRUE), "square", "circle")
tg <- as_tbl_graph(g)
set.seed(55)
library(tidygraph)
library(ggraph)

# Assuming 'tg' is your tbl_graph object

# Calculate degree centrality and add as node attribute
tg <- tg %>% 
  mutate(deg = centrality_degree())

# Then use the stored numeric attribute in aes(size=deg)
ggraph(tg, layout = "lgl") +
  geom_edge_link(aes(width = abs(weight)), alpha = 0.3, color = "gray50") +
  geom_node_point(aes(color = type, size = deg, shape = type)) + 
  geom_node_text(aes(label = name), repel = TRUE, size = 5) +
  scale_edge_width(range = c(0.2, 2)) +
  scale_color_manual(values = c("Metabolite" = "tomato", "Gene" = "dodgerblue")) +
  scale_shape_manual(values = c("Metabolite" = 15, "Gene" = 16)) +
  theme_void() +
  ggtitle("Metabolite-Gene Correlation Network")

####条件画图
ggraph(tg, layout = "lgl") +
  geom_edge_link(aes(width = abs(weight)), alpha = 0.3, color = "gray50") +
  geom_node_point(aes(color = type, size = deg, shape = type)) + 
  geom_node_text(aes(label = ifelse(type == "Metabolite", name,
                                   ifelse(deg < 4, name, ""))),
                 repel = TRUE, size = 2) +
  scale_edge_width(range = c(0.2, 2)) +
  scale_color_manual(values = c("Metabolite" = "tomato", "Gene" = "dodgerblue")) +
  scale_shape_manual(values = c("Metabolite" = 15, "Gene" = 16)) +
  theme_void() +
  ggtitle("Metabolite-Gene Correlation Network")



high_cor_mat
genes_107<-colnames(high_cor_mat)

library(openxlsx)
file_path <-"/Users/chaochen/Desktop/Acanthus/MWXS-25-2560-b/05.Quantitation/Unigene.count.annot.xlsx"
all_sample_data_RNA0<- read.xlsx(file_path,1)
head(all_sample_data_RNA0,5)
genes_107_on<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%genes_107,]

dim(genes_107_on)

write.csv(genes_107_on,"genes_107_on.csv")



######################################### 从significant_rna 找到转录因子################################3
## ---- 从significant_rna 找到转录因子----

print(neighbor_names)

length(neighbor_names)


significant_rna
dim(significant_rna)
dim(high_corr_rna)
high_corr_rna$rna <- sub("\\.rho$", "", high_corr_rna$rna)

####我这里做切换为我的neighbor_names。
high_corr_rna_tranfac<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%high_corr_rna$rna,]
head(high_corr_rna_tranfac)
dim(high_corr_rna_tranfac)

# high_corr_rna_tranfac<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%neighbor_names,]
# head(high_corr_rna_tranfac)
# dim(high_corr_rna_tranfac)


library(dplyr)

unique_nonempty <- high_corr_rna_tranfac %>%
  filter(!is.na(`TF-Family`) & `TF-Family` != "" & `TF-Family` != "--") %>%
  distinct(`TF-Family`) %>%
  pull(`TF-Family`)
print(unique_nonempty)
unique_nonempty

length(unique_nonempty)
high_corr_rna_TF<-high_corr_rna_tranfac[high_corr_rna_tranfac$`TF-Family`%in%unique_nonempty,]
head (high_corr_rna_TF)

dim(high_corr_rna_TF)




###RNA_matrix_107 <- genes_107_on[,c(1:10)]
#这里换

RNA_matrix_107 <- genes_107_on[,c(1:10)]
RNA_matrix_107 <-RNA_matrix_107[RNA_matrix_107$ID%in%neighbor_names,]

TF_matrix<-high_corr_rna_TF[,c(1:10)]
names(high_corr_rna_TF)


print(intersect(RNA_matrix_107$ID,TF_matrix$ID))  #这个基因要注意下。 是rna 也是tf




# 假设RNA_matrix_107和TF_matrix已读入，格式如示例所示
# 去ID列，只保留表达count矩阵
rna_mat <- as.matrix(RNA_matrix_107[, -1])


tf_mat <- as.matrix(TF_matrix[, -1])

rownames(rna_mat)<-RNA_matrix_107$ID
rownames(tf_mat)<-TF_matrix$ID



gene_info1 <- data.frame(
  name = genes_107_on$ID,
  gene_name = trimws(gsub("\\[.*?\\]", "", genes_107_on$NR))
)

gene_info2 <- data.frame(
  name = high_corr_rna_TF$ID,
  gene_name = high_corr_rna_TF$`TF-Family`
)

merged_gene_info <- bind_rows(gene_info1, gene_info2)
merged_gene_info



# 行基因的Z-score标准化（对每个基因做归一化）
rna_scaled <- t(scale(t(rna_mat), center = TRUE, scale = TRUE))
tf_scaled <- t(scale(t(tf_mat), center = TRUE, scale = TRUE))

dim(rna_scaled)  # 
dim(tf_scaled)

# 计算Pearson相关矩阵，行是rna基因，列是TF

cor_matrix <- cor(t(rna_scaled), t(tf_scaled), method = "pearson", use = "pairwise.complete.obs")

dim(cor_matrix)
threshold <- 0.95
# 查看筛选后的矩阵
#￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥
# 筛选高相关对（阈值可调整）


high_cor_indices <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)  # 我这里在调整下。

#high_cor_indices <- which(cor_matrix > threshold, arr.ind = TRUE)  # 我这里在调整下。

high_cor_indices
# 构建筛选的相关性子矩阵
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)

# 根据索引提取基因及转录因子名和对应相关系数
high_cor_df <- data.frame(
  RNA_gene = rownames(cor_matrix)[high_cor_indices[,1]],
  TF_gene = colnames(cor_matrix)[high_cor_indices[,2]],
  Correlation = cor_matrix[high_cor_indices]
)
dim(high_cor_df )

high_cor_df
write.csv(high_cor_df,"high_cor_df_tfrna.csv")


head(high_cor_df)

dim(high_cor_df)



# 构建边表
edges <- high_cor_df %>%
  #filter(Correlation > 0) %>%  # 若只保留正相关
  dplyr::select(from = RNA_gene, to = TF_gene, weight = Correlation)

# 如果有节点类型信息，也可以构建节点表
nodes <- data.frame(
  name = unique(c(edges$from, edges$to)),
  type_TF = ifelse(unique(c(edges$from, edges$to)) %in% rownames(rna_scaled), "RNA", "TF")
)
nodes
merged_gene_info

merged_gene_info <- distinct(merged_gene_info, name, .keep_all = TRUE)

anyDuplicated(nodes$name)
dim(nodes)

# 假设gene_info中包含name和gene_name列 为了作图有名字

nodes <- merge(nodes, merged_gene_info, by = "name", all.x = TRUE)

head(nodes)
##  -----node 进一步kegg分析----
nodes_name<-all_sample_data_RNA0[all_sample_data_RNA0$ID%in%nodes$name,]

nodes_name
head(nodes_name)
names(nodes_name)


write.xlsx(nodes_name[, c("ID", "NR","TF-Family")], file = "TableS4.tg_nodes.xlsx", sheetName = "Sheet1", names = FALSE)
write.xlsx(nodes_name[, c("ID", "NR","TF-Family","KEGG")], file = "Table.kegg.xlsx", sheetName = "Sheet1", names = FALSE)

write.xlsx(nodes_name, file = "Table.kegg.xlsx", sheetName = "Sheet1", names = FALSE)


nodes_name[,1:10]
tail(nodes_name[,1:10])
# Extract the first 10 columns from nodes_name
subset_nodes <- nodes_name[, 2:10]
dim(subset_nodes)
subset_nodes <- subset_nodes[complete.cases(subset_nodes), ]
dim(subset_nodes)
dim(subset_nodes)
# Extract the name column (column 1)
rownames(subset_nodes)<- nodes_name[, 1]
subset_nodes
subset_nodes_no_zero <- subset_nodes[rowSums(subset_nodes == 0) == 0, ]
dim(subset_nodes_no_zero)

subset_nodes_no_zero <- subset_nodes_no_zero[complete.cases(subset_nodes_no_zero), ]
dim(subset_nodes_no_zero)
# Calculate row means
row_means <- rowMeans(subset_nodes_no_zero)
row_medians <- apply(subset_nodes, 1, median)
# Find the cutoff value for the bottom 25% (leaving top 75%)
cutoff <- quantile(row_medians, probs = 0.9, na.rm = TRUE)
# 检查 row_medians 中是否含有 NA

# 输出结果
# Select only rows with row means above or equal to cutoff (top 75%)
filtered_subset <- subset(subset_nodes_no_zero, row_medians >= cutoff)
filtered_subset_no_empty <- filtered_subset[complete.cases(filtered_subset), ]
#ids<-rownames(filtered_subset)
ids<-rownames(filtered_subset_no_empty)
length(ids)
cat(paste(ids, collapse = " "), "\n")

##-----------网页的分析------------

# 安装并加载ggplot2包（如果未安装）
# install.packages("ggplot2")
library(ggplot2)

# 构建数据框
data <- data.frame(
  PathwayID = c("map01110", "map01100", "map01120", "map01200", "map01230",
                "map01240", "map00010", "map01210", "map00500", "map04075",
                "map00270", "map00564", "map00620", "map00020", "map00900",
                "map00260", "map00630", "map00030", "map00710", "map00940"),
  PathwayName = c("Biosynthesis of secondary metabolites", "Metabolic pathways", 
                  "Microbial metabolism in diverse environments", "Carbon metabolism", 
                  "Biosynthesis of amino acids", "Biosynthesis of cofactors", 
                  "Glycolysis / Gluconeogenesis", "2-Oxocarboxylic acid metabolism", 
                  "Starch and sucrose metabolism", "Plant hormone signal transduction", 
                  "Cysteine and methionine metabolism", "Glycerophospholipid metabolism", 
                  "Pyruvate metabolism", "Citrate cycle (TCA cycle)", "Terpenoid backbone biosynthesis", 
                  "Glycine, serine and threonine metabolism", "Glyoxylate and dicarboxylate metabolism", 
                  "Pentose phosphate pathway", "Carbon fixation by Calvin cycle", "Phenylpropanoid biosynthesis"),
  Count = c(294, 266, 69, 56, 55, 32, 28, 25, 20, 19, 18, 16, 16, 16, 16, 16, 15, 14, 14, 14)
)

# 绘制柱状图
ggplot(data, aes(x = reorder(PathwayName, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # 横向柱状图
  labs(title = "KEGG Pathway Gene Count Distribution",
     x = "Pathway Name",
     y = "Gene Count")+
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10)) +
  scale_fill_gradient(low = "skyblue", high = "steelblue")




head(nodes_name)
nodes_name
#Peroxisome (7)
map_node<-c(
"K00031",
"K00232",
"K00830",
"K00869",
"K03781",
"K07513",
"K11517"
)

#Glucagon signaling pathway (8)
map_node <- list(
  "K00016",
  "K00688",
  "K00850",
  "K00873",
  "K01834",
  "K03841",
  "K04498",
  "K07198"
)

#Glutathione metabolism (5)
map_node <- list(
  "K00031",
  "K00033",
  "K00036",
  "K01581",
  "K11188"
)
#Flavonoid biosynthesis (5)
map_node <- list(
  "K00588",
  "K01859",
  "K09754",
  "K13065",
  "K23179"
)


pattern <- paste(map_node, collapse = "|")
map_node_id <- nodes_name[grepl(pattern, nodes_name$KEGG), ]
map_node_id 
dim(map_node_id)
idsmap<-map_node_id$ID
cat(paste(idsmap, collapse = " "), "\n")

map_indexes <- which(V(g)$name %in% map_node_id$ID)
#######暂停下。 

map_node_neighbors_list <- adjacent_vertices(g, map_indexes, mode = "all")
connected_nodes <- unique(unlist(map_node_neighbors_list))
neighbor_names <- V(g)$name[connected_nodes]
# 打印所有邻居名字
print(neighbor_names)

# 创建igraph对象
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
g
# #write_graph(g, file = "network.graphml", format = "graphml")
library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g, title = "My Network rna_tf_616all0.8", collection = "MyCollection")
#上传节点信息
nodes_df <- data.frame(name = V(g)$name, gene_name = V(g)$gene_name,  TF = V(g)$type_TF)
loadTableData(nodes_df, data.key.column = "name", table = "node")

length(V(g)$name)       # 应该是659
length(V(g)$gene_name)  # 如果为0或不等于上面长度，说明缺少数据
length(V(g)$type_TF)    # 同上



# 转为tidygraph对象
V(g)$degree <- degree(g)
tg <- as_tbl_graph(g)
tg <- tg %>% mutate(degree = centrality_degree())




library(ggraph)
options(ggrepel.max.overlaps = Inf)
set.seed(42)
ggraph(tg, layout = "lgl") +
  geom_edge_link(aes(width = weight), alpha = 0.005) +
  geom_node_point(aes(color = type_TF, shape = type_TF, size = degree)) +  # size by degree
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  #geom_node_text(aes(label = degree), vjust = 1.5, size = 2, color = "red") +  # add degree labels
  scale_color_manual(values = c(RNA = "dodgerblue", TF = "tomato")) +
  scale_shape_manual(values = c(RNA = 16, TF = 15)) +
  theme_void() +
  ggtitle("High Correlation RNA - TF Network")

E(tg)$weight
tg
tg_nodes <- tg %>% activate(nodes) %>% as_tibble()
tg_nodes
dim(tg_nodes)





################取个核心网络
library(dplyr)
library(tidygraph)
# 假设 tg 是一个 tbl_graph 网络对象
# 过滤度大于某阈值的节点

# 基于过滤的节点和边重新构建子图
core_tg <- tg %>%
  activate(nodes) %>%
  filter(degree > 4) %>%
  activate(edges) %>%
  filter(weight > 0.95) %>%
  activate(nodes) %>% 
  filter(degree > 0)  # 删除孤立节点

core_nodes <- core_tg %>% activate(nodes) %>% as_tibble()
print(core_nodes)
core_nodes
dim(core_nodes)



dim(core_nodes)
core_net<-diff_geneall[diff_geneall$ID%in%core_nodes$name,]

dim(core_net)
write.csv(core_net,"core_net.csv")

# 绘图
ggraph(core_tg, layout = "lgl") +  # kk  fr graphopt  fr lgl
  geom_edge_link(aes(width = weight), alpha = 0.05) +
  geom_node_point(aes(color = type_TF, shape = type_TF, size = degree)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(RNA = "dodgerblue", TF = "tomato")) +
  scale_shape_manual(values = c(RNA = 16, TF = 17)) +
  theme_void() +
  ggtitle("Core High Correlation RNA - TF Network")

core_tg
metabolites
metabolites[c(5,6,7,8)]
all_sample_data_metabolite0
metabolites_fenlei<-all_sample_data_metabolite0[all_sample_data_metabolite0$Compounds%in%metabolites,]
table(metabolites_fenlei$物质二级分类)
metabolites_fenlei[,c("Compounds","物质二级分类")]
metabolites_fenlei

# library(RCy3)
# cytoscapePing()
# createNetworkFromIgraph(g, title = "My Network RNA", collection = "MyCollection")


#https://www.genome.jp/kegg/ko.html  这个可以做kegg