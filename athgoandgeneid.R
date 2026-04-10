# 载入必要的包
library(httr)

# 目标URL
url <- "https://rest.kegg.jp/link/ath/ko"

# 发送GET请求并获取内容
response <- GET(url)

# 检查请求是否成功
if (status_code(response) == 200) {
  # 获取文本内容
  content_text <- content(response, "text", encoding = "UTF-8")
  
  # 将内容按行拆分为向量
  data_lines <- strsplit(content_text, "\n")[[1]]
  
  # 将数据转换为数据框，以制表符分隔
  data_df <- do.call(rbind, strsplit(data_lines, "\t"))
  data_df <- as.data.frame(data_df, stringsAsFactors = FALSE)
  colnames(data_df) <- c("KO", "ATH")
  
  # 展示前几行
  print(head(data_df))
} else {
  cat("请求失败，状态码：", status_code(response), "\n")
}