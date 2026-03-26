# =========================================================================
# 脚本名称: 03a_TCGA_Expression_Fetcher.R
# 目标: 严谨获取原发肿瘤 ASAH1 表达量，生成 GitHub 标准化数据矩阵
# =========================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(stringr)

# 清理可能损坏的旧缓存 
if(dir.exists("GDCdata")) {
  unlink("GDCdata", recursive = TRUE)
  print("🗑️ 已清理旧的 GDCdata 缓存，准备重新干净下载。")
}

print("🚀 正在向 GDC 请求原发肿瘤 (Primary Tumor) RNA-Seq 数据...")
query_expr <- GDCquery(
  project = c("TCGA-COAD", "TCGA-READ"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

# 增加防弹下载参数 (分块下载，强制使用稳定API)
print("⬇️ 正在分块下载数据，已开启高稳定性模式，请稍候...")
tryCatch({
  GDCdownload(query_expr, method = "api", files.per.chunk = 20)
}, error = function(e) {
  print("⚠️ 捕捉到网络波动，尝试第二次重连下载...")
  GDCdownload(query_expr, method = "api", files.per.chunk = 10)
})

print("⚙️ 正在组装表达矩阵...")
expr_data_se <- GDCprepare(query_expr)

print("🎯 正在提取 ASAH1 基因 (ENSG00000104763) 的 TPM 表达量...")
tpm_matrix <- assay(expr_data_se, "tpm_unstrand")
gene_annotation <- rowData(expr_data_se)

# 双重保险定位基因
asah1_index <- which(gene_annotation$gene_name == "ASAH1")
if(length(asah1_index) == 0) {
  asah1_index <- which(str_detect(gene_annotation$gene_id, "ENSG00000104763"))
}

asah1_expression_values <- tpm_matrix[asah1_index, ]
patient_ids <- str_sub(names(asah1_expression_values), 1, 12)

# 构建数据框
asah1_df <- data.frame(
  Patient_ID = patient_ids,
  ASAH1_TPM = as.numeric(asah1_expression_values),
  stringsAsFactors = FALSE
)

# 剔除重复取样，保留单个患者最大表达量
asah1_df <- asah1_df %>%
  group_by(Patient_ID) %>%
  slice_max(ASAH1_TPM, n = 1, with_ties = FALSE) %>%
  ungroup()

# 按中位数严格二分类
median_expr <- median(asah1_df$ASAH1_TPM, na.rm = TRUE)
asah1_df <- asah1_df %>%
  mutate(
    ASAH1_Group = ifelse(ASAH1_TPM >= median_expr, "High", "Low"),
    ASAH1_Group = factor(ASAH1_Group, levels = c("Low", "High")) 
  )

# =========================================================================
github_filename <- "InputData_03_TCGA_CRC_ASAH1_TPM_Stratified.csv"

write.csv(asah1_df, github_filename, row.names = FALSE)
print(paste("✅ 成功！高质量数据已保存为:", github_filename))

# =========================================================================
# 脚本名称: 03b_Clinical_Cox_Survival.R
# 目标: 下载TCGA临床数据，合并ASAH1表达量，执行多因素Cox回归，绘制森林图
# =========================================================================

# 1. 加载必要的包
library(TCGAbiolinks)
library(survival)
library(survminer)
library(dplyr)
library(stringr)

print("🚀 Step 1: 正在从 GDC 下载最新的 COAD/READ 临床表型数据...")
clin_coad <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
clin_read <- GDCquery_clinic(project = "TCGA-READ", type = "clinical")
clin_crc <- bind_rows(clin_coad, clin_read)

print("🧹 Step 2: 正在执行极度严苛的临床数据清洗与变量重构...")
clin_clean <- clin_crc %>%
  select(submitter_id, vital_status, days_to_death, days_to_last_follow_up, 
         age_at_index, gender, ajcc_pathologic_stage) %>%
  mutate(
    OS.time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    OS.status = ifelse(vital_status == "Dead", 1, 0),
    Age = as.numeric(age_at_index),
    Stage_Raw = as.character(ajcc_pathologic_stage),
    Stage = case_when(
      str_detect(Stage_Raw, "Stage I$|Stage II$|Stage IA|Stage IIA|Stage IIB|Stage IIC") ~ "Early (I-II)",
      str_detect(Stage_Raw, "Stage III|Stage IV|Stage IIIA|Stage IIIB|Stage IIIC|Stage IVA|Stage IVB") ~ "Advanced (III-IV)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(OS.time) & OS.time > 0) %>%
  filter(!is.na(Age) & !is.na(Stage)) %>%
  mutate(
    Gender = factor(gender, levels = c("female", "male")),
    Stage = factor(Stage, levels = c("Early (I-II)", "Advanced (III-IV)"))
  )
clin_clean$submitter_id <- str_sub(clin_clean$submitter_id, 1, 12)

print("🔗 Step 3: 读取表达量数据并进行合并...")
# 读取上一步生成的干净表达量表
asah1_expr <- read.csv("InputData_03_TCGA_CRC_ASAH1_TPM_Stratified.csv")
asah1_expr$ASAH1_Group <- factor(asah1_expr$ASAH1_Group, levels = c("Low", "High"))

final_data <- inner_join(clin_clean, asah1_expr, by = c("submitter_id" = "Patient_ID"))

# 为 GitHub 保存一份最终的、干净的用于生存分析的表格
write.csv(final_data, "InputData_04_TCGA_CRC_Clinical_Survival.csv", row.names = FALSE)
print("✅ 临床清洗完毕！GitHub 数据已保存为: InputData_04_TCGA_CRC_Clinical_Survival.csv")

print("📊 Step 4: 正在构建多因素 Cox 比例风险模型 (Multivariate Cox)...")
# 模型公式：生存 ~ ASAH1表达分组 + 年龄 + 性别 + TNM分期
res.cox <- coxph(Surv(OS.time, OS.status) ~ ASAH1_Group + Age + Gender + Stage, data = final_data)

# 【核心任务】：请把这里打印出来的 summary 结果复制发给我！
print(summary(res.cox))

print("🎨 Step 5: 正在绘制符合顶刊审美的森林图 (Fig 6D)...")
p_forest <- ggforest(res.cox, 
                     data = final_data, 
                     main = "Multivariate Cox Regression in TCGA-CRC", 
                     fontsize = 1.0, 
                     refLabel = "Reference")

pdf("Fig6D_Multivariate_Cox_Forest.pdf", width = 8, height = 5)
print(p_forest)
dev.off()
print("🎉 恭喜！Fig 6D 森林图已导出至当前文件夹。")
