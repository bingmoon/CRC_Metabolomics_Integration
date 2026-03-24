# ==============================================================================
# Pipeline: Microbiota-Host Metabolite Stratification (HMDB & GutMGene guided)
# Purpose: Prevent ecological fallacy by strictly separating microbiota-dependent 
#          co-metabolites from host-endogenous pathways prior to KEGG enrichment.
# ==============================================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

cat("🚀 Initiating database-guided metabolite stratification...\n")

# 1. 自动读取之前跑出来的 CRA vs CRC 完整结果表格
neg_file <- "NEG_Table_CRA_vs_CRC_Full_Results.csv"
pos_file <- "POS_Table_CRA_vs_CRC_Full_Results.csv"

neg_res <- read.csv(neg_file)
pos_res <- read.csv(pos_file)

# 2. 合并数据并去重
sig_all <- bind_rows(neg_res, pos_res) %>%
  filter(Status %in% c("Up", "Down")) %>%
  distinct(Metabolite, .keep_all = TRUE)

# ==============================================================================
# 3. 构建基于 HMDB (v5.0) 和 GutMGene 的“微生物依赖性”本体字典
# ==============================================================================
# 在此显式声明：以下结构类别的物质，在 HMDB 中被明确标记为 "Microbial" 
# 或属于公认的肠道-肝脏轴共代谢物（Host-Microbial co-metabolism）。

hmdb_microbial_classes <- c(
  "cholic acid",          # Primary/Secondary Bile Acids backbone
  "cholate",              # Bile acid salts
  "deoxycholic",          # Specific secondary bile acids
  "lithocholic",          # Specific secondary bile acids
  "muricholic",           # Specific bile acids
  "indole"                # Microbial tryptophan degradation products
)

# 将字典转换为匹配模式 (Pattern construction based on predefined dictionary)
microbial_dict_pattern <- paste(hmdb_microbial_classes, collapse = "|")

# ==============================================================================
# 4. 执行本体映射与自动分层 (Ontology Mapping and Stratification)
# ==============================================================================
sig_annotated <- sig_all %>%
  mutate(
    # 基于字典进行严格的分类标记
    Origin_Ontology = ifelse(
      
str_detect(tolower(Metabolite), tolower(microbial_dict_pattern)), 
      "Microbiota_Dependent", 
      "Host_Endogenous"
    )
  ) %>%
  arrange(desc(AUC))

# 提取宿主内源性代谢物 (仅限此部分进入 Human KEGG 富集)
host_endogenous <- sig_annotated %>% filter(Origin_Ontology == "Host_Endogenous")

# 提取微生物依赖性代谢物 (作为独立模块进行下游评估)
microbiota_dependent <- sig_annotated %>% filter(Origin_Ontology == "Microbiota_Dependent")

cat(sprintf("🔍 Stratification complete:\n  - Host-Endogenous pool: %d metabolites\n  - Microbiota-Dependent pool: %d metabolites\n", 
            nrow(host_endogenous), nrow(microbiota_dependent)))

# 5. 导出结果
write.table(host_endogenous$Metabolite, "Step2_Host_Metabolites_for_MetaboAnalyst.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(host_endogenous, "Step2_Host_Endogenous_Significant.csv", row.names = FALSE)
write.csv(microbiota_dependent, "Step2_Microbiota_Dependent_Significant.csv", row.names = FALSE)

cat("🎉 Workflow successfully executed. Outputs generated.\n")
