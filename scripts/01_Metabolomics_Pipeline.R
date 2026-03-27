# ==========================================
# 结直肠癌 (ST002787) 终极双离子模式全自动分析流水线 
# (顶刊级美颜版 - 解决重叠 & 统一配色)
# ==========================================

# 1. 加载所有依赖包
suppressMessages({
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(impute)
  library(ggplot2)
  library(scales)
  library(ropls)        # 代谢组金标准：带置换检验的PLS-DA
  library(pROC)         # 计算诊断效能 AUC 及 CI
  library(ggVennDiagram)# 绘制韦恩图
  library(ggrepel)      # 火山图标签防重叠
  library(pheatmap)     # 画热图
})

# ==========================================
# 【美颜配置区】：定义全局统一的高级配色方案
# ==========================================
col_hc  <- "#1B9E77"   # 绿色 (Healthy Control)
col_cra <- "#D95F02"   # 橙色 (Colorectal Adenoma)
col_crc <- "#E7298A"   # 红色 (Colorectal Cancer)
col_ns  <- "#CCCCCC"   # 灰色 (Not Significant)

# ==========================================
# 第一部分：精确清洗分组信息（从云端拉取，只需运行一次）
# ==========================================
cat("📥 正在从数据库拉取临床分组信息...\n")
factors_url <- "https://www.metabolomicsworkbench.org/rest/study/study_id/ST002787/factors"
metadata_json <- fromJSON(factors_url)

sample_info <- bind_rows(metadata_json) %>%
  mutate(
    SampleID = trimws(local_sample_id),
    Group = str_extract(factors, "(?<=Group type:)[^|]+") %>% trimws()
  ) %>%
  dplyr::select(SampleID, Group) %>%
  mutate(Group = factor(Group, levels = c("Heathy control", "Colorectal adenoma", "Colorectal cancer")))

cat("✅ 分组信息准备完毕！\n")

# ==========================================
# 第二部分：定义“自动化分析工厂”函数（包含所有清洗、统计和出图）
# ==========================================
run_metabolomics_pipeline <- function(data_path, mode_prefix, sample_info) {
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("🚀 机器启动：正在全自动处理 [%s] 模式数据...\n", mode_prefix))
  cat(sprintf("======================================================\n"))
  
  # 1. 读取与对齐
  local_data <- read.delim(data_path, sep = "\t", check.names = FALSE)
  rownames(local_data) <- make.unique(as.character(local_data[[1]]))
  local_data <- local_data[, -1] 
  if("RefMet_name" %in% colnames(local_data)) local_data$RefMet_name <- NULL
  
  data_matrix <- t(local_data) %>% as.data.frame() %>% rownames_to_column("SampleID") %>% mutate(SampleID = trimws(SampleID))
  analysis_dataset <- inner_join(sample_info, data_matrix, by = "SampleID")
  
  # 2. 严谨的数据预处理 (应对审稿人致命伤)
  expr_mat <- analysis_dataset %>% dplyr::select(-SampleID, -Group) %>% as.matrix()
  group_factor <- analysis_dataset$Group
  
  # (1) 0值替换为NA，并剔除缺失率 > 50% 的代谢物
  expr_mat[expr_mat == 0] <- NA
  missing_rate <- colMeans(is.na(expr_mat))
  expr_mat <- expr_mat[, missing_rate < 0.5] 
  
  # (2) TIC 归一化 (总峰面积校正，消除基质效应)
  tic <- rowSums(expr_mat, na.rm = TRUE)
  expr_mat_tic <- sweep(expr_mat, 1, tic, "/") * median(tic, na.rm = TRUE)
  
  # (3) KNN 填充缺失值与 Pareto 缩放
  expr_imputed <- t(impute.knn(t(expr_mat_tic))$data)
  expr_pareto <- scale(expr_imputed, center = TRUE, scale = sqrt(apply(expr_imputed, 2, sd)))
  rownames(expr_pareto) <- analysis_dataset$SampleID
  
  # 3. 输出图1：三组全局 PCA 轨迹图 (采用统一配色)
  pca_result <- prcomp(expr_pareto, scale. = FALSE)
  pca_df <- as.data.frame(pca_result$x[, 1:2]) %>% mutate(Group = group_factor)
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3.5, alpha = 0.8) + 
    stat_ellipse(level = 0.95, linewidth = 1.2, linetype = "dashed") + 
    scale_color_manual(values = c("Heathy control" = col_hc, 
                                  "Colorectal adenoma" = col_cra, 
                                  "Colorectal cancer" = col_crc)) +
    theme_bw(base_size = 15) + 
    labs(title = paste(mode_prefix, "- Global PCA Trajectory"))
  
  ggsave(paste0(mode_prefix, "_Figure_1_PCA.pdf"), p_pca, width = 8, height = 6, dpi = 300)
  
  # ==============================================
  # 内部函数：专门用于执行两两对比 (非参数检验 + 置换检验PLS-DA)
  # ==============================================
  run_comparison <- function(g1, g2, comp_name) {
    idx <- group_factor %in% c(g1, g2)
    X_par <- expr_pareto[idx, ]
    X_raw <- expr_imputed[idx, ]
    Y_fac <- factor(group_factor[idx], levels = c(g1, g2))
    
    # Wilcoxon 非参数检验计算 P 值与 FDR
    pvals <- apply(X_raw, 2, function(x) wilcox.test(x ~ Y_fac, exact=FALSE)$p.value)
    fdr <- p.adjust(pvals, method = "BH")
    
    # 计算 Fold Change (g2 / g1)
    log2fc <- apply(X_raw, 2, function(x) {
      m1 <- mean(x[Y_fac == g1], na.rm = TRUE)
      m2 <- mean(x[Y_fac == g2], na.rm = TRUE)
      log2((m2 + 1e-5) / (m1 + 1e-5))
    })
    
    # 带有 200次置换检验的 PLS-DA
    cat(sprintf("正在计算 [%s] 的 PLS-DA 及置换检验...\n", comp_name))
    pdf(paste0(mode_prefix, "_Figure_S1_PLSDA_Permutation_", comp_name, ".pdf"), width=6, height=5)
    plsda_mod <- ropls::opls(X_par, Y_fac, predI = 2, permI = 200)
    plot(plsda_mod, typeVc = "permutation")
    dev.off()
    
    # 提取 VIP 值
    vip_values <- getVipVn(plsda_mod)
    vip_df <- data.frame(Metabolite = names(vip_values), VIP = as.numeric(vip_values))
    
    # 计算 ROC 的 AUC 值及其 95% 置信区间
    cat(sprintf("正在计算 [%s] 的 AUC 及 95%% CI...\n", comp_name))
    roc_results <- lapply(as.data.frame(X_raw), function(x) {
      r <- roc(Y_fac, x, quiet = TRUE)
      ci_val <- ci.auc(r)
      return(c(AUC = as.numeric(ci_val[2]), 
               CI_lower = as.numeric(ci_val[1]), 
               CI_upper = as.numeric(ci_val[3])))
    })
    roc_mat <- do.call(rbind, roc_results)
    
    # 整理结果并合并
    res <- data.frame(
      Metabolite = colnames(X_raw),
      Log2FC = log2fc,
      P_value = pvals,
      FDR = fdr,
      AUC = roc_mat[, "AUC"],
      AUC_95CI_Lower = roc_mat[, "CI_lower"],
      AUC_95CI_Upper = roc_mat[, "CI_upper"]
    ) %>% left_join(vip_df, by = "Metabolite") %>%
      mutate(
        Status = case_when(FDR < 0.05 & VIP > 1 & Log2FC > 0.58 ~ "Up",
                           FDR < 0.05 & VIP > 1 & Log2FC < -0.58 ~ "Down",
                           TRUE ~ "Not Sig")
      ) %>% arrange(desc(VIP))
    
    write.csv(res, paste0(mode_prefix, "_Table_", comp_name, "_Full_Results.csv"), row.names = FALSE)
    return(list(res = res, X_par = X_par, Y_fac = Y_fac))
  }
  
  # 4. 执行两组关键对比
  cat("\n---> 分析 1: 健康 vs 腺瘤 (HC vs CRA) ...\n")
  comp_hc_cra <- run_comparison("Heathy control", "Colorectal adenoma", "HC_vs_CRA")
  
  cat("---> 分析 2: 腺瘤 vs 肠癌 (CRA vs CRC) ...\n")
  comp_cra_crc <- run_comparison("Colorectal adenoma", "Colorectal cancer", "CRA_vs_CRC")
  
  # 5. 输出图2：韦恩图 (配色过渡到 CRC 的红色)
  sig_cra <- comp_hc_cra$res %>% filter(Status != "Not Sig") %>% pull(Metabolite)
  sig_crc <- comp_cra_crc$res %>% filter(Status != "Not Sig") %>% pull(Metabolite)
  
  venn_list <- list("Early Adenoma Phase" = sig_cra, "Malignant Transformation (CRC)" = sig_crc)
  p_venn <- ggVennDiagram(venn_list) + 
    scale_fill_gradient(low = "white", high = col_crc) +
    labs(title = paste(mode_prefix, "- Phase Specific Alterations")) + 
    theme(legend.position = "none")
  ggsave(paste0(mode_prefix, "_Figure_2_Venn.pdf"), p_venn, width = 6, height = 5, dpi = 300)
  
  # ==============================================
  # 重点关注恶性转化阶段 (CRA vs CRC)，绘制核心三件套
  # ==============================================
  res_crc <- comp_cra_crc$res
  sig_crc_full <- res_crc %>% filter(Status != "Not Sig")
  
  if(nrow(sig_crc_full) > 0) {
    # 6. 输出图3：火山图 (终极美颜版：解决重叠，统一配色)
    p_volcano <- ggplot(res_crc, aes(x = Log2FC, y = -log10(FDR), color = Status)) +
      geom_point(alpha = 0.8, size = 2.5) +
      scale_color_manual(values = c("Up" = col_crc, "Down" = col_cra, "Not Sig" = col_ns)) +
      geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", alpha = 0.5, color="grey30") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5, color="grey30") +
      geom_text_repel(
        data = head(sig_crc_full, 20),  # 提取Top 20的核心代谢物加标签
        aes(label = Metabolite), 
        size = 3.5, 
        fontface = "bold",
        box.padding = 0.8,              # 增加文本框周围的排斥力
        point.padding = 0.4,            # 增加文本与点之间的距离
        force = 5,                      # 强力推开重叠标签
        segment.color = "grey50",       # 引导线颜色
        segment.size = 0.5,             # 引导线粗细
        min.segment.length = 0,         # 强制画出所有引导线
        max.overlaps = Inf              # 核心：无限重叠也强制全部显示
      ) +
      theme_classic(base_size = 15) +
      labs(title = paste(mode_prefix, "- Volcano: CRA vs CRC"), x = "Log2(FC)", y = "-Log10(FDR)") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    ggsave(paste0(mode_prefix, "_Figure_3_Volcano.pdf"), p_volcano, width = 8, height = 7, dpi = 300)
    
    # 7. 输出图4：Top 15 VIP 排位图
    top_vip <- head(sig_crc_full, 15) %>% arrange(VIP) %>% mutate(Metabolite = factor(Metabolite, levels = Metabolite))
    p_vip <- ggplot(top_vip, aes(x = Metabolite, y = VIP, color = Status)) +
      geom_segment(aes(x = Metabolite, xend = Metabolite, y = 0, yend = VIP), linewidth = 1.2) +
      geom_point(size = 4) +
      scale_color_manual(values = c("Up" = col_crc, "Down" = col_cra)) +
      coord_flip() + theme_bw(base_size = 14) +
      labs(title = paste(mode_prefix, "- Top VIP Drivers"), x = "", y = "VIP Score")
    ggsave(paste0(mode_prefix, "_Figure_4_VIP.pdf"), p_vip, width = 7, height = 6, dpi = 300)
    
    # 8. 输出图5：Top 25 差异代谢物聚类热图
    if(nrow(sig_crc_full) >= 2) {
      top_metabolites <- head(sig_crc_full$Metabolite, 25)
      heatmap_data <- t(comp_cra_crc$X_par[, top_metabolites]) 
      annotation_col <- data.frame(Group = comp_cra_crc$Y_fac)
      rownames(annotation_col) <- rownames(comp_cra_crc$X_par)
      
      # 热图注释条的颜色严格对应分组全局颜色
      ann_colors = list(Group = c("Colorectal adenoma" = col_cra, "Colorectal cancer" = col_crc))
      
      pheatmap(heatmap_data, 
               annotation_col = annotation_col, 
               annotation_colors = ann_colors,
               scale = "row", cluster_cols = TRUE, show_colnames = FALSE, fontsize_row = 8,
               color = colorRampPalette(c("#4DBBD5", "white", col_crc))(50), # 热图表达量色阶
               main = paste(mode_prefix, "- Metabolic Reprogramming (CRA vs CRC)"),
               filename = paste0(mode_prefix, "_Figure_5_Heatmap.pdf"), width = 8, height = 7)
    }
  } else {
    cat(sprintf("⚠️ [%s] 模式在 CRA vs CRC 对比中未找到符合严格标准的差异代谢物。\n", mode_prefix))
  }
  
  cat(sprintf("🎉 [%s] 模式运行完毕！所有图表和CSV已经保存在你的文件夹中了。\n", mode_prefix))
}

# ==========================================
# 第三部分：把数据放进机器里执行！（请务必修改这里的路径）
# ==========================================
# ⚠️ 注意：请把下面两行的路径，换成你电脑上真实的 txt 文件路径
neg_path <- "/Users/bing/ST002787/st002787_negative_clean.txt"
pos_path <- "/Users/bing/ST002787/st002787_positive_clean.txt"

# 开始处理负离子模式！
run_metabolomics_pipeline(data_path = neg_path, mode_prefix = "NEG", sample_info = sample_info)

# 开始处理正离子模式！
run_metabolomics_pipeline(data_path = pos_path, mode_prefix = "POS", sample_info = sample_info)
