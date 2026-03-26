# =========================================================================
# 脚本名称: 04_ML_Diagnostic_Panel_Auto.R
# 目标: 自动转置矩阵、智能匹配代谢物名、构建多组学联合诊断 Panel
# =========================================================================

library(randomForest)
library(caret)
library(pROC)
library(dplyr)

print("📂 Step 1: 读取正负离子矩阵...")
pos_data <- read.delim("st002787_positive_clean.txt", stringsAsFactors = FALSE, check.names = FALSE)
neg_data <- read.delim("st002787_negative_clean.txt", stringsAsFactors = FALSE, check.names = FALSE)

print("🔍 Step 2: 智能模糊匹配核心代谢物的精确名称...")
# 在 Metabolite_name 这一列中寻找目标
ido1_exact <- grep("Kynurenic", pos_data$Metabolite_name, ignore.case = TRUE, value = TRUE)[1]
asah1_exact <- grep("Hex3Cer", pos_data$Metabolite_name, ignore.case = TRUE, value = TRUE)[1]
# 胆汁酸防缩写策略：寻找全称或常见缩写 (UDCA/DCA/CDCA等)
fxr_exact <- grep("Ursodeoxycholic|UDCA|Deoxycholic|DCA", neg_data$Metabolite_name, ignore.case = TRUE, value = TRUE)[1]

print(paste("找到色氨酸轴代表物:", ido1_exact))
print(paste("找到神经酰胺代表物:", asah1_exact))
print(paste("找到胆汁酸轴代表物:", fxr_exact))

if(is.na(ido1_exact) | is.na(asah1_exact) | is.na(fxr_exact)) {
  stop("❌ 严重错误：未找到上述某个代谢物，请检查数据表！")
}

print("🔄 Step 3: 转置矩阵，将样本作为行，代谢物作为列...")
# 提取所需的行，并转置
pos_sub <- pos_data[pos_data$Metabolite_name %in% c(ido1_exact, asah1_exact), ]
rownames(pos_sub) <- pos_sub$Metabolite_name
pos_t <- as.data.frame(t(pos_sub[, -c(1,2)])) # 剔除前两列非数值注释
pos_t$SampleID <- rownames(pos_t)

neg_sub <- neg_data[neg_data$Metabolite_name == fxr_exact, ]
rownames(neg_sub) <- neg_sub$Metabolite_name
neg_t <- as.data.frame(t(neg_sub[, -c(1,2)]))
neg_t$SampleID <- rownames(neg_t)

# 合并正负离子数据
merged_data <- merge(pos_t, neg_t, by = "SampleID")

print("🏷️ Step 4: 从样本ID中自动提取分组信息 (CRA vs CRC)...")
merged_data$Group <- ifelse(grepl("^CRA", merged_data$SampleID), "CRA",
                            ifelse(grepl("^CRC", merged_data$SampleID), "CRC",
                                   ifelse(grepl("^HC", merged_data$SampleID), "HC", NA)))

# 过滤数据，仅保留 CRA 和 CRC，转换为因子，必须将数值列转化为 numeric
ml_data <- merged_data %>% filter(Group %in% c("CRA", "CRC"))
ml_data$Group <- factor(ml_data$Group, levels = c("CRA", "CRC"))

# 确保所有特征列都是数值型
ml_data[[ido1_exact]] <- as.numeric(ml_data[[ido1_exact]])
ml_data[[asah1_exact]] <- as.numeric(ml_data[[asah1_exact]])
ml_data[[fxr_exact]] <- as.numeric(ml_data[[fxr_exact]])

print("🧹 Step 5: 清理列名以适应 R 语言模型公式 (消除括号和空格)...")
var_ido1 <- make.names(ido1_exact)
var_asah1 <- make.names(asah1_exact)
var_fxr <- make.names(fxr_exact)
colnames(ml_data) <- make.names(colnames(ml_data))

print("🤖 Step 6: 运行 5折交叉验证 Random Forest...")
set.seed(123)
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
formula_str <- as.formula(paste("Group ~", var_ido1, "+", var_asah1, "+", var_fxr))

rf_model <- train(formula_str, data = ml_data, 
                  method = "rf", metric = "ROC",
                  trControl = train_control, importance = TRUE)

print("🌟 随机森林各标志物重要性:")
print(varImp(rf_model))

print("🧮 Step 7: 构建 Logistic 联合诊断模型并计算 AUC...")
glm_model <- glm(formula_str, data = ml_data, family = binomial)
ml_data$Panel_Prob <- predict(glm_model, type = "response")

roc_ido1 <- roc(ml_data$Group, ml_data[[var_ido1]], quiet = TRUE)
roc_panel <- roc(ml_data$Group, ml_data$Panel_Prob, quiet = TRUE)

auc_ido1 <- ci.auc(roc_ido1)
auc_panel <- ci.auc(roc_panel)

print("================ 【核心结果汇报】 ================")
print(paste(ido1_exact, "单标 AUC:", round(auc_ido1[2], 3), 
            " (95% CI:", round(auc_ido1[1], 3), "-", round(auc_ido1[3], 3), ")"))
print(paste("3轴联合 Panel AUC:", round(auc_panel[2], 3), 
            " (95% CI:", round(auc_panel[1], 3), "-", round(auc_panel[3], 3), ")"))

test_p <- roc.test(roc_panel, roc_ido1, method = "delong")$p.value
print(paste("DeLong 检验 P-value:", signif(test_p, 3)))
print("==================================================")

print("🎨 Step 8: 导出高颜值 Fig 7 ROC 曲线...")
pdf("Fig7_ML_Diagnostic_ROC.pdf", width = 6, height = 6)
plot(roc_panel, col = "#D55E00", lwd = 3, main = "Diagnostic Performance (CRA vs CRC)",
     print.auc = FALSE, legacy.axes = TRUE)
plot(roc_ido1, col = "#0072B2", lwd = 2, lty = 2, add = TRUE)

legend_panel <- sprintf("3-Axis Panel (AUC=%.2f, 95%% CI: %.2f-%.2f)", auc_panel[2], auc_panel[1], auc_panel[3])
legend_ido1 <- sprintf("%s (AUC=%.2f)", ido1_exact, auc_ido1[2])
legend_p <- sprintf("DeLong P = %.3g", test_p)

legend("bottomright", 
       legend = c(legend_panel, legend_ido1, "", legend_p),
       col = c("#D55E00", "#0072B2", NA, NA),
       lwd = c(3, 2, NA, NA), lty = c(1, 2, NA, NA),
       bty = "n", cex = 0.9)
dev.off()
print("🎉 恭喜！Fig 7 PDF 文件已生成。")
