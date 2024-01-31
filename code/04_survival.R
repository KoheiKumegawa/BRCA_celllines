#----------------------------------------------------------------------------
# 04_survival.R
#----------------------------------------------------------------------------
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
se <- readRDS("rds/TCGA_BRCA_RNAExp.rds")
assays(se)$tpm_unstrand_log2 <- log2(assays(se)$tpm_unstrand+1)

#----- subtype -----#
se_subtype <- list(LumA = se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "LumA")],
                   LumB = se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "LumB")],
                   Her2 = se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "Her2")],
                   Basal = se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "Basal")],
                   Normal = se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "Normal")])

#----- ENTREZ ID -----#
ensg_FOXA1 <- rownames(se)[which(rowData(se)$gene_name == "FOXA1")]
ensg_GRHL2 <- rownames(se)[which(rowData(se)$gene_name == "GRHL2")]

#----- survival analysis -----#
p1 <- list()
for(i in seq_along(se_subtype)){
  x <- names(se_subtype)[i]
  tmp_se <- se_subtype[[x]]

  sorted_exp1 <- assays(tmp_se)$tpm_unstrand_log2[ensg_FOXA1,] %>% sort()
  sorted_exp2 <- assays(tmp_se)$tpm_unstrand_log2[ensg_FOXA1,] %>% sort(.,decreasing = T)
  
  num <- round(ncol(tmp_se)*0.33)
  low_group <- sorted_exp1[c(1:num)] %>% names()
  high_group <- sorted_exp2[c(1:num)] %>% names()
    
  surv <- data.frame(sample = colnames(tmp_se),
                     daysFO = tmp_se$days_to_last_follow_up,
                     daysDE = tmp_se$days_to_death,
                     status = tmp_se$vital_status)
  surv$group <- NA
  surv$group[surv$sample %in% high_group] <- "High"
  surv$group[surv$sample %in% low_group] <- "Low"
  
  surv <- surv[!is.na(surv$group),]
  
  surv$days <- surv$daysFO
  surv$days[is.na(surv$days)] <- surv$daysDE[is.na(surv$days)]
  surv$status2 <- 0
  surv$status2[surv$status == "Dead"] <- 1 
  sf <- survfit(Surv(surv$days, surv$status2)~surv$group)
  p <- ggsurvplot(fit = sf, data = surv,
                  pval = TRUE, pval.method = TRUE,
                  risk.table = TRUE, conf.int = FALSE,
                  ncensor.plot = FALSE, size = 1.5, #linetype = c(1, 3),
                  title = x,
                  legend.title = "Expression Pattern",
                  log.rank.weights = "1",
                  risk.table.title = "",
                  risk.table.y.text.col = TRUE,
                  risk.table.y.text = FALSE)
  p1[[i]] <- p
}

p2 <- list()
for(i in seq_along(se_subtype)){
  x <- names(se_subtype)[i]
  tmp_se <- se_subtype[[x]]
  
  sorted_exp1 <- assays(tmp_se)$tpm_unstrand_log2[ensg_GRHL2,] %>% sort()
  sorted_exp2 <- assays(tmp_se)$tpm_unstrand_log2[ensg_GRHL2,] %>% sort(.,decreasing = T)
  
  num <- round(ncol(tmp_se)*0.33)
  low_group <- sorted_exp1[c(1:num)] %>% names()
  high_group <- sorted_exp2[c(1:num)] %>% names()
  
  surv <- data.frame(sample = colnames(tmp_se),
                     daysFO = tmp_se$days_to_last_follow_up,
                     daysDE = tmp_se$days_to_death,
                     status = tmp_se$vital_status)
  surv$group <- NA
  surv$group[surv$sample %in% high_group] <- "High"
  surv$group[surv$sample %in% low_group] <- "Low"
  
  surv <- surv[!is.na(surv$group),]
  
  surv$days <- surv$daysFO
  surv$days[is.na(surv$days)] <- surv$daysDE[is.na(surv$days)]
  surv$status2 <- 0
  surv$status2[surv$status == "Dead"] <- 1 
  sf <- survfit(Surv(surv$days, surv$status2)~surv$group)
  p <- ggsurvplot(fit = sf, data = surv,
                  pval = TRUE, pval.method = TRUE,
                  risk.table = TRUE, conf.int = FALSE,
                  ncensor.plot = FALSE, size = 1.5, #linetype = c(1, 3),
                  title = x,
                  legend.title = "Expression Pattern",
                  log.rank.weights = "1",
                  risk.table.title = "",
                  risk.table.y.text.col = TRUE,
                  risk.table.y.text = FALSE)
  p2[[i]] <- p
}

pdf("output/Plots/06_SurvPlot_FOXA1_33perc.pdf")
p1
dev.off()
pdf("output/Plots/06_SurvPlot_GRHL2_33perc.pdf")
p2
dev.off()
