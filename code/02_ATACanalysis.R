#------------------------------------------------------------------------------
# 02_ATACanalysis.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(ComplexHeatmap)
library(circlize)
library(viridis)
source("code/edgeR_PairwiseFunction.R")
se <- readRDS("rds/atac_se.rds")

#--------- Annotate Peaks ---------#
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(rowRanges(se), tssRegion=c(-1000, 100), TxDb = txdb, annoDb = "org.Hs.eg.db")
mcols(se) <- mcols(as.GRanges(peakAnno))

pdf("output/Plots/02_peakAnnoPie.pdf", width = 7, height = 4)
plotAnnoPie(peakAnno) 
dev.off()

#saveRDS(se, "rds/atac_se.rds")

#--------- PCA ---------#
pca1 <- prcomp(t(assays(se)$normcounts))
summary(pca1) # PC1: 0.2295, PC2: 0.1219
df <- pca1$x[, c(1:2)] %>% data.frame
df$sample <- se$sample
df$subtype <- se$subtype

p1 <- ggplot(df, aes(x = PC1, y = PC2, color = subtype, label = sample)) + geom_point() + ArchR::theme_ArchR() + 
  scale_color_manual(values = c("TN" = "red", "ER+/HER2-" = "blue", "ER-/HER2+" = "orange", "ER+/HER2+" = "purple")) +
  geom_text_repel(size = 2) + labs(x = "PC1 (23.0% variance)", y = "PC2 (12.2% variance)")
pdf("output/Plots/02_PCA.pdf", width = 4, height = 4.5)
p1
dev.off()

#--------- Correlation analysis ---------#
corr_promoter <- cor(assays(se[mcols(se)$annotation == "Promoter", ])$normcounts)
corr_distal   <- cor(assays(se[mcols(se)$annotation != "Promoter", ])$normcounts)
corr_total    <- cor(assays(se)$normcounts)

ha1 <- HeatmapAnnotation(Subtype = colData(se)$subtype, col = list(Subtype = c("TN" = "red", "ER+/HER2-" = "blue", "ER+/HER2+" = "purple", "ER-/HER2+" = "orange")))
fh = function(x) hclust(dist(x), method="ward.D2")
col_fun1 <- colorRamp2(c(0, 0.5, 1.0), c("white","orange", "red"))
ht1 <- Heatmap(corr_promoter, name = "Pearson's correlation", col = col_fun1, cluster_rows = fh, cluster_columns = fh,
               show_column_names = F, show_column_dend = F, show_row_names = T, show_row_dend = T, row_names_gp = gpar(fontsize = 10), 
               column_title = paste0("Promoter elements (N = ", length(which(mcols(se)$annotation == "Promoter")),")"),
               column_title_side = "top", row_title = "", top_annotation = ha1)
p2 <- draw(ht1)
ht2 <- Heatmap(corr_distal, name = "Pearson's correlation", col = col_fun1, cluster_rows = fh, cluster_columns = fh,
               show_column_names = F, show_column_dend = F, show_row_names = T, show_row_dend = T, row_names_gp = gpar(fontsize = 10), 
               column_title = paste0("Distal elements (N = ", length(which(mcols(se)$annotation != "Promoter")),")"),
               column_title_side = "top", row_title = "", top_annotation = ha1)
p3 <- draw(ht2)
ht3 <- Heatmap(corr_total, name = "Pearson's correlation", col = col_fun1, cluster_rows = fh, cluster_columns = fh,
               show_column_names = F, show_column_dend = F, show_row_names = T, show_row_dend = T, row_names_gp = gpar(fontsize = 10), 
               column_title = "Total",column_title_side = "top", row_title = "", top_annotation = ha1)
p4 <- draw(ht3)
pdf("output/Plots/02_corHeatmap.pdf", width = 6, height = 3)
p2
p3
p4
dev.off()

#--------- Most variable 5000 regions ---------#
var_se <- rowVars(assays(se)$normcounts)
var5000_idx <- order(var_se, decreasing = T)[c(1:5000)]

mtx <- assays(se)$normcounts[var5000_idx,]
mtx_z <- t(scale(t(mtx)))
col_fun2 <- colorRamp2(c(-2,-1,0,1,2,3), viridis::viridis(6))
ht4 <- Heatmap(mtx_z, 
               name = "z-score[Normalized ATAC counts]", col = col_fun2, 
               cluster_rows = fh, cluster_columns = fh, top_annotation = ha1,
               show_column_names = T, show_column_dend = T, show_row_names = F, show_row_dend = F, 
               #row_names_gp = gpar(fontsize = 10), 
               row_title_gp = gpar(fontsize = 8), 
               column_title = "", column_title_side = "top", row_title = "", use_raster = T)
p5 <- draw(ht4)
pdf("output/Plots/02_Heatmap_var5000_Zscore.pdf", width = 6, height = 6)
p5
dev.off()

#--------- classify cell lines ---------#
group <- c(rep("B", 9), rep("M",6), rep("D",7), NA)
names(group) <- colData(se)$sample[row_order(p3)] 
colData(se)$group <- group[colnames(se)]

#--------- Identify differential accessible regions ---------#
#differential analysis
diff_ls <- list(D = list("D", c("M", "B")), 
                B = list("B", c("D", "M")), 
                M = list("M", c("D", "B")))
DiffTest <- lapply(diff_ls, function(x) edgeR_pairwise(se, compareCol = "group", topGroup = x[[1]], bottomGroup = x[[2]]))

DAR_MB <- rownames(DiffTest[[1]])[which(assay(DiffTest[[1]])[,"log2FC"] < -1 & assay(DiffTest[[1]])[,"FDR"] < 0.01)]
DAR_DM <- rownames(DiffTest[[2]])[which(assay(DiffTest[[2]])[,"log2FC"] < -1 & assay(DiffTest[[2]])[,"FDR"] < 0.01)]
DAR_DB <- rownames(DiffTest[[3]])[which(assay(DiffTest[[3]])[,"log2FC"] < -1 & assay(DiffTest[[3]])[,"FDR"] < 0.01)]

DAR_D <- rownames(DiffTest[[1]])[which(assay(DiffTest[[1]])[,"log2FC"] > 1 & assay(DiffTest[[1]])[,"FDR"] < 0.01)]
DAR_B <- rownames(DiffTest[[2]])[which(assay(DiffTest[[2]])[,"log2FC"] > 1 & assay(DiffTest[[2]])[,"FDR"] < 0.01)]
DAR_M <- rownames(DiffTest[[3]])[which(assay(DiffTest[[3]])[,"log2FC"] > 1 & assay(DiffTest[[3]])[,"FDR"] < 0.01)]

DAR_D_sp <- DAR_D[DAR_D %ni% unique(c(intersect(DAR_D, DAR_DM), intersect(DAR_D, DAR_DB)))]
DAR_B_sp <- DAR_B[DAR_B %ni% unique(c(intersect(DAR_B, DAR_MB), intersect(DAR_B, DAR_DB)))]
DAR_M_sp <- DAR_M[DAR_M %ni% unique(c(intersect(DAR_M, DAR_MB), intersect(DAR_M, DAR_DM)))]

#check exclusiveness
#gplots::venn(list(DAR_D_sp, DAR_B_sp, DAR_M_sp))

#export
peaks <- rowRanges(se)
gr <- GRangesList(DAR_D_specific = peaks[DAR_D_sp], 
                  DAR_B_specific = peaks[DAR_B_sp], 
                  DAR_M_specific = peaks[DAR_M_sp],
                  DAR_MB = peaks[DAR_MB], 
                  DAR_DM = peaks[DAR_DM],
                  DAR_DB = peaks[DAR_DB])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/DAR/", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#DAR heatmap
DAR_ls <- list(D_specific = DAR_D_sp, 
               B_specific = DAR_B_sp, 
               M_specific = DAR_M_sp,
               MB_shared = DAR_MB, 
               DM_shared = DAR_DM,
               DB_shared = DAR_DB)
mtx <- assays(se)$normcounts[unlist(DAR_ls), ]
mtx_z <- t(scale(t(mtx)))

col_fun4 <- colorRamp2(c(-3,-2,-1,0,1,2,3), viridis::viridis(7))
ht6 <- Heatmap(mtx_z[,idy], 
               name = "z-score[Normalized ATAC counts]", col = col_fun4, 
               cluster_rows = F, cluster_columns = fh, top_annotation = ha2,
               show_column_names = T, show_column_dend = T, show_row_names = F, show_row_dend = F, 
               column_names_gp = gpar(fontsize = 8), 
               row_title_gp = gpar(fontsize = 8), row_title_rot = 0,
               row_split = unlist(lapply(DAR_ls, length)) %>% rep(names(.), .),
               column_split = colData(se)[idy,]$group,
               column_title = "", column_title_side = "top", row_title = "", use_raster = T)
p8 <- draw(ht6)

pdf("output/Plots/02_HeatmapDAR.pdf", width = 6, height = 5)
p8
dev.off()

#--------- GREAT visualization ---------#
great_file <- list.files("output/great_go/DAR/", pattern = "GREAT_")
GREAT_df <- lapply(great_file, function(i){
  out <- data.table::fread(paste0("output/great_go/DAR/", i), header = F, skip = 1)[, c(1,4)] %>% 
    data.frame %>% `colnames<-`(., c("Term", "FDR"))
  out$log10FDR <- -log10(out$FDR)
  return(out)
})
names(GREAT_df) <- great_file

p9 <- lapply(names(GREAT_df), 
             function(i) ggplot(GREAT_df[[i]], aes(x = log10FDR, y = reorder(Term, log10FDR))) + 
               geom_bar(stat = "identity", fill = "gray") + ArchR::theme_ArchR() + ggtitle(i))

pdf("output/Plots/02_GREAT_DAR.pdf", height = 5, width = 10)
p9
dev.off()

#----------------- chromVAR -----------------#
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
motifs <- chromVARmotifs::homer_pwms
motif_ix <- matchMotifs(motifs, se, genome = BSgenome.Hsapiens.UCSC.hg38)
bg <- getBackgroundPeaks(object = se)
dev <- computeDeviations(object = se, annotations = motif_ix, background_peaks = bg)

saveRDS(dev, "rds/chromvar_dev.rds")

df <- data.frame(sample = colnames(dev), group = factor(colData(dev)$group, levels = c("D","B","M")),
                 ERE = assays(dev)$z[grep("ERE", rownames(dev)),], 
                 GRHL2 = assays(dev)$z[grep("GRHL2", rownames(dev)),],
                 FOXA1 = assays(dev)$z["FOXA1(Forkhead)/MCF7-FOXA1-ChIP-Seq(GSE26831)/Homer",],
                 GATA3 = assays(dev)$z["GATA3(Zf)/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer",])
df <- df[!is.na(df$group),]
p10 <- ggplot(df, aes(x = group, y = ERE, fill = group)) + geom_boxplot(alpha = 0.5) + geom_point() + 
       scale_fill_manual(values = c("D"="#377EB8","B"="#E41A1C","M"="#4DAF4A")) + ArchR::theme_ArchR() +
       labs(y = "ERE ChromVAR motif score", x = "Group") + 
       geom_signif(comparisons = list(c("D","B"), c("B","M"), c("D","M")), test = "t.test", textsize = 2)
p11 <- ggplot(df, aes(x = group, y = GRHL2, fill = group)) + geom_boxplot(alpha = 0.5) + geom_point() + 
       scale_fill_manual(values = c("D"="#377EB8","B"="#E41A1C","M"="#4DAF4A")) + ArchR::theme_ArchR() +
       labs(y = "GRHL2 ChromVAR motif score", "Group") +
       geom_signif(comparisons = list(c("D","B"), c("B","M"), c("D","M")), test = "t.test", textsize = 2)
p12 <- ggplot(df, aes(x = group, y = FOXA1, fill = group)) + geom_boxplot(alpha = 0.5) + geom_point() + 
  scale_fill_manual(values = c("D"="#377EB8","B"="#E41A1C","M"="#4DAF4A")) + ArchR::theme_ArchR() +
  labs(y = "FOXA1 ChromVAR motif score", "Group") +
  geom_signif(comparisons = list(c("D","B"), c("B","M"), c("D","M")), test = "t.test", textsize = 2)
p13 <- ggplot(df, aes(x = group, y = GATA3, fill = group)) + geom_boxplot(alpha = 0.5) + geom_point() + 
  scale_fill_manual(values = c("D"="#377EB8","B"="#E41A1C","M"="#4DAF4A")) + ArchR::theme_ArchR() +
  labs(y = "GATA3 ChromVAR motif score", "Group") +
  geom_signif(comparisons = list(c("D","B"), c("B","M"), c("D","M")), test = "t.test", textsize = 2)

pdf("output/Plots/02_ChromVARscore_v1.pdf", height = 5, width = 3)
p10
p11
p12
p13
dev.off()

#heatmap
idx <- na.omit(colData(se)$group)

mtx <- assays(dev)$z
var_dev <- rowVars(mtx)
mtx <- mtx[order(var_dev, decreasing = T)[c(1:50)],names(idx)]

col_fun5 <- colorRamp2(c(-20,-10,0,10,20), c("#1984c5","#63bff0","#e2e2e2","#de6e56", "#c23728"))
ht7 <- Heatmap(mtx,
               name = "ChromVAR score", col = col_fun5, 
               cluster_rows = fh, cluster_columns = fh, top_annotation = ha2,
               show_column_names = T, show_column_dend = T, show_row_names = T, show_row_dend = F, 
               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
               row_title_gp = gpar(fontsize = 8),
               column_split = idx,
               column_title = "", column_title_side = "top", row_title = "Most variable 30 enriched motifs", use_raster = F)
p14 <- draw(ht7)

pdf("output/Plots/02_Heatmap_var50_chromVAR.pdf", width = 10, height = 6)
p14
dev.off()
