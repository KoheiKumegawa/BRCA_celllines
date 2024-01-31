#------------------------------------------------------------------------------
# 05_RNAanalysis.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)
countToTpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#---- make se -----#
counts <- fread("data/counts_exon_uns.txt")
len <- counts$Length
gene <- counts$Geneid
counts <- counts[, -c(1:6)] %>% 
  `colnames<-`(., gsub("_hg38_Aligned.sortedByCoord.out.bam", "", colnames(counts[, -c(1:6)]))) %>%
  as.matrix %>% `rownames<-`(., gene) 
tpms <- countToTpm(counts, len)
se <- SummarizedExperiment(assays = list(counts = counts, log2tpm = log2(tpms+1)),
                           rowData = DataFrame(symbol = gene),
                           colData = DataFrame(sample = colnames(counts)))

colData(se)$sample2 <- c("HCC38_siCtrl", "HCC38_siGRHL2", "T47D_siCtrl", "T47D_siFOXA1", "T47D_siGRHL2")
colnames(se) <- colData(se)$sample2 
saveRDS(se, "rds/05_rna_se.rds")

#---- MA plot -----#
T47D_siFOXA1_MA <- data.frame(log2FC = assays(se)$log2tpm[,"T47D_siFOXA1"] - assays(se)$log2tpm[,"T47D_siCtrl"],
                              Average = (assays(se)$log2tpm[,"T47D_siFOXA1"] + assays(se)$log2tpm[,"T47D_siCtrl"])/2)
T47D_siGRHL2_MA <- data.frame(log2FC = assays(se)$log2tpm[,"T47D_siGRHL2"] - assays(se)$log2tpm[,"T47D_siCtrl"],
                              Average = (assays(se)$log2tpm[,"T47D_siGRHL2"] + assays(se)$log2tpm[,"T47D_siCtrl"])/2)
HCC38_siGRHL2_MA <- data.frame(log2FC = assays(se)$log2tpm[,"HCC38_siGRHL2"] - assays(se)$log2tpm[,"HCC38_siCtrl"],
                              Average = (assays(se)$log2tpm[,"HCC38_siGRHL2"] + assays(se)$log2tpm[,"HCC38_siCtrl"])/2)

SigGenes <- list(T47D_siFOXA1_UP = rownames(T47D_siFOXA1_MA)[T47D_siFOXA1_MA$log2FC > 1 & T47D_siFOXA1_MA$Average > 2] %>% sort,
                 T47D_siFOXA1_DN = rownames(T47D_siFOXA1_MA)[T47D_siFOXA1_MA$log2FC < -1 & T47D_siFOXA1_MA$Average > 2] %>% sort,
                 T47D_siGRHL2_UP = rownames(T47D_siGRHL2_MA)[T47D_siGRHL2_MA$log2FC > 1 & T47D_siGRHL2_MA$Average > 2] %>% sort,
                 T47D_siGRHL2_DN = rownames(T47D_siGRHL2_MA)[T47D_siGRHL2_MA$log2FC < -1 & T47D_siGRHL2_MA$Average > 2] %>% sort,
                 HCC38_siGRHL2_UP = rownames(HCC38_siGRHL2_MA)[HCC38_siGRHL2_MA$log2FC > 1 & HCC38_siGRHL2_MA$Average > 2] %>% sort,
                 HCC38_siGRHL2_DN = rownames(HCC38_siGRHL2_MA)[HCC38_siGRHL2_MA$log2FC < -1 & HCC38_siGRHL2_MA$Average > 2] %>% sort)

T47D_siFOXA1_MA$Sig <- "X"
T47D_siFOXA1_MA$Sig[T47D_siFOXA1_MA$log2FC > 1 & T47D_siFOXA1_MA$Average > 2] <- "UP"
T47D_siFOXA1_MA$Sig[T47D_siFOXA1_MA$log2FC < -1 & T47D_siFOXA1_MA$Average > 2] <- "DN"
T47D_siFOXA1_MA <- T47D_siFOXA1_MA[order(T47D_siFOXA1_MA$Sig, decreasing = T),]

T47D_siGRHL2_MA$Sig <- "X"
T47D_siGRHL2_MA$Sig[T47D_siGRHL2_MA$log2FC > 1 & T47D_siGRHL2_MA$Average > 2] <- "UP"
T47D_siGRHL2_MA$Sig[T47D_siGRHL2_MA$log2FC < -1 & T47D_siGRHL2_MA$Average > 2] <- "DN"
T47D_siGRHL2_MA <- T47D_siGRHL2_MA[order(T47D_siGRHL2_MA$Sig, decreasing = T),]

HCC38_siGRHL2_MA$Sig <- "X"
HCC38_siGRHL2_MA$Sig[HCC38_siGRHL2_MA$log2FC > 1 & HCC38_siGRHL2_MA$Average > 2] <- "UP"
HCC38_siGRHL2_MA$Sig[HCC38_siGRHL2_MA$log2FC < -1 & HCC38_siGRHL2_MA$Average > 2] <- "DN"
HCC38_siGRHL2_MA <- HCC38_siGRHL2_MA[order(HCC38_siGRHL2_MA$Sig, decreasing = T),]

T47D_siFOXA1_MA$label <- ""
T47D_siFOXA1_MA$label[which(rownames(T47D_siFOXA1_MA) == "FOXA1")] <- "FOXA1"
T47D_siFOXA1_MA$label[which(rownames(T47D_siFOXA1_MA) == "GRHL2")] <- "GRHL2"

T47D_siGRHL2_MA$label <- ""
T47D_siGRHL2_MA$label[which(rownames(T47D_siGRHL2_MA) == "FOXA1")] <- "FOXA1"
T47D_siGRHL2_MA$label[which(rownames(T47D_siGRHL2_MA) == "GRHL2")] <- "GRHL2"

HCC38_siGRHL2_MA$label <- ""
HCC38_siGRHL2_MA$label[which(rownames(HCC38_siGRHL2_MA) == "FOXA1")] <- "FOXA1"
HCC38_siGRHL2_MA$label[which(rownames(HCC38_siGRHL2_MA) == "GRHL2")] <- "GRHL2"

p1 <- ggplot(T47D_siFOXA1_MA, aes(x = Average, y = log2FC, color = Sig, label = label)) + geom_point_rast() + 
  geom_text_repel() + ArchR::theme_ArchR() + scale_color_manual(values = c("UP" = "red", "DN" = "blue", "X" = "darkgray")) +
  ggtitle("T47D FOXA1-KD 48h") + geom_hline(yintercept = 0, lty = "dotted")
p2 <- ggplot(T47D_siGRHL2_MA, aes(x = Average, y = log2FC, color = Sig, label = label)) + geom_point_rast() + 
  geom_text_repel() + ArchR::theme_ArchR() + scale_color_manual(values = c("UP" = "red", "DN" = "blue", "X" = "darkgray")) +
  ggtitle("T47D GRHL2-KD 48h") + geom_hline(yintercept = 0, lty = "dotted")
p3 <- ggplot(HCC38_siGRHL2_MA, aes(x = Average, y = log2FC, color = Sig, label = label)) + geom_point_rast() + 
  geom_text_repel(max.overlaps = 100) + ArchR::theme_ArchR() + scale_color_manual(values = c("UP" = "red", "DN" = "blue", "X" = "darkgray")) +
  ggtitle("HCC38 GRHL2-KD 48h") + geom_hline(yintercept = 0, lty = "dotted")

pdf("output/Plots/05_MAplot.pdf", width = 6, height = 6)
p1
p2
p3
dev.off()

lapply(names(SigGenes), function(i) write.table(SigGenes[[i]], paste0("output/Tables/05_", i, ".txt"), quote = F, row.names = F, col.names = F))

#---- GO -----#
library(clusterProfiler)
library(stringr)

#gene_list
SigGenes

#gene set
gmt_files <- paste0("ref/", list.files("ref", pattern = ".txt"))
gmt_files

#function
EnrichmentAnalysis <- function(gmt_file, genelist_file, DIR) {
  # gmt_file: an item in gmt_files
  # genelist_file: gene list of interest
  # DIR: directory to export csv and pdf files
  
  # import gmt file
  enrichr_gmt <- read.gmt(gmt_file)
  # get gmt file name
  gmt_file_strip <- str_split(gmt_file, "/", simplify = T)
  gmt_file_name <- str_split(gmt_file_strip[length(gmt_file_strip)], "\\.", simplify = T)[1]
  print(paste0("calc enrichment with ", gmt_file_name, " ..."))
  
  # calc. enrichment score by culsterProfiler
  CPfile <- compareCluster(geneCluster = genelist_file, fun='enricher', TERM2GENE=enrichr_gmt, pvalueCutoff = 0.2) #pvalueCutoff=0.05
  
  # export graphs showCategory
  for (j in c(1,2)) {
    w <- myWidthAlgorithm(CPfile, j)
    g <- dotplot(CPfile,  showCategory=j) + ggtitle(gmt_file_name) + 
      theme(axis.text.y = element_text(size = 6, lineheight = 0.6),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
    ggsave(plot=g, filename=paste0(DIR, "/", gmt_file_name, "_show_", j, ".pdf"), 
           width=w, h=6)
    # showCategory parameter: most significant categories for each Cluster
    # see https://guangchuangyu.github.io/2016/11/showcategory-parameter-for-visualizing-comparecluster-output/
  }
  
  # export table
  # write_csv(as.data.frame(CPfile), paste0(DIR, "/Tables/", gmt_file_name,".csv"))
}

myWidthAlgorithm <- function(CPfile, n) {
  df <- CPfile@compareClusterResult
  df_show <- df %>% group_by(Cluster) %>% top_n(-qvalue, n = n)
  maxnchar <- max(nchar(df_show$ID))
  nclu <- df$Cluster %>% unique() %>% length()
  return(nclu*0.8 + maxnchar/20) 
}

EnrichmentAnalysis(gmt_file = gmt_files[1], genelist_file = SigGenes, DIR = "output/Plots/")
EnrichmentAnalysis(gmt_file = gmt_files[2], genelist_file = SigGenes, DIR = "output/Plots/")
