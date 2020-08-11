setwd("M:/")

RosmapCov <- read.table("covSexAgeBatchPCs.txt", header = TRUE)

# SNPfile <- readRDS("../mqtlAnalysis_Feb2020/rosmap12Merged_CpGoverlap_MAF1_dosage_AllChrs.rds")
# 
# resRosmap <- readRDS("C:/Users/lxg255/Meta_analysis_code/DATASETS/ROSMAP/step10_residuals/ROSMAP_QNBMIQ_PCfiltered_mvalResiduals.RDS")
# 
# overlap <- readRDS("../mqtlAnalysis_Feb2020/CpGimpVarOverlap.rds")

## Filter for AD associated sNPs
# adSNP <- read.table("glmCogdx_Rosmap12_CpGoverlap_MAF1/glmCogdx_ResultsFreqAll_CpGs.txt",
#                     header = TRUE)
# adSNPsign <- adSNP[which(adSNP$P < 0.05),]
# write.csv(adSNPsign, "ADcogdx_associatedSNPs_CpGs.csv", row.names = FALSE)
# 
# SNPfile_ADsign <- SNPfile[which(SNPfile$SNP %in% adSNPsign$ID), ]
# saveRDS(SNPfile_ADsign, "SNPfile_ADsign.rds")
SNPfile_ADsign <- readRDS("SNPfile_ADsign.rds")

# overlap_ADsign <- overlap[which(overlap$SNP %in% adSNPsign$ID), ]
# saveRDS(overlap_ADsign, "overlap_ADsign.rds")
overlap_ADsign <- readRDS("overlap_ADsign.rds")
CpGs_ADsign <- unique(overlap_ADsign$cpg)

sampleKey <- read.table("../mqtlAnalysis_Feb2020/sampleKey.txt", header = TRUE)

## Filter residual matrix for 688 individuals 
## included in the mqtl analysis and 3679 sign CpGs
# resCpGs_ADsign <- resRosmap[which(rownames(resRosmap)%in%CpGs_ADsign),
#                             which(colnames(resRosmap)%in%sampleKey$methySampleID)]
# saveRDS(resCpGs_ADsign, "residuals_3679CpGs_688individuals_ADsign.rds")
resCpGs_ADsign <- readRDS("residuals_3679CpGs_688individuals_ADsign.rds")

## Match Ids in methy and SNP datasets
methyIDs <- data.frame(IDs = colnames(resCpGs_ADsign))
methyIDsKey <- merge(methyIDs, 
                     sampleKey,
                     by.x = "IDs",
                     by.y = "methySampleID",
                     sort = FALSE)
identical(methyIDs$IDs, methyIDsKey$IDs)
colnames(resCpGs_ADsign) <- methyIDsKey$SNPsampleID

resMethyCpGs_ADsign <- resCpGs_ADsign[, as.character(RosmapCov$FID)]

CpGs <- rownames(resMethyCpGs_ADsign)

identical(as.vector(RosmapCov$FID), colnames(SNPfile_ADsign)[-1])
identical(as.vector(RosmapCov$FID), colnames(resMethyCpGs_ADsign))

## mQTL analysis
## Methylation_resid ~ SNP dosage + batch + PC1 + PC2 + PC3
lmodel1 <- function(SNP) {
  
  f <- lm(CpGmethy ~ SNP + RosmapCov$batch + RosmapCov$PC1 + RosmapCov$PC2 + RosmapCov$PC3)
  
  res <- t(coef(summary(f))[2, c(1, 2, 4)])
  
  res
  
}


all_m1 <- data.frame(matrix(ncol=6,nrow=0))
for (i in CpGs[1:length(CpGs)]){
  
  CpGmethy <- as.numeric(as.vector(resMethyCpGs_ADsign[i,]))
  CpG_SNPs <- subset(overlap_ADsign, cpg == i, select = SNP)
  CpG_SNPs <- CpG_SNPs$SNP
  SNPdosage <- SNPfile_ADsign[which(SNPfile_ADsign$SNP%in%CpG_SNPs), ]
  CpG_SNP <- cbind(paste0(i, "-", SNPdosage$SNP), i, SNPdosage$SNP)
  
  m1 <- cbind(CpG_SNP, as.data.frame(t(apply(SNPdosage[,-1], 1, lmodel1))))
  all_m1 <- rbind(all_m1, m1)
  
}


colnames(all_m1) <- c("CpG-SNP", "CpG", "SNP", "Estimate_SNP", "StdError_SNP", "P_SNP")
all_m1$FDR <- p.adjust(all_m1$P_SNP, "fdr")
write.csv(all_m1, "ADmqtlAnalysisResults_CpGs_mQTLsALL_ADcogdxSign.csv")
saveRDS(all_m1, "ADmqtlAnalysisResults_CpGs_mQTLsALL_ADcogdxSign.rds")
