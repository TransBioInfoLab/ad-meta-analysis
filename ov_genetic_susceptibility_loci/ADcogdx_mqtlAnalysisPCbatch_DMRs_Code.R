setwd("M:/")

RosmapCov <- read.table("covSexAgeBatchPCs.txt", header = TRUE)

SNPfile <- readRDS("../mqtlAnalysis_Feb2020/rosmap12Merged_DMRoverlap_MAF1_dosage_AllChrs.rds")

# resRosmap <- readRDS("C:/Users/lxg255/ROSMAP_QNBMIQ_PCfiltered_mvalResiduals.RDS")

regionCpGall <- readRDS("../mqtlAnalysis_Feb2020/DMRs_Probes_list.rds")
DMRs <- unique(as.character(regionCpGall$region))

# sampleKey <- read.table("../mqtlAnalysis_Feb2020/sampleKey.txt", header = TRUE)

## Filter residual matrix for 688 individuals 
## included in the mqtl analysis and 608 DMR cpgs
# resDMRcpgs <- resRosmap[which(rownames(resRosmap)%in%regionCpGall$cpg),
#                         which(colnames(resRosmap)%in%sampleKey$methySampleID)]
# saveRDS(resDMRcpgs, "residuals_608DMRcpgs_688individuals.rds")
# resDMRcpgs <- readRDS("../mqtlAnalysis_Feb2020/residuals_608DMRcpgs_688individuals.rds")

## Match Ids in methy and SNP datasets
# methyIDs <- data.frame(IDs = colnames(resDMRcpgs))
# methyIDsKey <- merge(methyIDs, sampleKey, by.x =  "IDs", by.y = "methySampleID", sort = FALSE)
# identical(methyIDs$IDs, methyIDsKey$IDs)
# colnames(resDMRcpgs) <- methyIDsKey$SNPsampleID


## Calculate median (Methylation_resid) 
# DMRmedian118<- data.frame(matrix(ncol = 688, nrow = 0))
# for (i in DMRs[1:118]){
#   
#   CpGs <- regionCpGall[which(regionCpGall$region%in%i), "cpg"]
#   CpGresidual <- resDMRcpgs[which(rownames(resDMRcpgs)%in%CpGs),]
#   DMRmedian <- apply(CpGresidual, 2, median)
#   DMRmedian118 <- rbind(DMRmedian118, DMRmedian)
#   
# }
# 
# rownames(DMRmedian118) <- DMRs
# colnames(DMRmedian118) <- names(DMRmedian)
# saveRDS(DMRmedian118, "Rosmap_118DMRs_residualMedians.rds")
DMRmedian118 <- readRDS("../mqtlAnalysis_Feb2020/Rosmap_118DMRs_residualMedians.rds")
RosmapMvalMediansDMRs <- DMRmedian118[,as.character(RosmapCov$FID)]

overlap <- readRDS("../mqtlAnalysis_Feb2020/DMRimpVarOverlap.rds")

## Filter for AD associated sNPs
adSNP <- read.table("glmCogdx_Rosmap12_DMRoverlap_MAF1/glmCogdx_ResultsFreqAll_DMRs.txt",
                    header = TRUE)
adSNPsign <- adSNP[which(adSNP$P < 0.05),]
SNPfile_ADsign <- SNPfile[which(SNPfile$SNP %in% adSNPsign$ID), ]
overlap_ADsign <- overlap[which(overlap$SNP %in% adSNPsign$ID), ]
write.csv(adSNPsign, "ADcogdx_associatedSNPs.csv", row.names = FALSE)


## Remove "chr20:62198872-62199190" from the list of DMRs 
## because it doesn't overlap with any SNP
DMRs <- DMRs[DMRs != "chr20:62198872-62199190"]


## mQTL analysis
## Methylation_resid ~ SNP dosage + batch + PC1 + PC2 + PC3
lmodel1 <- function(SNP) {
  
  f <- lm(DMRmethy ~ SNP + RosmapCov$batch + RosmapCov$PC1 + RosmapCov$PC2 + RosmapCov$PC3)
  
  res <- t(coef(summary(f))[2, c(1, 2, 4)])
  
  res
  
}


identical(as.vector(RosmapCov$FID), colnames(SNPfile)[-1])
identical(as.vector(RosmapCov$FID), colnames(RosmapMvalMediansDMRs))

all_m1 <- data.frame(matrix(ncol=6,nrow=0))
for (i in DMRs[1:length(DMRs)]){
  DMRmethy <- as.numeric(as.vector(RosmapMvalMediansDMRs[i,]))
  DMR_SNPs <- subset(overlap_ADsign, DMR == i, select = SNP)
  DMR_SNPs <- DMR_SNPs$SNP
  SNPdosage <- SNPfile_ADsign[which(SNPfile_ADsign$SNP%in%DMR_SNPs), ]
  DMR_SNP <- cbind(paste0(i, "-", SNPdosage$SNP), i, SNPdosage$SNP)
  
  m1 <- cbind(DMR_SNP, as.data.frame(t(apply(SNPdosage[,-1], 1, lmodel1))))
  all_m1 <- rbind(all_m1, m1)
  
}


colnames(all_m1) <- c("DMR-SNP", "DMR", "SNP", "Estimate_SNP", "StdError_SNP", "P_SNP")
all_m1$FDR <- p.adjust(all_m1$P_SNP, "fdr")
write.csv(all_m1, "ADmqtlAnalysisResults_DMRs_mQTLsALL_ADcogdxSign.csv")
saveRDS(all_m1, "ADmqtlAnalysisResults_DMRs_mQTLsALL_ADcogdxSign.rds")

