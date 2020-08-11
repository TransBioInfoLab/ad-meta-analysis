## mQTL analysis for co-localisation test 
## Rosmap dataset (n=688)
## 24 CpGs that overlap with AD GWAS loci
## Test association of each CpG with all the SNPs in the GWAS locus in cis with the CpG
## Lissette Gomez, July 2020


setwd("M:/")

### Import input files

## Covariates file
# RosmapPhenoAll <- read.csv("RosmapCovTransposed_caseControls.csv")
# rownames(RosmapPhenoAll)<-RosmapPhenoAll[,1]
# RosmapPhenoAll<-RosmapPhenoAll[,-1]
# saveRDS(RosmapPhenoAll, "RosmapPhenoAll.rds")

RosmapPheno <- readRDS("RosmapPhenoAll.rds")
Ids <- colnames(RosmapPheno)

## SNP file
# SNPfile <- read.table(
#     "M:/rosmap12Merged_GWASlociSNPs_MAF1_dosage_AllChr.gen",
#     header = TRUE
# )
# SNPfile <- SNPfile[, Ids]
# saveRDS(SNPfile, "rosmap12Merged_GWASlociSNPs_MAF1_dosage_AllChrs.rds")
SNPfile <- readRDS("rosmap12Merged_GWASlociSNPs_MAF1_dosage_AllChrs.rds")

SNPids <- read.table("SNP_alleles.txt")
colnames(SNPids) <- c("SNPchrPos", "SNPid")

## SNPs/CpGs location
GWASlociSNPs_CpG <- readRDS("GWASlociSNPs_CpGs_pairs.rds")
#24 AD CpGs that overlap with GWAS loci
CpG24 <- read.table("overlapGWASloci_24CpGs.txt", header = TRUE) 
GWASlociSNPs_CpG24 <- GWASlociSNPs_CpG[which(GWASlociSNPs_CpG$cpg %in% CpG24$CpG),]
GWASlociSNPsPosAlleles <- merge(GWASlociSNPs_CpG24, SNPids, by = "SNPchrPos")
saveRDS(GWASlociSNPsPosAlleles, "GWASlociSNPs_CpGs24_chrPos_region.rds")

CpGloc <- unique(
  data.frame(
     "CpG" = GWASlociSNPsPosAlleles$cpg,
     "chr" = GWASlociSNPsPosAlleles$Chromosome,
     "left" = GWASlociSNPsPosAlleles$start,
     "right"= GWASlociSNPsPosAlleles$end
     )
   )


SNPloc <- unique(
  data.frame(
   "snpid" = GWASlociSNPsPosAlleles$SNPid,
   "chr" = GWASlociSNPsPosAlleles$Chromosome,
   "pos" = GWASlociSNPsPosAlleles$Position
 )
)


## Mvalue residuals file
# resRosmap <- readRDS("C:/ROSMAP_QNBMIQ_PCfiltered_mvalResiduals.RDS")
# sampleKey <- read.table("sampleKey.txt", header = TRUE)


## Filter residual matrix for 688 individuals included in the mqtl analysis and 24 CpGs
# resCpGs <- resRosmap[which(rownames(resRosmap) %in% CpGloc$CpG), which(colnames(resRosmap)%in%sampleKey$methySampleID)]
# sampleKey <- sampleKey[match(Ids, sampleKey$SNPsampleID),]
# identical(colnames(resCpGs), sampleKey$methySampleID)
# colnames(resCpGs) <- sampleKey$SNPsampleID
# saveRDS(resCpGs, "MvalueResiduals_colocCpGs24.rds")
resCpGs <- readRDS("MvalueResiduals_colocCpGs24.rds")

identical(colnames(RosmapPheno), colnames(SNPfile))
identical(colnames(SNPfile), colnames(resCpGs))

### matrixEQTL analysis
library(MatrixEQTL)

snps = SlicedData$new(as.matrix(SNPfile))

CpGs = SlicedData$new(as.matrix(resCpGs))

cvrt = SlicedData$new(as.matrix(RosmapPheno))

outFile = "matrixEQTL_Result_LDregionsCpG24.txt"

me = Matrix_eQTL_main(snps = snps,
                      gene = CpGs,
                      cvrt = cvrt,
                      output_file_name.cis = outFile,
                      pvOutputThreshold = 0,
                      useModel = modelLINEAR,
                      pvOutputThreshold.cis = 1,
                      snpspos = SNPloc,
                      genepos = CpGloc,
                      cisDist = 1157150
)

# > max(abs(GWASlociSNPsPosAlleles$Position-GWASlociSNPsPosAlleles$start))
# [1] 1157141


