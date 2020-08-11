## Colocalisation analysis, Rosmap dataset (n=688)
## mQTLs (24 AD CpGs - SNPs in AD GWAS locus in cis with CpG) / GWAS results 
## Lissette Gomez, July 2020


setwd("M:/")

GWASresults <- read.table("Kunkle_etal_Stage1_results.txt", header = TRUE)
mQTLresults <- read.table("matrixEQTL_Result_LDregionsCpG24.txt", header = TRUE)
resCpGs <- readRDS("MvalueResiduals_colocCpGs24.rds")

library(coloc)
GWASresults$chrPos <- paste0(GWASresults$Chromosome, ":", GWASresults$Position)

spl <-strsplit(as.character(mQTLresults$SNP), "_")
mQTLresults$A1 <- sapply(spl, "[", 2)
mQTLresults$A2 <- sapply(spl, "[", 3)
mQTLresults$chrPos <- sapply(spl, "[", 1)

SNPoverlap <- merge(mQTLresults,
                    GWASresults,
                    by = "chrPos")
SNPs <- unique(SNPoverlap$chrPos) #n=5358

#Subset GWAS results for SNPs in both datasets
GWASresultsOverlap <- GWASresults[which(GWASresults$chrPos %in% SNPs),]

#Find the SNPs with betas for the same allele in both datasets
effectAlleleOverlap <- merge(mQTLresults,
                             GWASresultsOverlap,
                             by.x = c("chrPos", "A1"),
                             by.y = c("chrPos", "Effect_allele"))

overlapChrPos <- unique(effectAlleleOverlap$chrPos) #n=3785

#Find SNPs with different effect allele in GWAS and mQTL results
diffEffectAllesChrPos <- SNPs[which(!SNPs%in%overlapChrPos)] #n=4375

#Change Beta and alleles for SNPs with different effect allele in each dataset
noOverlapGWASres <- GWASresultsOverlap[which(GWASresultsOverlap$chrPos %in% diffEffectAllesChrPos),]
noOverlapGWASres$Beta <- -noOverlapGWASres$Beta
nonEffectAlleleFix <- noOverlapGWASres$Effect_allele
EffectAlleleFix <- noOverlapGWASres$Non_Effect_allele 
noOverlapGWASres$Effect_allele <- EffectAlleleFix
noOverlapGWASres$Non_Effect_allele <- nonEffectAlleleFix

#Concat SNPs with same and different effect allele after changing the beta 
GWASresults_fixEffectAllele <- rbind(
  GWASresultsOverlap[which(GWASresultsOverlap$chrPos %in% overlapChrPos),
                     c("chrPos", "Effect_allele", "Non_Effect_allele", "Beta", "SE")],
  noOverlapGWASres[, c("chrPos", "Effect_allele", "Non_Effect_allele", "Beta", "SE")]
)

GWASresults_fixEffectAllele$SNP <- paste0(GWASresults_fixEffectAllele$chrPos,
                                          "_",
                                          GWASresults_fixEffectAllele$Effect_allele,
                                          "_",
                                          GWASresults_fixEffectAllele$Non_Effect_allele)



##Run co-localisation test by CpG
CpG <- unique(mQTLresults$gene)

dataset2<-list(beta= GWASresults_fixEffectAllele$Beta,
               varbeta= (GWASresults_fixEffectAllele$SE^2),
               type = "cc",
               snp = GWASresults_fixEffectAllele$SNP,
               s=0.34,
               N=63926)


summaryRes <- data.frame(matrix(ncol=7,nrow=0))
allRes <- data.frame(matrix(ncol=12,nrow=0))
for (i in CpG[1:length(CpG)]){
  
  mQTLcpg <- mQTLresults[which(mQTLresults$gen == i),] 
  
  methyResCpG <-as.numeric(resCpGs[i,])
  
  dataset1<-list(beta=mQTLcpg$beta,
                 varbeta=(mQTLcpg$beta/mQTLcpg$t.stat)^2,
                 type = "quant",
                 snp = mQTLcpg$SNP,
                 sdY=sd(methyResCpG),
                 N=688)
  
  my.res<-coloc.abf(dataset1, dataset2)
  
  allRes <- rbind(allRes, cbind(i,my.res$results))
  summaryRes <- rbind(summaryRes, cbind(i, as.data.frame(t(my.res$summary))))
  
  
}

colnames(allRes)[1] <- "CpG"
colnames(summaryRes)[1] <- "CpG"

write.csv(allRes, "colocalisationResults_all_LDregionCpG24.csv", row.names = FALSE)
write.csv(summaryRes, "colocalisationResults_summary_LDregionCpG24.csv", row.names = FALSE)

annot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19::Other)
write.csv(annot, "IlluminaMethy450kannot.csv")
