# Epigenome-wide meta-analysis of DNA methylation changes in prefrontal cortex highlights the immune processes in Alzheimer’s Disease 
Lanyu Zhang, Tiago Chedraoui Silva, Juan Young, Lissette Gomez, Michael Schmidt, Kara Lyn Hamilton-Nelson, Brian Kunkle, Xi Chen, Eden R. Martin, Lily Wang     
https://www.nature.com/articles/s41467-020-19791-w

### Description

Given the modest effect sizes of DNA methylation changes in Alzheimer’s disease (AD) and the inconsistencies often observed in different studies, we conducted a meta-analysis of more than one thousand prefrontal cortex brain samples, to prioritize the most consistent methylation changes in multiple cohorts. Using a uniform analysis pipeline, we identified 3751 CpGs and 119 DMRs significantly associated with Braak stage. Our analysis nominated many new differentially methylated genes such as MAMSTR, AGAP2, AZU1 and provided new insights into epigenetics of AD. For example, the most significant DMR is located on the MAMSTR gene, which encodes a cofactor that stimulates MEF2C. Notably, MEF2C cooperates with another transcription factor PU.1, a central hub in AD gene network. Our enrichment analysis highlighted the particular relevant roles of the immune system and polycomb repressive complex 2 in AD. These results will help facilitate future mechanistic and biomarker discovery studies in AD.

### Single cohort analysis

This section includes scripts for cohort specific analysis. 
The association between CpG methylation levels and Braak stage was assessed using 
linear statistical models for each cohort. We adjusted for potential confounding 
factors including age at death, sex, batch effects, and estimated proportions of neurons in the samples. 

Each of the files has the following sections:

1. Data retrieval 
2. Data Pre-processing
    1. Probes QC
    2. Samples QC
3. Outliers detection - PCA analysis
4. Summary after QC steps
5. Compute neuron proportion
6. Linear regression by cpgs Methylation
7. Linear regression by  median Methylation in regions


| File                 | Dataset | HTML |
|----------------------|-------------|-------------| 
| single_cohort_analysis/Gasparoni.Rmd        |   Gasparoni (Gasparoni, 2018) | [Link to compiled report](https://rpubs.com/tiagochst/Supplemental_AD_Gasparoni_dataset)|
| single_cohort_analysis/London.Rmd           |   London (Lunnon, 2014)    | [Link to compiled report](https://rpubs.com/tiagochst/604982)|
| single_cohort_analysis/MtSinai.Rmd          |   Mt. Sinai (Smith, 2018)  | [Link to compiled report](https://rpubs.com/tiagochst/Supplemental_AD_MtSinai_dataset)|
| single_cohort_analysis/ROSMAP.Rmd           |   ROSMAP (PMID: 29865057)    | [Link to compiled report](https://rpubs.com/tiagochst/AD_supplemental_ROSMAP_dataset)|


### Meta-analysis 

To meta-analyze individual CpG results across different cohorts, we used the meta R package. 

For region based meta-analysis, we used two complementary analytical pipelines, 
the comb-p approach and the coMethDMR approach. The significant DMRs were selected by both approaches.  

(1) comb-p appraoch - we used meta-analysis p-values of the four brain samples cohorts as input for comb-p. 

(2) coMethDMR approach - we performed cohort specific analysis for genomic regions first. First, coMethDMR selects co-methylated sub-regions within the contiguous genomic regions. Next, we summarized methylation M values within these co-methylated sub-regions using medians and tested them against AD Braak stage. 

We adjusted for potential confounding factors including age at death, sex, batch, and estimated proportions of neurons. The cohort specific p-values for each contiguous genomic region were then combined across cohorts using fixed effects meta-analysis model (or random effects model if test of heterogeneity p-value was less than 0.05). 

Note that the coMethDMR approach also allowed us to assess between cohort heterogeneities for genomic regions. 

| File                 | HTML |
|----------------------|----------------------|
| meta_analysis/Meta-analysis.Rmd | [Link to compiled report](https://rpubs.com/tiagochst/Supplemental_AD_Meta_analysis)|


### Enrichment analysis of significant DNA methylation changes 

Meta-analysis results were divided into two groups, methylation changes with positive estimates 
(hypermethylation in AD compared to control) and negative estimate (hypomethylation in AD compared to control). 
For each group, we performed an enrichment analysis (Fisher's test) separately for DMRs and CpGs. 

The main regions/probe annotation used were: 

- Relation_to_Island: `IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC`
- UCSC_RefGene_Group_hierarchy: `IlluminaHumanMethylation450kanno.ilmn12.hg19::Other`
- ChmmModels: `https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E073_15_coreMarks_segments.bed`

| File                 |                      | 
|----------------------|----------------------|
| enrichment_analysis/enrichment-analysis.Rmd | [Link to compiled report](https://rpubs.com/tiagochst/Supplemental_AD_Enrichment_analysis)|



### Meta-analysis of age and sex matched samples

To prioritize methylation changes most likely to be affected by the 
AD pathogenesis process, we performed additional analysis using a sample 
matching strategy to reduce confounding effects due to age and sex. 
More specifically, we first matched each case with control samples with the same age at death (in years) and sex in the same cohort using matchControls function in e1071 R package. The age and sex matched samples were then analyzed in the same way as described above, except by removing age at death and sex effects in the linear models. 

| File                 | HTML |
|----------------------|----------------------|
| meta_analysis/6_matched-analysis-by-both-age-and-sex.Rmd | [Link to compiled report](https://rpubs.com/tiagochst/ad_meta_matched_sex_age)|

### Correlation of methylation changes in brain and blood samples

Using the London cohort which consisted of 69 pairs of samples with matched PFC and blood samples, 
we compared brain-blood methylation levels in significant CpGs and those CpGs mapped within significant DMRs using Spearman correlations. 
Two approaches were used to quantify methylation levels: using beta values, or using corrected methylation levels (after removing effects of cell type, batch, age at death and sex). In addition, we also conducted look up analysis using the BeCon tool, which compared 
brain-blood methylation levels of Broadmann areas 7, 10 and 20 in postmortem samples of 16 subjects. 

| File                 |    HTML                  | 
|----------------------|----------------------|
| single_cohort_analysis/London_blood.Rmd     | [Link to compiled report](https://rpubs.com/tiagochst/london_blood)|
| cor_methylation_changes_brain_and_blood/London_blood_brain_correlation.Rmd | [Link to compiled report](https://rpubs.com/tiagochst/Supplemental_AD_london_brain_blood_cor)|

### Correlation of significant DMRs with expression of nearby genes

We used 529 samples from the ROSMAP study with matched DNA methylation and RNA-seq data for this analysis. 
More specifically, normalized FPKM (Fragments Per Kilobase of transcript per Million mapped reads) 
gene expression values for ROSMAP study were downloaded from AMP-AD Knowledge Portal (Synapse ID: syn3388564). First, we removing potential confounding effects (cell type, batch, sample plate, ageAtDeath, sex) in DNA methylation data and RNAseq data, separately. We then tested for associations between residual methylation levels and residual gene expression levels in genes found ± 250 kb away from the CpG, or start/end of the DMR. 

| File                 |                      | 
|----------------------|----------------------|
| cor_sig_DMRs_exp_nearby_genes/gene-expression-analysis-for-all-not-by-case-control.Rmd| [Link to compiled report](https://rpubs.com/tiagochst/cor_ad_met_gene)|

### Correlation and co-localization with genetic susceptibility loci

To identify methylation quantitative trait loci (mQTLs) for the significant DMRs and CpGs, we tested associations between the methylation levels with nearby SNPs, using the ROSMAP study dataset, which had matched genotype data and DNA methylation data for 688 samples. To evaluate if the significant methylation changes are located closely to genetic risk loci, we tested enrichment of significant CpGs and DMRs identified in this study with the 24 LD blocks of genetic variants reaching genome-wide significance in a recent AD meta-analysis (PMID: 30820047). In addition, we also performed a co-localization analysis using the method of Giambartolomei et al. (2014) (PMID: 24830394), to determine whether the association signals at the GWAS loci (variant to AD status, and variant to CpG methylation levels) are due to a single shared causal variant or to distinct causal variants close to each other. 

| File                 |                      | 
|----------------------|----------------------|
| ov_genetic_susceptibility_loci/ov_genetic_susceptibility_loci.Rmd | [Link to compiled report](https://rpubs.com/tiagochst/Supplemental_AD_ov_with_genetic_susc_loc)|
| ov_genetic_susceptibility_loci/ADcogdx_mqtlAnalysiPCbatch_CpGs_Code.R  | [Link to script](https://github.com/TransBioInfoLab/ad-meta-analysis/blob/master/ov_genetic_susceptibility_loci/ADcogdx_mqtlAnalysiPCbatch_CpGs_Code.R)|  
| ov_genetic_susceptibility_loci/ADcogdx_mqtlAnalysisPCbatch_DMRs_Code.R |[Link to script](https://github.com/TransBioInfoLab/ad-meta-analysis/blob/master/ov_genetic_susceptibility_loci/ADcogdx_mqtlAnalysisPCbatch_DMRs_Code.R)|
| ov_genetic_susceptibility_loci/colocAnalysis24CpGs_code.R |[Link to script](https://github.com/TransBioInfoLab/ad-meta-analysis/blob/master/ov_genetic_susceptibility_loci/colocAnalysis24CpGs_code.R)|
| ov_genetic_susceptibility_loci/matrixEQTLcodeForColocalisation_24CpG.R  |[Link to script](https://github.com/TransBioInfoLab/ad-meta-analysis/blob/master/ov_genetic_susceptibility_loci/matrixEQTLcodeForColocalisation_24CpG.R)|

### Sensitivity analysis

To assess the potential inflation in our results, we estimated genomic inflation factors using both the conventional and the _bacon_ method (PMID: 28129774), specifically proposed for EWAS. In addition, we also conducted sensitivity analyses to confirm main findings from our enrichment analysis, using inflation corrected effect sizes computed using the _bacon_ method.  

This section includes scripts for:
1. Estimation of genomic inflation factors using both the conventional and _bacon_ method for results previously described in the "Single cohort analysis" section. 
2. Meta-analyzing inflation corrected individual CpG results across different cohorts.
3. Enrichment analysis of significant CpGs obtained in step 2., using the same approach as previously described in section "Enrichment analysis of significant DNA methylation changes" above.
4. Enrichment analysis using a logistic mixed effects regression model that accounts for correlations between CpGs in the same chromosome (file: hpmixed_glimmix.sas).

| File                 | HTML |
|----------------------|-------------| 
| sensitivity_analysis/single-cohort-analysis.Rmd  | [Link to compiled report](https://rpubs.com/tiagochst/Sensitivity_Analysis_estimation_genomic_inflation)  |  
| sensitivity_analysis/Meta-analysis.Rmd  | [Link to compiled report](https://rpubs.com/tiagochst/Meta_analysis_using_bacon_inflation)  |  
| sensitivity_analysis/enrichment-analysis-bacon.Rmd  | [Link to compiled report](https://rpubs.com/tiagochst/enrichment_analysis_of_sigCpGs_after_inflation_correction_bacon) |
| sensitivity_analysis/hpmixed_glimmix  |[Link to script](https://github.com/TransBioInfoLab/ad-meta-analysis/blob/master/sensitivity_analysis/hpmixed_glimmix.sas)|

### Genome-wide summary statistics

Analysis results after bacon correction [Link to file](https://www.dropbox.com/s/ntdmobyl1w1hvvy/meta_analysis_single_cpg_bacon_df.csv?st=7rz5r03u&dl=0)

### Acknowledgement
All datasets used in this study are publicly available. The Mt. Sinai, London, Gasparoni and ROSMAP datasets were obtained from GEO (accessions GSE80970, GSE59685, GSE66351) and Synapse (accession syn3157275). The ROSMAP study data were provided by the Rush Alzheimer’s Disease Center, Rush University Medical Center, Chicago. Data collection was supported through funding by NIA grants P30AG10161, R01AG15819, R01AG17917, R01AG30146, R01AG36836, U01AG32984, U01AG46152, the Illinois Department of Public Health, and the Translational Genomics Research Institute. 


