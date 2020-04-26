# Epigenome-wide meta-analysis of DNA methylation changes in prefrontal cortex highlights the immune processes in Alzheimer’s Disease 
Lanyu Zhang, Tiago Chedraoui Silva, Juan Young, Lissette Gomez, Michael Schmidt, Kara Lyn Hamilton-Nelson, Brian Kunkle, Xi Chen, Eden R. Martin, Lily Wang

### Description

Given the small effect sizes of DNA methylation changes in Alzheimer’s disease (AD) and the inconsistencies often observed in different studies, we conducted a meta-analysis of more than one thousand prefrontal cortex brain samples, to prioritize the most consistent methylation changes in multiple cohorts. Using a uniform analysis pipeline, we identified 3751 CpGs and 119 DMRs significantly associated with Braak stage. Our analysis nominated many new differentially methylated genes such as MAMSTR, AGAP2, AZU1 and provided new insights into epigenetics of AD. For example, the most significant DMR is located on the MAMSTR gene, which encodes a cofactor that stimulates MEF2C. Notably, MEF2C cooperates with another transcription factor PU.1, a central hub in AD gene network. Our enrichment analysis highlighted the particular relevant roles of the immune system and polycomb repressive complex 2 in AD. These results will help facilitate future mechanistic and biomarker discovery studies in AD.

### Single cohort analysis

This section includes scripts for cohort specific analysis. 
The association between CpG methylation levels and Braak stage was assessed using 
linear statistical models for each cohort. We adjusted for potential confounding 
factors including age at death, sex, batch effects, and proportion of different 
cell types in the samples estimated by the CETS R package. 

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
| single_cohort_analysis/Gasparoni.Rmd        |   Gasparoni (Gasparoni, 2018) | [Link to compiled report](https://www.dropbox.com/s/1nfwh6i73rq8836/Gasparoni.html?dl=0)|
| single_cohort_analysis/London.Rmd           |   London (Lunnon, 2014)    | [Link to compiled report](https://www.dropbox.com/s/t56dakqkzxygbhm/London.html?dl=0)|
| single_cohort_analysis/MtSinai.Rmd          |   Mt. Sinai (Smith, 2018)  | [Link to compiled report](https://www.dropbox.com/s/s2r4143u3yp94yq/MtSinai.html?dl=0)|
| single_cohort_analysis/ROSMAP.Rmd           |   ROSMAP (PMID: 29865057)    | [Link to compiled report](https://www.dropbox.com/s/68qrntnh2knskgz/ROSMAP.html?dl=0)|

### Meta-analysis 

To meta-analyze individual CpG results across different cohorts, we used the meta R package. 

For region based meta-analysis, we used two complementary analytical pipelines, 
the comb-p approach and the coMethDMR approach. The significant DMRs were selected by both approaches.  

(1) comb-p appraoch - we used meta-analysis p-values of the four brain samples cohorts as input for comb-p. 

(2) coMethDMR approach - we performed cohort specific analysis for genomic regions first. First, coMethDMR selects co-methylated sub-regions within the contiguous genomic regions. Next, we summarized methylation M values within these co-methylated sub-regions using medians and tested them against AD Braak stage. 

We adjusted for potential confounding factors including age at death, sex, batch, and proportion of different cell types. The cohort specific p-values for each contiguous genomic region were then combined across cohorts using fixed effects meta-analysis model (or random effects model if test of heterogeneity had p-value was less than 0.05). 

Note that the coMethDMR approach allowed us to assess between cohort heterogeneities for genomic regions. 

| File                 | HTML |
|----------------------|----------------------|
| meta_analysis/Meta-analysis.Rmd | [Link to compiled report](https://www.dropbox.com/s/bxmhizaz11tyog7/Meta-analysis.html?dl=0)|


### Enrichment analysis of significant DNA methylation changes 

Meta-analysis results were divided into two groups, methylation changes with positive estimates 
(hypermethylation in AD compared to control) and negative estimae (hypomethylation in AD compared to control). 
For each group, we performed an enrichment analysis (Fisher's test) separately for DMRs and CpGs. 

The main regions/probe annotation used were: 

- Relation_to_Island: `IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC`
- UCSC_RefGene_Group_hierarchy: `IlluminaHumanMethylation450kanno.ilmn12.hg19::Other`
- ChmmModels: `https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E073_15_coreMarks_segments.bed`

| File                 |                      | 
|----------------------|----------------------|
| enrichment_analysis/enrichment-analysis.Rmd | [Link to compiled report](https://www.dropbox.com/s/5adfukotn6gaxfu/enrichment-analysis.html?dl=0)|

### Matched Meta-analysis

To prioritize methylation changes most likely to be affected by the 
AD pathogenesis process, we performed additional analysis using a sample 
matching strategy to reduce confounding effects due to age. 
More specifically, we first matched each case with a control sample using matchControls 
function in e1071 R package. 
The matched samples were then analyzed in the same way as described above, 
except by removing age at death effect in the linear models. 

| File                 | HTML |
|----------------------|----------------------|
| meta_analysis/Matched-analysis.Rmd | [Link to compiled report](https://www.dropbox.com/s/zl0lavf5njx27t1/Matched-analysis.html?dl=0)|

### Correlation of methylation changes in brain and blood samples

Using the London cohort which consisted of 69 pairs of samples with matched PFC and blood samples, 
we compared brain-blood methylation levels in significant CpGs and those CpGs mapped within significant DMRs using Spearman correlations. 
Two approaches were used to quantify methylation levels: using beta values, or using corrected methylation levels (after removing effects of cell type, batch, age at death and sex). In addition, we also conducted look up analysis using the BeCon tool, which compared 
brain-blood methylation levels of Broadmann areas 7, 10 and 20 in postmortem samples of 16 subjects. 

| File                 |    HTML                  | 
|----------------------|----------------------|
| single_cohort_analysis/London_blood.Rmd     | [Link to compiled report](https://www.dropbox.com/s/xgv7err1j9sxkfa/London_blood.html?dl=0)|
| cor_methylation_changes_brain_and_blood/London_blood_brain_correlation.Rmd | [Link to compiled report](https://www.dropbox.com/s/ske24n7nt7lcplw/London_blood_brain_correlation.html?dl=0)|

### Correlation of significant DMRs with expression of nearby genes

We used 529 samples from the ROSMAP study with matched DNA methylation and RNA-seq data for this analysis. 
More specifically, normalized FPKM (Fragments Per Kilobase of transcript per Million mapped reads) 
gene expression values for ROSMAP study were downloaded from AMP-AD Knowledge Portal (Synapse ID: syn3388564). First, we removing potential confounding effects (cell type, batch, sample plate, ageAtDeath, sex) in DNA methylation data and RNAseq data, separately. We then tested for associations between residual methylation levels and residual gene expression levels in genes found ± 250 kb away from the CpG, or start/end of the DMR. 

| File                 |                      | 
|----------------------|----------------------|
| cor_sig_DMRs_exp_nearby_genes/DMR_gene_expression_analysis.Rmd | [Link to compiled report](https://rpubs.com/tiagochst/Supplemental_AD_cor_met_gene_exp)|

### Overlap with genetic susceptibility loci

To evaluate if the significant methylation changes are located closely to genetic risk loci, 
we tested enrichment of significant CpGs and DMRs identified in this study with the 
24 LD blocks of genetic variants reaching genome-wide significance in a 
recent AD meta-analysis102 (PMID: 30820047), using one-sided Fisher’s test. 

| File                 |                      | 
|----------------------|----------------------|
|  ov_genetic_susceptibility_loci/ ov_genetic_susceptibility_loci.Rmd | [Link to compiled report](https://rpubs.com/tiagochst/Supplemental_AD_ov_with_genetic_susc_loc)|



