# DNA methylation changes associated with Alzheimer’s Disease pathology: an epigenome-wide meta-analysis of prefrontal cortex brain samples

Authors:  
- Lanyu Zhang
- Tiago Chedraoui Silva
- Juan Young
- Lissette Gomez
- Michael Schmidt
- Kara Lyn Hamilton-Nelson
- Brian Kunkle
- Gabriel J. Odom
- Xi Chen
- Eden R. Martin
- Lily Wang

## Auxliary files

### Single cohort analysis

#### Description

First, we performed cohort specific analysis for individual CpGs. 

The association between CpG methylation levels and Braak stage was assessed using linear statistical models in each cohort. 

Given that methylation M-values (logit transformation of methylation beta values) 
has better statistical properties (i.e. homoscedasticity) for linear regression models, 
we used the M-values as the outcome variable in our statistical models. 

We adjusted for potential confounding factors including 
age at death, sex, methylation slide effects, and proportion of different cell types (i.e. neurons vs. glia cells)
in the samples estimated by the CETS R package. 

For ROSMAP cohort, we additionally included variable "batch" that was available 
in the dataset to adjust for technical batches occurred during data generation.    



Each of the files has the following structure:

1. Data retrieval 
2. Data Pre-processing
    1. Probes QC
    2. Samples QC
3. Outliers detection - PCA analysis
4. Summary after QC steps
5. Compute neuron proportion
6. Linear regression by cpgs Methylation
7. Linear regression by regions median Methylation


#### Files
| File                 | Dataset | HTML |
|----------------------|-------------|-------------| 
| single_cohort_analysis/Gasparoni.Rmd        |   Gasparoni (Gasparoni, 2018) | [Link to compiled report](https://www.dropbox.com/s/1nfwh6i73rq8836/Gasparoni.html?dl=0)|
| single_cohort_analysis/London.Rmd           |   London (Lunnon, 2014)    | [Link to compiled report](https://www.dropbox.com/s/yd74s21mssbo0xq/London.html?dl=0)|
| single_cohort_analysis/London_blood.Rmd     |   London (Lunnon, 2014)     | [Link to compiled report](https://www.dropbox.com/s/yf9vih7dkpdw06r/London_blood.html?dl=0)|
| single_cohort_analysis/MtSinai.Rmd          |   Mt. Sinai (Smith, 2018)  | [Link to compiled report](https://www.dropbox.com/s/tnc12y3myfrx53w/MtSinai.html?dl=0)|
| single_cohort_analysis/ROSMAP.Rmd           |   ROSMAP (PMID: 29865057)    | [Link to compiled report](https://www.dropbox.com/s/8am2p72xlbn0kja/ROSMAP.html?dl=0)|
| single_cohort_analysis/SanchexMut.Rmd       |   SanchexMut (Sanchez-Mut et al. 2016)| [Link to compiled report](https://www.dropbox.com/s/mzrb6vc0c7dmti1/SanchexMut.html?dl=0)|
| single_cohort_analysis/Semick.Rmd           |   Semick (Semick et al. 2019)   | [Link to compiled report](https://www.dropbox.com/s/ubnede70grp5e9e/Semick.html?dl=0)|

### Meta-analysis 

#### Description
To meta-analyze individual CpG results across different cohorts, we used the meta R package. 
The evidence for heterogeneity of study effects was tested using Cochran’s Q statistic. 
Fixed effects model (also referred to as inverse variance-weighted) meta-analysis was applied to synthesize statistical significance from individual cohorts. 
Although the fixed effects model for meta-analysis does not require the assumption of homogeneity, for those regions with nominal evidence for heterogeneity (raw Pheterogeneity < 0.05), 
we also applied random effects meta-analysis and assigned final meta-analysis p-value based on random effects model. 

For region based meta-analysis, we used two complementary analytical pipelines, the comb-p approach and the coMethDMR approach. 
Briefly, comb-p takes single CpG p-values and locations of CpG sites, to scan the genome for regions enriched with series of adjacent low P-values. 
In our analysis, we used meta-analysis p-values of the four brain samples discovery cohorts as input for comb-p. 
As Comb-p uses Sidak method to account for multiple comparisons, we considered DMRs with Sidak p-values less than 0.05 to be significant. 

In the coMethDMR approach, we performed cohort specific analysis for genomic regions first. 
We define “contiguous genomic regions” to be genomic regions on the Illumina array covered with clusters of contiguous 
CpGs where the maximum separation between any two consecutive probes is 200 base pairs. 
Instead of testing all CpGs within a genomic region, coMethDMR carries out an additional step that selects 
co-methylated sub-regions within the contiguous genomic regions first. 
 
Next, we summarized methylation M values within these co-methylated sub-regions using medians and tested them against AD Braak stage. 
In the same way as in single CpG analyses, we adjusted for potential confounding factors including 
age at death, sex, methylation slide, and proportion of different cell types in the samples estimated by the CETS R package. 
 
For ROSMAP cohort, we additionally included variable "batch" that was available in the dataset to adjust for technical batches occurred during data generation.    
To meta-analyze coMethDMR results across different cohorts, first, we assigned co-methylated clusters from each cohort to the non-overlapping contiguous genomic regions that they overlap.
The cohort specific p-values for each contiguous genomic region were then combined across cohorts using fixed effects meta-analysis model (or random effects model if test of heterogeneity had p-value was less than 0.05) as described above. 
Co-methylated DMRs with FDR less than 5% were considered to be significant.

#### Files
| File                 | HTML |
|----------------------|----------------------|
| meta_analysis/Meta-analysis.Rmd | [Link to compiled report](https://www.dropbox.com/s/bxmhizaz11tyog7/Meta-analysis.html?dl=0)|

### Matched Meta-analysis

#### Description
To prioritize methylation changes most likely to be affected by the AD pathogenesis process, 
we performed additional analysis using an alternative strategy to control for confounding effects for age. 
More specifically, we first matched each case with a control sample using matchControls function in e1071 R package. 
The matched samples were analyzed in the same way as described above, except by removing age at death effect in the linear models. 

#### Files

| File                 | HTML |
|----------------------|----------------------|
| meta_analysis/Matched-analysis.Rmd | [Link to compiled report](https://www.dropbox.com/s/1q3srp6g40r9817/Matched-analysis.html?dl=0)|

### Correlation of methylation changes in brain and blood samples


#### Description

Using the London cohort which consisted of 69 samples with matched PFC and blood samples, 
we compared brain-blood methylation levels in significant CpGs and those CpGs mapped within significant DMRs using Spearman correlations. 
Two approaches were used to quantify methylation levels: using beta values, 
or using corrected methylation levels. 
In addition, we also conducted look up analysis using the BeCon tool, which compared brain-blood methylation levels of 
Broadmann areas 7, 10 and 20 in postmortem samples of 16 subjects. 

#### Files
| File                 |    HTML                  | 
|----------------------|----------------------|
| cor_methylation_changes_brain_and_blood/London_blood_brain_correlation.Rmd | [Link to compiled report](https://www.dropbox.com/s/kswu3xa7lzk4g61/London_blood_brain_correlation.html?dl=0)|

### Correlation of significant DMRs with expression of nearby genes

#### Description

The ROSMAP study also generated RNA-seq data for a subset of samples with available DNA methylation data. 
We used 529 samples with matched DNAm and gene expression data for this analysis. 
More specifically, normalized FPKM (Fragments Per Kilobase of transcript per Million mapped reads) 
gene expression values for ROSMAP study were downloaded from AMP-AD Knowledge Portal (Synapse ID: syn3388564). 

Next, for each significant DMR identified in meta-analysis, we first removed confounding effects in DNA methylation data 
by fitting model `median methylation M value ~ neuron.proportions + batch + sample plate + ageAtDeath + sex` 
and extracting residuals from this model, these are the methylation residuals. 

Similarly, we also removed potential confounding effects in RNA-seq data by fitting model
`log2(normalized FPKM values) ~ ageAtDeath + sex  + markers for cell types`. 

The last term “markers for cell types” included multiple covariate variables, to adjust for the multiple types of cells in the brain samples. 

More specifically, we estimated expression levels of genes that are specific 
for the main five cell types present in the CNS, ENO2 for neurons, GFAP for 
astrocytes, CD68 for microglia, OLIG2 for oligodendrocytes and CD34 for endothelial 
cells and included these as variables in the above linear regression model, 
as was done in a previous large study of AD samples. 
The residuals extracted from this model are the gene expression residuals. 
For each gene expression and DMR pair, we then tested the association between 
gene expression residuals and methylation residuals using a linear model 
gene expression residuals ~ methylation residuals for AD cases (Braak stage ≥3) and controls (Braak stage < 3) separately. 
For significant DMS, this analysis was repeated, except by replacing median methylation level
in the DMR with methylation level of the CpG, and correlating with expression values of genes found ± 250 kb away from the CpG. 

To compare the DNA methylation to RNA associations (DNAm-RNA) for DMRs vs. CpGs, 
we used a generalized estimating equations (GEE) model where `-log10Pvalues` from each DMR or CpG were treated as a cluster. 
The GEE model included log p-value of the DNAm-RNA associations as the outcome variable, isDMR (yes/no) and isCaseAssociation (yes/no) as independent variables. 
We assumed an exchangeable working correlation structure for the clusters of correlated observations, a
long with log link and gamma distribution for the outcome variable. 

#### Files
| File                 |                      | 
|----------------------|----------------------|
| cor_sig_DMRs_exp_nearby_genes/DMR_gene_expression_analysis.Rmd | [Link to compiled report](https://www.dropbox.com/s/wo7rn5177g3lgkn/code.html?dl=0)|


### Enrichment analysis of significant DNA methylation changes 

#### Description

Meta-analysis results were slipt into two groups, the ones with positive estimate (hypermethylation compared to control), and the ones with negative estimae (hypomethylation compared to control). 
For each group, we performed an enrichment analysis (fisher test) comparing the significant regions/cpgs (foreground) to all the cpgs/regions used in the single analysis evaluation (background).
The main regions/probe annotation used were: 

- Relation_to_Island: `IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC`
- UCSC_RefGene_Group_hierarchy: `IlluminaHumanMethylation450kanno.ilmn12.hg19::Other`
- ChmmModels: `https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E073_15_coreMarks_segments.bed`

#### Files
| File                 |                      | 
|----------------------|----------------------|
| enrichment_analysis/enrichment-analysis.Rmd | [Link to compiled report](https://www.dropbox.com/s/mco3j8rq8j4as71/enrichment-analysis.html?dl=0)|

