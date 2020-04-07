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

#### Descritption

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


| File                 | Dataset |
|----------------------|-------------|
| Gasparoni.Rmd        |   Gasparoni (Gasparoni, 2018) |
| London.Rmd           |   London (Lunnon, 2014)    |
| London_blood.Rmd     |   London (Lunnon, 2014)     |
| MtSinai.Rmd          |   Mt. Sinai (Smith, 2018)  |
| ROSMAP.Rmd           |   ROSMAP (PMID: 29865057)    |
| SanchexMut.Rmd       |   SanchexMut (Sanchez-Mut et al. 2016)|
| Semick.Rmd           |   Semick (Semick et al. 2019)   |

### Meta-analysis 

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


| File                 | Description |
|----------------------|-------------|
| Meta-analysis.Rmd    |             |


### Matched Meta-analysis

To prioritize methylation changes most likely to be affected by the AD pathogenesis process, 
we performed additional analysis using an alternative strategy to control for confounding effects for age. 
More specifically, we first matched each case with a control sample using matchControls function in e1071 R package. 
The matched samples were analyzed in the same way as described above, except by removing age at death effect in the linear models. 


| File                 | Description |
|----------------------|-------------|
| Matched-analysis.Rmd |             |

### Correlation of methylation changes in brain and blood samples

| File                 | Description |
|----------------------|-------------|

