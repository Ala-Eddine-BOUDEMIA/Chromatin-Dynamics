# Chromatin-Dynamics
Integrated Analysis of Co-expression Patterns of Histone Variants and their Chaperones Across Tissues from Publicly Available RNA-Seq Data

# 21/06/03

- Most of TCGA attributes are Null or redundunt
- Selected attributes are :

__TCGA__ | __GTEx__
---------|---------
reads_downloaded | run 
mapped_read_count | sample
auc | read_count_as_reported_by_sra
gdc_cases.case_id | reads_downloaded
gdc_cases.demographic.gender | proportion_of_reads_reported_by_sra_downloaded
gdc_cases.demographic.year_of_birth | mapped_read_count	
gdc_cases.demographic.race | auc
gdc_cases.project.primary_site | avg_read_length
gdc_cases.project.project_id | sampid
gdc_cases.tissue_source_site.project | smrin
gdc_cases.diagnoses.tumor_stage | smts
gdc_cases.diagnoses.age_at_diagnosis | smtsd
gdc_cases.diagnoses.vital_status | smntrart
gdc_cases.samples.sample_type | smmaprt
cgc_filename | smexncrt
cgc_file_upload_date | smgnsdtc
cgc_case_year_of_diagnosis | smmncpb
cgc_case_tumor_status | smestlbs	
cgc_case_age_at_diagnosis | smmppd
 / | smnterrt	
 / | smrrnanm	
 / | smrdttl	
 / | smvqcfl	
 / | smtrscpt	
 / | smexpeff
 / | smrrnart

# 21/06/04

# EdgeR

## Library Normalization

EdgeR adjust for :
- Sequencing depth and gene length 
- Library composition 

### Step 1:
- Remove all untranscribed genes (genes with 0 read ciunts in **all samples**)

### Step 2:
- Pick one sample to be the reference that will be used to normalize all of the other samples against :

1. Scale each sample by its total read counts.
2. For each sample, determine the value such that 75% of the scaled data are equal or samller than it.
3. Average the 75th quantiles.
4. The reference sample is the one who's 75th quantile is closest to the average.

### Step 3:
- Select the genes for calculating the scaling factors. This is done for each sample relative to the reference sample.

1. Filter biased genes: 
2. Look at log fold differences of the scaled read counts between the sample and the reference sample (Log of the ratio between the two samples).
3. Remove gene with Infinite values.

4. Identify the genes that are highly and lowly expressed in both samples:
5. Calculate the **geometric mean** for the scaled read counts (Mean of logs).
6. Remove gene with Infinite values.

7. Sort both tables.

8. Filter out the top 30% and the bottom 30% biassed genes.

9. Filter out the top 5% and the bottom 5% of the highly and lowly transcribed genes.

10. Genes that are still in both tables are used to calculate the scaling factor.

### Step 4:
- Calculate the weighted average of the remaining log2 ratios between the two samples. 

### Step 5:
- Convert the weighted average of log2 values to normal numbers.

### Step 6:
- Center the scaling factors around 1.