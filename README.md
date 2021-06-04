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