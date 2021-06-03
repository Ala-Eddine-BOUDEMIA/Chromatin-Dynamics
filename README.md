# Chromatin-Dynamics
Integrated Analysis of Co-expression Patterns of Histone Variants and their Chaperones Across Tissues from Publicly Available RNA-Seq Data

# 21/06/03

- Most of TCGA attributes are Null or redundunt
- Selected attributes are :

__TCGA__ | __GTEx__
---------|---------
Project  | Project
reads_downloaded | sample
mapped_read_count | read_count_as_reported_by_sra
auc | reads_downloaded
gdc_cases.case_id | proportion_of_reads_reported_by_sra_downloaded
gdc_cases.demographic.gender | mapped_read_count
gdc_cases.demographic.year_of_birth | auc	
gdc_cases.demographic.race | avg_read_length
gdc_cases.project.primary_site | sampid
gdc_cases.project.project_id | smrin
gdc_cases.tissue_source_site.project | smts
gdc_cases.diagnoses.tumor_stage | smtsd
gdc_cases.diagnoses.age_at_diagnosis | smntrart
gdc_cases.diagnoses.vital_status | smmaprt
gdc_cases.samples.sample_type | smexncrt
cgc_filename | smgnsdtc
cgc_file_upload_date | smmncpb
cgc_case_year_of_diagnosis | smestlbs
cgc_case_tumor_status | smmppd	
cgc_case_age_at_diagnosis | smnterrt
 / | smrrnanm	
 / | smrdttl	
 / | smvqcfl	
 / | smtrscpt	
 / | smexpeff	
 / | smrrnart