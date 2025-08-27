# Cell-and-Tissue-specific-cRE-pipeline

## Set up conda environment 

## Pipeline
Write more detailed about the piepline in Pipeline_Testing Testes.ipynb, writing from a both techincal and boiological perspectives

## Results

## Technology    

Instructions
B. PIPELINE
1. Download FASTQ Files
•	By the time I’m writing this documentation, not all the files are available
•	Once the files are ready to install, please put them under /group/tran3/gchahal/scATAC_analysis_Vento_raw/fastq_path
•	INPUT: Fastq files such as FCAGND8768482_S1_L001_I1_001.fastq.gz
•	OUTPUT: Put the fastq files under /group/tran3/gchahal/scATAC_analysis_Vento_raw/fastq_path
2. Peak Calling (Cellranger ATAC)
•	Each sample like FCAGND8768482 requires 4 FastQ files
•	Perform peak calling with cellranger-atac count
For each sample data (FCAGND8768482, FCA_GND8768484, …)
o	INPUT: 4 FastQ files for this sample data
o	Code: /group/tran3/duytran/PBS_scripts/ scATAC_fastq_pipeline_single_sample1.pbs
o	How to run: qsub scATAC_fastq_pipeline_single_sample1.pbs (Remember to change the –sample argument below to your desired sample)
o	OUTPUT: /group/tran3/gchahal/scATAC_analysis_Vento_raw/FCAGND8768482/ (e.g.)
o	Continue perform peak calling for other samples
•	Note: If there are N samples, then there should be N peak-calling scripts 



3. SIGNAC PIPELINE
o	Perform QC Filtering, UMAP visualizations, Q75 Variants Filtering
o	2 separate pipelines for male and female: 
o	Ovary pipeline: /group/tran3/duytran/cre_pipeline/Gonads_Project_Handover/Pipeline_Testing.ipynb
o	Testes pipeline: tran3/duytran/cre_pipeline/Gonads_Project_Handover/Pipeline_Testing Testes.ipynb

How to run the pipeline (demo on Pipeline_Testing – ovary data)
o	1. Move generated outs folders from step 2 - Peak Calling to their corresponding folder either 
o	/group/tran3/gchahal/Garcia_et_all_FASTqs/male/
o	Or /group/tran3/gchahal/Garcia_et_all_FASTqs/female/ 
o	For example, FCA_GND8046539 is ovary data, then it should be moved to female folder and the path should be: "/group/tran3/gchahal/Garcia_et_all_FASTqs/female/FCAGND8046539/outs/",
o	(INPUT) 2. After moving, open the Jupyter Notebook Pipeline_Testing.ipynb, and ¬edit the first block of notebook with ALL the folder destinations, and their corresponding timepoints. 
 

o	(INPUT) 3. Choose your configuration setup
o	The example below is interpreted as we want to run on ovary data from week 8.6, 8.8, and 9.
o	Note that gender is either ‘female’ (Pipeline_Testing.ipynb) or ‘male’ (Pipeline_Testing Testes.ipynb)
o	The chosen_weeks values (8.6,8.8, 9) should match the values you declared for timepoints earlier in the previous step.
 
o	(OUTPUT) 4. Run the pipeline sequentially all the way from top to down for both notebooks:
o	Visualizations, which could be viewed directly in the Jupyter Notebook
o	BED files for input to Hieu’s Pipeline: /group/tran3/gchahal/other_tissues/CardiacNetworkComponentPredictor/Duy_CardiacNetworkComponentPredictor/data/scATAc_data_other_Domke_cellatlas/q75_filtered_week_8.6_8.8_9_female_bed_file.bed
/group/tran3/gchahal/other_tissues/CardiacNetworkComponentPredictor/Duy_CardiacNetworkComponentPredictor/data/scATAc_data_other_Domke_cellatlas/q75_filtered_week_7_9_male_bed_file.bed
o	BED files for cell-type-specific CREs (with count, used later in Step 5)
/group/tran3/gchahal/other_tissues/CardiacNetworkComponentPredictor/Duy_CardiacNetworkComponentPredictor/data/scATAc_data_other_Domke_cellatlas/scores_q75_filtered_week_8.6_8.8_9_female_bed_file.bed
/group/tran3/gchahal/other_tissues/CardiacNetworkComponentPredictor/Duy_CardiacNetworkComponentPredictor/data/scATAc_data_other_Domke_cellatlas/scores_q75_filtered_week_8.6_8.8_9_female_bed_file.bed
 

4. HIEU CRE PIPELINE (TISSUE-SPECIFIC CRES)
o	INPUT: Organ regions: lifted_heart.bed, lifted_eye.bed, q75_filtered_week_8.6_8.8_9_female_bed_file.bed, …
o	How to run:
o	cd /group/tran3/gchahal/other_tissues/CardiacNetworkComponentPredictor/Duy_CardiacNetworkComponentPredictor/data/scATAc_data_other_Domke_cellatlas/
o	python ../../scripts/LD_getSubsets.py q75_filtered_week_8.6_8.8_9_female_bed_file.bed  q75_filtered_week_7_9_male_bed_file.bed  lifted_heart.bed lifted_eye.bed lifted_cerebrum.bed lifted_stomach.bed ../../out/q75_ovary_testes_from_week_7_to_9/
o	OUTPUT:
o	/group/tran3/gchahal/other_tissues/CardiacNetworkComponentPredictor/Duy_CardiacNetworkComponentPredictor/out/q75_ovary_testes_from_week_7_to_9/
o	This folder contains tissue-specific CREs, which could be submitted to GREAT Tool.
o	The files you want are: testes.bed and ovary.bed
o	 
5. CELL-TYPE SPECIFIC CRES

o	2 separate pipelines for male and female: 
o	Ovary pipeline: /group/tran3/duytran/cre_pipeline/Gonads_Project_Handover/investigation_ovary.ipynb
o	Testes pipeline: tran3/duytran/cre_pipeline/Gonads_Project_Handover/investigation_testes.ipynb
o	The main focus is the row-wise approach, so just run the code inside the row-wise part
 
o	INPUT:
o	Testes-specific and Ovary-specific CREs obtained from Hieu’s pipeline
o	Scores obtained from SIGNAC Pipeline
o	OUTPUT:
o	Ovary folder: /group/tran3/duytran/cre_pipeline/Results/cell_type_specific_regions/row_wise/ovary 
o	Testes folder: /group/tran3/duytran/cre_pipeline/Results/cell_type_specific_regions/row_wise/testes 
o	These regions could be submitted to GREAT tool





BIOLOGICAL BACKGROUND
 
Introduction and Aim:

Disorders of sex development (DSD) are common congenital defects that affect infants worldwide. The underlying cause of DSD is often genetic and thus, requires a comprehensive understanding of the genes involved in gonadal development and disorders. While performing genome sequencing for diagnosing DSD, numerous variants calls are made as “variant of unknown significance” (VUS) 1, 2. These VUS often lie in the non-coding genome or novel genes, where the role in gonadal development remains unexplored. To assign a functional role to these VUS variants, there is an urgent need to decipher the gene regulatory networks involved in gonadal development. Therefore, we aim to identify functional cell-specific regulatory regions and novel genes involved in testis development that when altered result in DSD. 

Rationale:
 
There are several cell lineages involved in development of the testis. We hypothesise that the identification of regulatory region specific for each testis cell type can in turn facilitate the identification of novel genes essential in cell specification and development of the testis. To capture the functional regulatory regions of each cell type, we identify the open chromatin regions using single cell ATAC sequencing (scATAC-seq) experiments. These nucleosome-deprived open chromatin regions represent cis-regulatory regions that are involved in modulating gene expression. 
 
Method: 

We utilise the cis-regulatory directed pipeline 3 to systematically 1) map testis-specific cis-regulatory elements and 2) identify all the genes that are essential for male gonad development and disorders. We develop a male gonad-specific version of the pipeline by first generating a repertoire of cis-regulatory elements essential for male gonad development. The uniquely identified cREs are then used to identify their corresponding gene that are expressed in various cell types of male gonads. These gene are used to construct the gene regulatory networks involved in gonadal development and disorders (Fig.1).
To identify cis-regulatory elements (cREs) that are specific to male gonadal development, we utilise available organ-specific scATAC-seq datasets for testis and other tissues including stomach, heart, eyes and cerebrum (Table 1) 4, 5. 
Table 1: Details of the scATAC datasets utilised for the cRE pipeline.
Organ	Days post conception	Cell type
Testis	84	Coelomic epithelium
 	 	Early somatic cells
 	 	Sertoli
 	 	Leydig cells
 	 	Early supporting gonadal cells
Primordial germ cells
Cerebrum	89	Arythrocytes
 	 	Vascular endothelial cells
 	 	Inhibitory neurons
 	 	Limbic system neurons
 	 	Excitatory neurons
Eye	91	Ganglion cells
 	 	Photoreceptor cells
 	 	Retinal pigment cells
 	 	Ganglion cells
 	 	Photoreceptor cells
 	 	Retinal pigment cells
 	 	Retinal progenitors and Muller glia
Stomach	91	Vascular endothelial cells
 	 	ENS glia cells
 	 	Goblet cells 
 	 	Parietal and chief cells
 	 	Stromal cells 
 	 	Vascular endothelial cells
Heart	94	Cardiomyocytes
 	 	Endocardial cells
 	 	Vascular endothelial cells
 	 	Epicardial fat cells
 	 	Erythroblasts
 	 	Lymphatic endothelial cells
 	 	Lymphoid cells
 	 	Myeloid cells
 	 	Schwann cells
 	 	Smooth muscle cells
 	 	Stromal cells

The open chromatin regions identified in scATAC-seq represent the cis-regulatory elements (cREs) governing the transcriptional regulation in different cell lineages. Comparison of these cREs identified in the developing testis at various cell stages to the cREs from other developing tissues will be performed to identify those that are unique to testicular development (Fig. 1, Step 1 and 2). These unique cREs are essential in understanding the epigenetic mechanisms underlying male gonadal development and its associated disorders. 
The testis-specific cREs will be further used to identify the corresponding genes involved in male gonadal development using the GREAT tool 6 (Fig. 1, Step 3). These genes can be attributed to different cell lineages that are involved in male gonadal development such as coelomic epithelial cells, sertoli cells, leydig cells, early somatic cells etc. Next, we will filter the identified genes with cell-specific gene expression in developing testis using scRNA-seq data for testis development (Fig. 1 Step 4). Genes that are strongly expressed in each cell-type will be retained as crucial components of the gene regulatory network governing that specific cell-type. The identified genes will be then used to reconstruct the cell type-specific gene regulatory networks involved in testis development using STRING database [https://string-db.org/] (Fig.1 Step 5). For validation of the network, we will RNAi screening in mouse models. The gene list will also be compared with the genes identified in developmental sexual disorder patients to confirm the clinical relevance of the identified novel gene markers.


