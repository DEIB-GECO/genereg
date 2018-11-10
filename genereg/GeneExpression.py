# coding: utf-8

# Import libraries
import gmql as gl
import pandas as pd
from pandas import ExcelWriter
import pickle
import collections


def extract_expression(tumor, platform, gencode_version):

	"""
	The EXTRACT_EXPRESSION operation extracts expression values from TCGA for all the genes of interest and their candidate regulatory genes. Intermediate results files are exported locally during the execution of the function, while the final dataframes are returned as Pandas dataframes and exported locally in the Excel files 'Gene Expression - InterestGenes.xlsx' and 'Gene Expression - RegulatoryGenes.xlsx'.

	:param tumor: full name of the tumor of interest, encoded as a string (e.g. 'Ovarian Serous Cystadenocarcinoma', 'Breast Invasive Carcinoma', ...)
	:param platform: number identifying the sequencing platform (either 27 for the 27k probes sequencing platform or 450 for the 450k probes sequencing platform)
	:param gencode_version: number representing the GENCODE genomic annotations to use (currently, for assembly GRCh38, versions 22, 24 and 27 can be used)
	:return: two Pandas dataframes

	Example::
	
		import genereg as gr
		expr_interest_df, expr_regul_df = gr.GeneExpression.extract_expression(tumor='Ovarian Serous Cystadenocarcinoma', platform=27, gencode_version=22)
	"""

	# Check input parameters
	tcga_tumors = ["Acute Myeloid Leukemia","Adrenocortical Carcinoma","Bladder Urothelial Carcinoma","Brain Lower Grade Glioma" ,"Breast Invasive Carcinoma","Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma","Cholangiocarcinoma","Colon Adenocarcinoma","Esophageal Carcinoma","Glioblastoma Multiforme","Head and Neck Squamous Cell Carcinoma","Kidney Chromophobe","Kidney Renal Clear Cell Carcinoma","Kidney Renal Papillary Cell Carcinoma","Liver Hepatocellular Carcinoma","Lung Adenocarcinoma","Lung Squamous Cell Carcinoma","Lymphoid Neoplasm Diffuse Large B-cell Lymphoma","Mesothelioma","Ovarian Serous Cystadenocarcinoma","Pancreatic Adenocarcinoma","Pheochromocytoma and Paraganglioma","Prostate Adenocarcinoma","Rectum Adenocarcinoma","Sarcoma","Skin Cutaneous Melanoma","Stomach Adenocarcinoma","Testicular Germ Cell Tumors","Thymoma","Thyroid Carcinoma","Uterine Carcinosarcoma","Uterine Corpus Endometrial Carcinoma","Uveal Melanoma"]
	if tumor not in tcga_tumors:
		raise ValueError('PATHOLOGY NOT SUPPORTED! You can analyze one of these 33 types of TCGA tumors: '+(', '.join(tcga_tumors)))
	
	if platform not in [27, 450]:
		raise ValueError('PLATFORM NOT RECOGNIZED! Sequencing platforms available: 27 and 450')
	
	if gencode_version not in [22, 24, 27]:
		raise ValueError('GRCh38 GENCODE versions available are 22, 24 and 27')
	
	
	# Load the list of genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})
	
	# Create a list containing the Gene Symbols of the genes of interest
	genesSYM_of_interest = []
	for i, r in EntrezConversion_df.iterrows():
		sym = r['GENE_SYMBOL']
		if sym not in genesSYM_of_interest:
			genesSYM_of_interest.append(sym)

	# Import the dictionary of genes of interest with their candidate regulatory genes
	dict_RegulGenes = pickle.load(open('./2_Regulatory_Genes/dict_RegulGenes.p', 'rb'))

	# Import the gene-TFs mapping dataframe 
	Mapping_df = pd.read_excel('./0_Genes_Mapping/Genes_Mapping.xlsx',sheetname='Sheet1',header=0,converters={'ENTREZ_GENE_ID':str,'HGNC_ID':str})

	# Create a list containing the Gene Symbols of the regulatory genes of genes of interest
	regulatory_genesSYM = []
	for key, value in dict_RegulGenes.items():
		for gene in value:  
			if gene not in regulatory_genesSYM:
				regulatory_genesSYM.append(gene)

	# Extract the list of distinct Gene Symbols mapped in the mapping table
	mapped_gene_SYMs = []
	for index, row in Mapping_df.iterrows():
		sym = row['GENE_SYMBOL']
		if sym not in mapped_gene_SYMs:
			mapped_gene_SYMs.append(sym)


	# Execute the query for the extraction of gene expression values on the remote server, using the PyGMQL Python library
	gl.set_remote_address('http://gmql.eu/gmql-rest/')
	gl.login()
	gl.set_mode('remote')

	# Load the TCGA datasets to be used in the query
	methylation_dataset = gl.load_from_remote(remote_name='GRCh38_TCGA_methylation', owner='public')  
	expression_dataset = gl.load_from_remote(remote_name='GRCh38_TCGA_gene_expression', owner='public') 

	# Identify the sequencing platform to be used
	if platform == 27:
		seq_platform = 'Illumina Human Methylation 27'
	elif platform == 450:
	    seq_platform = 'Illumina Human Methylation 450'
	
	# Extract all the samples for the current tumor and platform
	all_methyl = methylation_dataset.meta_select((methylation_dataset['manually_curated__cases__disease_type'] == tumor) & (methylation_dataset['manually_curated__platform'] == seq_platform) & ((methylation_dataset['biospecimen__bio__sample_type'] == 'Primary Tumor') | (methylation_dataset['biospecimen__bio__sample_type'] == 'Recurrent Tumor')) & (methylation_dataset['clinical__shared__history_of_neoadjuvant_treatment'] == 'No'))
	all_expr = expression_dataset.meta_select((expression_dataset['manually_curated__cases__disease_type'] == tumor) & ((expression_dataset['biospecimen__bio__sample_type'] == 'Primary Tumor') | (expression_dataset['biospecimen__bio__sample_type'] == 'Recurrent Tumor')) & (expression_dataset['clinical__shared__history_of_neoadjuvant_treatment'] == 'No'))

	# Gene Expression:
	expr_0 = all_expr.reg_project(field_list=['ensembl_gene_id','entrez_gene_id','gene_symbol','fpkm'])
	expr = expr_0.meta_select(semiJoinDataset=all_methyl, semiJoinMeta=['biospecimen__bio__bcr_sample_barcode'])

	# Materialize the results into a GDataframe
	expr_Gdf = expr.materialize('./(MaterializeResults)')


	# The result dataset is loaded as a GDataframe, an object containing two pandas dataframes, one for the region data and one for the metadata.
	# Get the two pandas dataframes:
	expr_df_regs = expr_Gdf.regs
	expr_df_meta = expr_Gdf.meta
	n_regs = len(expr_df_regs)
	n_samples = len(expr_df_meta)

	# Rename 'chr', 'start', and 'stop' columns header
	expr_df_regs.rename(columns={'chr':'chrom','start':'left','stop':'right'}, inplace=True)
	# Change index into progressive integer numbers and store the name of the sample in another column
	expr_df_regs['sample_id'] = expr_df_regs.index
	expr_df_regs.index = range(n_regs)

	# Convert unknown values (NaN) to empty strings
	expr_df_regs = expr_df_regs.fillna('')

	# Convert all the metadata values into strings, since they're encode as lists in Python
	col_names = []
	for name, values in expr_df_meta.iteritems():
		col_names.append(name)
	for index, row in expr_df_meta.iterrows():
		for c in col_names:
			list_val = row[c] # it's encoded as a list
			str_val = ''.join(list_val)  # convert the value stored as a list in a string
			expr_df_meta.set_value(index,c,str_val)

		
	# Since we have to extract the expression values for each distinct sample barcode (aliquot), we create a list containing these distinct identifiers
	expr_sample_barcodes_all = []
	for index, row in expr_df_meta.iterrows():
		barcode = row['biospecimen__bio__bcr_sample_barcode']    
		if barcode not in expr_sample_barcodes_all: # get distinct values
			expr_sample_barcodes_all.append(barcode)
        
	# Check which are repeated aliquots, if present
	all_aliqouts = []
	for index, row in expr_df_meta.iterrows():
		barcode = row['biospecimen__bio__bcr_sample_barcode']  
		all_aliqouts.append(barcode)
	multiple_aliquots = [item for item, count in collections.Counter(all_aliqouts).items() if count > 1]

	samples_to_remove = []
	expr_sample_barcodes = []
	if len(multiple_aliquots) != 0:    
		# Among the repeated aliquots, keep only the most recent ones (of 2013)
		for index, row in expr_df_meta.iterrows():
			year = row['biospecimen__bio__year_of_shipment']
			barcode = row['biospecimen__bio__bcr_sample_barcode']  
			if (barcode in multiple_aliquots) and year == '2011':
				expr_df_meta.drop(index, inplace=True)
				samples_to_remove.append(index)

		# Import the list of aliquots in the methylation dataset 
		text_file = open('./3_TCGA_Data/Common_Aliquots.txt', 'r')
		aliquots = text_file.read().split('\n')
		aliquots.remove('')
		text_file.close()
			
		# Extract the new list of distinct TCGA Aliquots to extract
		for index, row in expr_df_meta.iterrows():
			barcode = row['biospecimen__bio__bcr_sample_barcode'] 
			if barcode in aliquots:
				if barcode not in expr_sample_barcodes:
					expr_sample_barcodes.append(barcode)        
			else:
				expr_df_meta.drop(index, inplace=True)
				samples_to_remove.append(index)
			
		# Remove regions that corresponded to eliminated repeated aliquots
		expr_df_regs = expr_df_regs.loc[~(expr_df_regs['sample_id'].isin(samples_to_remove))].copy()

	else:
		expr_sample_barcodes = expr_sample_barcodes_all		

		
	# Export the metadata dataframe setting the TCGA aliquots as indexes.
	Metadata_df = expr_df_meta.copy()
	Metadata_df['id_sample'] = Metadata_df.index
	Metadata_df.set_index('biospecimen__bio__bcr_sample_barcode', inplace=True)
	writer = ExcelWriter('./3_TCGA_Data/Gene_Expression/EXPR_(Metadata).xlsx')
	Metadata_df.to_excel(writer,'Sheet1')
	writer.save()	


	# Extract from the expression dataset all the regions that belong to genes of interest
	expr_df_regs_interest = expr_df_regs.loc[expr_df_regs['gene_symbol'].isin(genesSYM_of_interest)].copy()
	# Extract from the expression dataset all the regions that belong to regulatory genes of genes of interest
	expr_df_regs_regulatory = expr_df_regs.loc[expr_df_regs['gene_symbol'].isin(regulatory_genesSYM)].copy()


	# Gene expression values for each gene of interest:

	# Create a dictionary for storing all the gene expression values for each gene of interest and for each aliquot TCGA
	from collections import defaultdict
	dict_expr_interest = defaultdict(dict)

	for key, value in dict_expr_interest.items():
		value = defaultdict(list)

	# The main dictionary has the Gene Symbols of the genes of interest as keys and each gene has another dictionary as value, which, in turn, has the different aliquots as keys and lists as values.
	# The idea is having a list, containing all the fpkm values, for each gene in each TCGA aliquot.

	# Set the Gene Symbol as keys of the main dictionary
	for name in genesSYM_of_interest:
		dict_expr_interest[name] = {}

	# Set the names of the samples barcodes as keys for each dictionary set as value of a specific key (genes)
	for sample in expr_sample_barcodes:
		for k, v in dict_expr_interest.items():
			v[sample] = []
			
	# Set the values by appending the expression values for each gene of interest: these expression values (fpkm) can be found in the 'expr_df_regs_interest' dataframe
	for index, row in expr_df_regs_interest.iterrows():   # iterating along the whole dataframe
		sym = row['gene_symbol']  # get the Gene Symbol of the gene
		fpkm = row['fpkm']  # get the gene expression value
		sample = row['sample_id']  # get the name of the sample
		# get the aliquot corresponding to current sample
		aliq = expr_df_meta.get_value(sample, 'biospecimen__bio__bcr_sample_barcode')  
		# add the value according to the correct gene ID and TCGA aliquot, rounding it to a float with maximum 6 decimal numbers,
		dict_expr_interest[sym][aliq].append(round(float(fpkm),6))
		

	# Convert the nested dictionary also into a dataframe

	# Create a dataframe whose row indexes are the different TCGA samples and the columns are the distinct genes of interest
	expr_interest_df1 = pd.DataFrame(index = expr_sample_barcodes, columns = [genesSYM_of_interest])

	# Add three additional columns for the name of the sample and the ID and barcode of the patient corresponding to each aliquot, in order to have them available if we will need it
	expr_interest_df2 = pd.DataFrame(index = expr_sample_barcodes, columns = ['Sample_ID','Tumor','Patient_ID'])

	# Create the final dataframe
	expr_interest_df = expr_interest_df1.join(expr_interest_df2)

	# Fill the previously created dataframe with the correct gene expression values, for each gene of interest and for each TCGA aliquot            
	for gene_sym, dict_value in dict_expr_interest.items():
		for tcga_aliq, exp_list in dict_value.items():
			if (len(exp_list) != 0):
				fpkm = exp_list[0]
				# add the expression value in the proper cell of the dataframe, rounding it to a float with maximum 6 decimal numbers
				expr_interest_df.set_value(tcga_aliq,gene_sym,round(fpkm,6))
				

	# Add to the dataframe the name of each sample, the tumor code and the patient's ID in correspondence of each TCGA aliquot
	for index, row in expr_df_meta.iterrows():
		aliquot = row['biospecimen__bio__bcr_sample_barcode']
		tumor_tag = row['clinical__admin__disease_code']
		patient_id = row['clinical__shared__patient_id']
		expr_interest_df.set_value(aliquot,'Sample_ID',index)
		expr_interest_df.set_value(aliquot,'Tumor',tumor_tag)
		expr_interest_df.set_value(aliquot,'Patient_ID',patient_id)
		
	# Add a row at the beginning of the dataframe to insert also the Entrez Gene ID of each gene of interest
	additional_index = ['ENTREZ_GENE_ID']
	expr_interest_df0_1 = pd.DataFrame(index = additional_index, columns = [genesSYM_of_interest])
	expr_interest_df0_2 = pd.DataFrame(index = additional_index, columns = ['Sample_ID','Tumor','Patient_ID'])
	expr_interest_df0 = expr_interest_df0_1.join(expr_interest_df0_2)

	frames = [expr_interest_df0, expr_interest_df]
	expr_interest_df = pd.concat(frames)

	# Add for each Gene Symbol of our genes of interest the corresponding Entrez Gene ID in the first row of the dataframe
	for i, r in EntrezConversion_df.iterrows():
		entrez_id = r['ENTREZ_GENE_ID']
		gene_name = r['GENE_SYMBOL']
		expr_interest_df.set_value('ENTREZ_GENE_ID',gene_name,entrez_id)

	# Set empty strings for NaN values in the 'GENE_SYMBOL' row
	expr_interest_df.set_value('ENTREZ_GENE_ID','Sample_ID',"")
	expr_interest_df.set_value('ENTREZ_GENE_ID','Tumor',"")
	expr_interest_df.set_value('ENTREZ_GENE_ID','Patient_ID',"")


	# Export the dataframe with the gene expression values for our genes of interest for each TCGA aliquot 
	writer = ExcelWriter('./3_TCGA_Data/Gene_Expression/Gene_Expression-InterestGenes.xlsx')
	expr_interest_df.to_excel(writer,'Sheet1')
	writer.save()


	# Gene expression values for each candidate regulatory gene of the genes of interest:

	# Create a dictionary for storing all the gene expression values for each gene of interest and for each aliquot TCGA
	from collections import defaultdict
	dict_expr_regulatory = defaultdict(dict)

	for key, value in dict_expr_regulatory.items():
		value = defaultdict(list)

	# The main dictionary has the Gene Symbols of the candidate regulatory genes as keys and each gene has another dictionary as value, which, in turn, has the different aliquots as keys and lists as values.
	# The idea is having a list, containing all the fpkm values, for each gene in each TCGA aliquot.

	# Set the Gene Symbols as keys of the main dictionary
	for name in regulatory_genesSYM:
		dict_expr_regulatory[name] = {}

	# Set the names of the samples barcodes as keys for each dictionary set as value of a specific key (genes)
	for sample in expr_sample_barcodes:
		for k, v in dict_expr_regulatory.items():
			v[sample] = []
        
	# Set the values by appending the expression values for each candidate regulatory gene: these expression values (fpkm) can be found in the "expr_df_regs_regulatory" dataframe
	for index, row in expr_df_regs_regulatory.iterrows():   # iterating along the whole dataframe
		sym = row['gene_symbol']  # get the Gene Symbol of the gene
		ens_id = row['ensembl_gene_id']  # get the Ensembl Gene ID
		fpkm = row['fpkm']  # get the gene expression value
		sample = row['sample_id']  # get the name of the sample
		# get the aliquot corresponding to current sample
		aliq = expr_df_meta.get_value(sample, 'biospecimen__bio__bcr_sample_barcode')
		# add the value according to the correct gene ID and TCGA aliquot, rounding it to a float with maximum 6 decimal numbers
		if (gencode_version == 22):
			if (ens_id not in ['ENSG00000277726.3','ENSG00000275895.3','ENSGR0000214717.8']):
				dict_expr_regulatory[sym][aliq].append(round(float(fpkm),6))
		else:
			dict_expr_regulatory[sym][aliq].append(round(float(fpkm),6))
	


	# Convert the nested dictionary also into a dataframe

	# Create a dataframe whose row indexes are the different TCGA samples and the columns are the distinct candidate regulatory genes
	expr_regulatory_df1 = pd.DataFrame(index = expr_sample_barcodes, columns = [regulatory_genesSYM])

	# Add three additional columns for the name of the sample and the ID and barcode of the patient corresponding to each aliquot, in order to have them available if we will need it
	expr_regulatory_df2 = pd.DataFrame(index = expr_sample_barcodes, columns = ['Sample_ID','Tumor','Patient_ID'])

	# Create the final dataframe
	expr_regulatory_df = expr_regulatory_df1.join(expr_regulatory_df2)

	# Fill the previously created dataframe with the correct gene expression values, for each candidate regulatory gene and for each TCGA aliquot            
	for gene_sym, dict_value in dict_expr_regulatory.items():
		for tcga_aliq, exp_list in dict_value.items():
			if (len(exp_list) != 0):
				fpkm = exp_list[0]
				# add the expression value in the proper cell of the dataframe, rounding it to a float with maximum 6 decimal numbers
				expr_regulatory_df.set_value(tcga_aliq,gene_sym,round(fpkm,6))
				

	# Add to the dataframe the name of each sample, the tumor code and the patient's ID in correspondence of each TCGA aliquot
	for index, row in expr_df_meta.iterrows():
		aliquot = row['biospecimen__bio__bcr_sample_barcode']
		tumor_tag = row['clinical__admin__disease_code']
		patient_id = row['clinical__shared__patient_id']
		expr_regulatory_df.set_value(aliquot,'Sample_ID',index)
		expr_regulatory_df.set_value(aliquot,'Tumor',tumor_tag)
		expr_regulatory_df.set_value(aliquot,'Patient_ID',patient_id)
		
	# Add a row at the beginning of the dataframe to insert also the Gene Symbols of each gene of interest
	additional_index = ['ENTREZ_GENE_ID']
	expr_regulatory_df0_1 = pd.DataFrame(index = additional_index, columns = [regulatory_genesSYM])
	expr_regulatory_df0_2 = pd.DataFrame(index = additional_index, columns = ['Sample_ID','Tumor','Patient_ID'])
	expr_regulatory_df0 = expr_regulatory_df0_1.join(expr_regulatory_df0_2)

	frames = [expr_regulatory_df0, expr_regulatory_df]
	expr_regulatory_df = pd.concat(frames)

	# Add for each Gene Symbol of the regulatory genes the corresponding Entrez Gene ID in the first row of the dataframe
	for i in regulatory_genesSYM:
		if i == 'PTRF':
			entrez_id = Mapping_df.loc[Mapping_df['GENE_SYMBOL'] == 'CAVIN1', 'ENTREZ_GENE_ID'].iloc[0]
		else:
			entrez_id = Mapping_df.loc[Mapping_df['GENE_SYMBOL'] == i, 'ENTREZ_GENE_ID'].iloc[0]
		expr_regulatory_df.set_value('ENTREZ_GENE_ID',i,entrez_id)

	# Set empty strings for NaN values in the 'GENE_SYMBOL' row
	expr_regulatory_df.set_value('ENTREZ_GENE_ID','Sample_ID',"")
	expr_regulatory_df.set_value('ENTREZ_GENE_ID','Tumor',"")
	expr_regulatory_df.set_value('ENTREZ_GENE_ID','Patient_ID',"")


	# Export the dataframe with the gene expression values for the regulatory genes of our genes of interest for each TCGA aliquot 
	writer = ExcelWriter('./3_TCGA_Data/Gene_Expression/Gene_Expression-RegulatoryGenes.xlsx')
	expr_regulatory_df.to_excel(writer,'Sheet1')
	writer.save()
	
	return expr_interest_df, expr_regulatory_df