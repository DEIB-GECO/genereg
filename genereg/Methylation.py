# coding: utf-8

# Import libraries
import gmql as gl
import pandas as pd
from pandas import ExcelWriter
import pickle
import numpy as np
import math
import collections
from collections import defaultdict
from sklearn.utils import shuffle


def extract_methylation(tumor, platform, gencode_version, methyl_upstream, methyl_downstream):

	"""
	The EXTRACT_METHYLATION operation extracts methylation values from TCGA for all the genes of interest. For each gene of interest, the mean value of all the beta_values associated to methylation sites that are localized within areas -methyl_upstream/+methyl_downstream bases from its TSSs are retrieved. Intermediate results files are exported locally during the execution of the function, while the final dataframe is returned as a Pandas dataframe and exported locally in the Excel file 'Methylation Values.xlsx'.

	:param tumor: full name of the tumor of interest, encoded as a string (e.g. 'Ovarian Serous Cystadenocarcinoma', 'Breast Invasive Carcinoma', ...)
	:param platform: number identifying the sequencing platform (either 27 for the 27k probes sequencing platform or 450 for the 450k probes sequencing platform)
	:param gencode_version: number representing the GENCODE genomic annotations to use (currently, for assembly GRCh38, versions 22, 24 and 27 can be used)
	:param methyl_upstream: number of bases upstream the gene TSS to consider for the extraction of methylation sites of interest
	:param methyl_downstream: number of bases downstream the gene TSS to consider for the extraction of methylation sites of interest
	:return: a Pandas dataframe
	
	Example::
	
		import genereg as gr
		methyl_df = gr.Methylation.extract_methylation(tumor='Ovarian Serous Cystadenocarcinoma', platform=27, gencode_version=22, methyl_upstream=4000, methyl_downstream=1000)
	"""

	# Check input parameters
	tcga_tumors = ["Acute Myeloid Leukemia","Adrenocortical Carcinoma","Bladder Urothelial Carcinoma","Brain Lower Grade Glioma" ,"Breast Invasive Carcinoma","Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma","Cholangiocarcinoma","Colon Adenocarcinoma","Esophageal Carcinoma","Glioblastoma Multiforme","Head and Neck Squamous Cell Carcinoma","Kidney Chromophobe","Kidney Renal Clear Cell Carcinoma","Kidney Renal Papillary Cell Carcinoma","Liver Hepatocellular Carcinoma","Lung Adenocarcinoma","Lung Squamous Cell Carcinoma","Lymphoid Neoplasm Diffuse Large B-cell Lymphoma","Mesothelioma","Ovarian Serous Cystadenocarcinoma","Pancreatic Adenocarcinoma","Pheochromocytoma and Paraganglioma","Prostate Adenocarcinoma","Rectum Adenocarcinoma","Sarcoma","Skin Cutaneous Melanoma","Stomach Adenocarcinoma","Testicular Germ Cell Tumors","Thymoma","Thyroid Carcinoma","Uterine Carcinosarcoma","Uterine Corpus Endometrial Carcinoma","Uveal Melanoma"]
	if tumor not in tcga_tumors:
		raise ValueError('PATHOLOGY NOT SUPPORTED! You can analyze one of these 33 types of TCGA tumors: '+(', '.join(tcga_tumors)))
	
	if platform not in [27, 450]:
		raise ValueError('PLATFORM NOT RECOGNIZED! Sequencing platforms available: 27 and 450')
	
	if gencode_version not in [22, 24, 27]:
		raise ValueError('GRCh38 GENCODE versions available are 22, 24 and 27')

	# Execute the query for the extraction of methylation values on the remote server, using the PyGMQL Python library
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

	# Methylation:
	methyl_0 = all_methyl.reg_project(field_list=['beta_value'])
	methyl = methyl_0.meta_select(semiJoinDataset=all_expr, semiJoinMeta=['biospecimen__bio__bcr_sample_barcode'])

	# Materialize the results into a GDataframe
	methyl_Gdf = methyl.materialize('./(MaterializeResults)')


	# The result dataset is loaded as a GDataframe, an object containing two pandas dataframes, one for the region data and one for the metadata.
	# Get the two pandas dataframes:
	methyl_df_regs = methyl_Gdf.regs
	methyl_df_meta = methyl_Gdf.meta
	n_regs = len(methyl_df_regs)
	n_samples = len(methyl_df_meta)

	# Change index into progressive integer numbers and store the name of the sample in another column
	methyl_df_regs['sample_id'] = methyl_df_regs.index
	methyl_df_regs.index = range(n_regs)

	# Convert all the metadata values into strings, since they're encode as lists in Python
	col_names = []
	for name, values in methyl_df_meta.iteritems():
		col_names.append(name)
	for index, row in methyl_df_meta.iterrows():
		for c in col_names:
			list_val = row[c] # it's encoded as a list
			str_val = ''.join(list_val)  # convert the value stored as a list in a string
			methyl_df_meta.set_value(index,c,str_val)
		
	# Export the metadata dataframe setting the TCGA aliquots as indexes.
	Metadata_df = methyl_df_meta.copy()
	Metadata_df['id_sample'] = Metadata_df.index
	Metadata_df.set_index('biospecimen__bio__bcr_sample_barcode', inplace=True)
	writer = ExcelWriter('./3_TCGA_Data/Methylation/METHYL (Metadata).xlsx')
	Metadata_df.to_excel(writer,'Sheet1')
	writer.save()

	# Extract the sample barcodes (TCGA Aliquots)
	methyl_sample_barcodes = []
	for index, row in methyl_df_meta.iterrows():
		barcode = row['biospecimen__bio__bcr_sample_barcode']    
		if barcode not in methyl_sample_barcodes: # get distinct values
			methyl_sample_barcodes.append(barcode)

	# Load the list of genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})

	# Create a list containing the Gene Symbols of the genes of interest
	genesSYM_of_interest = []
	for i, r in EntrezConversion_df.iterrows():
		name = r['GENE_SYMBOL']
		if name not in genesSYM_of_interest:
			genesSYM_of_interest.append(name)

	
	# Create a dictionary for storing all the methylation values for each gene of interest and for each aliquot TCGA
	dict_methyl_list = defaultdict(dict)

	for key, value in dict_methyl_list.items():
		value = defaultdict(list)
	# The main dictionary has the Gene Symbols of the genes of interest as keys and each gene has another dictionary as value, which, in turn, has the different aliquots as keys and lists as values.
	# The idea is having a list, containing all the beta_values, for each gene in each TCGA aliquot.

	# Set the Gene Symbols as keys of the main dictionary
	for name in genesSYM_of_interest:
		dict_methyl_list[name] = {}

	# Set the samples barcodes as keys for each dictionary set as value of a specific key (genes)
	for sample in methyl_sample_barcodes:
		for k, v in dict_methyl_list.items():
			v[sample] = []


	# Extract the methyl_areas dataset
	gl.set_remote_address('http://gmql.eu/gmql-rest/')
	gl.login()
	gl.set_mode('remote')

	ann_dataset = gl.load_from_remote(remote_name='GRCh38_ANNOTATION_GENCODE', owner='public')
	annotations_version = str(gencode_version)
	coding_transcripts_0 = ann_dataset.meta_select((ann_dataset['release_version'] == annotations_version) & (ann_dataset['annotation_type'] == 'transcript'))
	coding_transcripts = coding_transcripts_0.reg_select((coding_transcripts_0.transcript_type == 'protein_coding') & ((coding_transcripts_0.tag == 'basic') | (coding_transcripts_0.tag == 'CCDS')))
	methyl_areas_reg = coding_transcripts.reg_project(['gene_id','gene_name','entrez_gene_id'], new_field_dict = {'start': coding_transcripts.start - methyl_upstream, 'stop': coding_transcripts.start + methyl_downstream})
	gencode_grch38_methyl_areas = methyl_areas_reg.group(regs=['gene_name'], regs_aggregates = {'ensembl_gene_id': gl.BAG('gene_id'), 'gene_symbol': gl.BAG('gene_name'), 'entrez_gene_id': gl.BAG('entrez_gene_id')})

	# Materialize the results into a GDataframe
	Gencode_df_TSS_Gdf = gencode_grch38_methyl_areas.materialize('./(MaterializeResults)')

	# Get the regions dataframe
	Gencode_df_TSS = Gencode_df_TSS_Gdf.regs
	Gencode_df_TSS.rename(columns={'chr': 'chrom', 'start': 'methyl_left', 'stop': 'methyl_right'}, inplace=True)

	# Remove the transcripts that don't belong to genes of interest
	Gencode_df_TSS_interest = Gencode_df_TSS.loc[Gencode_df_TSS['gene_symbol'].isin(genesSYM_of_interest)].copy()

	# Extract from the TCGA data only useful columns for the following procedure
	methyl_df_regs.rename(columns={'chr': 'chrom', 'start': 'left', 'stop': 'right'}, inplace=True)
	methyl_df_regs_useful = methyl_df_regs[['chrom','left','right','strand','beta_value','sample_id']].copy()

	# Create a dictionary for storing all the methylation regions associated to each gene of interest
	dict_methyl_df = {}

	# Set the Gene Symbols of genes of interest as keys of the main dictionary and an empty dataframe as values (with the same columns of 'methyl_df_regs_useful')
	columns = ['left','right','strand','beta_value','sample_id']
	for i in genesSYM_of_interest:
		dict_methyl_df[i] = pd.DataFrame(columns=columns)

	# The dictionary has the Gene Symbols of the genes of interest as keys and a dataframe containing all the methylation regions with genomic coordinates that are within +methyl_upstream/-methyl_downstream bases from the TSS, for each gene of interest.
	# Fill the empty dataframes set as values in the dictionary.

	# Iterate along the GENCODE dataframe containing transcripts belonging to genes of interest
	for index, row in Gencode_df_TSS_interest.iterrows():
		# extract values of attributes we are interested in
		gencode_chrom = row['chrom']
		gencode_left = row['methyl_left']
		gencode_right = row['methyl_right']
		gene = row['gene_symbol']
		# create a list with 'int' elements in the range [gencode_left, gencode_right)
		methyl_area = list(range(gencode_left,gencode_right))
		# select the methylation regions that are inside the region selected (i.e. 'methyl_area')
		selected_methyl_regs = methyl_df_regs_useful.loc[(methyl_df_regs_useful['chrom'] == gencode_chrom) & (methyl_df_regs_useful['left'].isin(methyl_area)) & (methyl_df_regs_useful['right'].isin(methyl_area))].copy()
		# set the extracted dataframe as value of the corresponding key (gene) in the dictionary
		value_df = dict_methyl_df[gene]  # get the old dataframe
		# concatenate the old dataframe and the new one as value in the dictionary
		frames = [value_df, selected_methyl_regs]
		dict_methyl_df[gene] = pd.concat(frames)
		
	# For each dataframe set as value in the dictionary, remove duplicated rows, if present
	for key, value in dict_methyl_df.items():
		value.drop_duplicates(keep='first', inplace=True)
		
	# Store in a list the Entrez Gene IDs of the genes of interest for which no regions has been found
	gencode_missing_values_genes = []
	for key, value in dict_methyl_df.items():
		if len(value) == 0:
			gencode_missing_values_genes.append(key)

	# Extract the methylation values for each gene of interest and for each TCGA aliquot.
	# Set the values by appending the methylation values for each gene of interest: these methylation values (beta_values) can be found in the dataframes set as values in dictionary "dict_ov_methyl_df".
	for gene, value_df in dict_methyl_df.items():
		for index, row in value_df.iterrows():
			beta = row['beta_value']  # get the methylation value
			sample = row['sample_id']  # get the name of the sample
			# get the aliquot corresponding to current sample
			aliq = methyl_df_meta.get_value(sample, 'biospecimen__bio__bcr_sample_barcode')
			# add the value according to the correct Gene Symbol and TCGA aliquot
			dict_methyl_list[gene][aliq].append(round(float(beta),8))

	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------		
	# Extract in a list the names of the distinct sample barcodes (aliquots)
	methyl_sample_barcodes = list((list(dict_methyl_list.values()))[0].keys())
	# Export the list of common aliquots in a .txt file
	with open ('./3_TCGA_Data/Common_Aliquots.txt', 'w') as fp:
		for i in methyl_sample_barcodes:
			fp.write("%s\n" % i)
			
	# Shuffle and randomly splits into five different sets the TCGA aliquots to be analyzed.
	# This five sets of aliquots will be used in the feature selection procedure during the data analysis phase as five different test sets (with the remaining aliquots forming the corresponding training set), allowing the data analysis method to be trained and tested.
	# Thus, in order to reduce the bias, a cross-validation procedure is adopted and the feature selection is executed on each data matrix five times: the final set of features selected for that matrix 
	# is the intersection of the five different subsets returned by the five different feature selection sub-processes on that same data matrix.

	# Import the list of common TCGA aliquots to analyze
	aliquot_file = open('./3_TCGA_Data/Common_Aliquots.txt', 'r')
	aliquots = aliquot_file.read().split('\n')
	aliquots.remove('')
	aliquot_file.close()

	# Create a dataframe having the TCGA aliquots as indexes of its rows
	model_gene_df = pd.DataFrame(index=aliquots, columns=['C1','C2'])

	# Shuffle the rows of the model gene dataframe in a random way, in order to reduce the bias
	model_gene_df = shuffle(model_gene_df)

	# Split the dataframe into five dataframes that will be used as test sets
	model_gene_df_split = np.array_split(model_gene_df, 5)
	model_gene_df_test1 = model_gene_df_split[0]
	model_gene_df_test2 = model_gene_df_split[1]
	model_gene_df_test3 = model_gene_df_split[2]
	model_gene_df_test4 = model_gene_df_split[3]
	model_gene_df_test5 = model_gene_df_split[4]

	# Save the aliquots selected for each of the five test dataframes in a dictionary
	dict_test_split = defaultdict(list)
	dict_test_split['Test_1'] = list(model_gene_df_test1.index.values)
	dict_test_split['Test_2'] = list(model_gene_df_test2.index.values)
	dict_test_split['Test_3'] = list(model_gene_df_test3.index.values)
	dict_test_split['Test_4'] = list(model_gene_df_test4.index.values)
	dict_test_split['Test_5'] = list(model_gene_df_test5.index.values)

	# Export the dictionary 
	pickle.dump(dict_test_split, open('./5_Data_Analysis/dict_test_split.p', 'wb'))
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	
	# Convert the nested dictionary into a dataframe:

	# Create a dataframe whose row indexes are the different TCGA samples and the columns are the distinct genes of interest
	methyl_df1 = pd.DataFrame(index = methyl_sample_barcodes, columns = [genesSYM_of_interest])

	# Add three additional columns for the name of the sample, the ID of the patient and the tumor tag corresponding to each aliquot
	methyl_df2 = pd.DataFrame(index = methyl_sample_barcodes, columns = ['Sample_ID','Tumor','Patient_ID'])

	# Create the final dataframe
	methyl_list_df = methyl_df1.join(methyl_df2)

	# Add to the dataframe the name of each sample, the patient ID and the tumor tag in correspondence of each TCGA aliquot
	for index, row in Metadata_df.iterrows():
		sample = row['id_sample']
		tumor_tag = row['clinical__admin__disease_code']
		patient_id = row['clinical__shared__patient_id']
		methyl_list_df.set_value(index,'Sample_ID',sample)
		methyl_list_df.set_value(index,'Tumor',tumor_tag)
		methyl_list_df.set_value(index,'Patient_ID',patient_id)
    
	# Add a row at the beginning of the dataframe to insert also the Entrez Gene ID of each gene of interest
	additional_index = ['ENTREZ_GENE_ID']
	methyl_df0_1 = pd.DataFrame(index = additional_index, columns = [genesSYM_of_interest])
	methyl_df0_2 = pd.DataFrame(index = additional_index, columns = ['Sample_ID','Tumor','Patient_ID'])
	methyl_df0 = methyl_df0_1.join(methyl_df0_2)

	frames = [methyl_df0, methyl_list_df]
	methyl_list_df = pd.concat(frames)

	# Add for each Gene Symbol of our genes of interest the corresponding Entrez ID in the first row of the dataframe
	for i, r in EntrezConversion_df.iterrows():
		entrez_id = r['ENTREZ_GENE_ID']
		gene_name = r['GENE_SYMBOL']
		methyl_list_df.set_value('ENTREZ_GENE_ID',gene_name,entrez_id)

	# Set empty strings for NaN values in the 'GENE_SYMBOL' row
	methyl_list_df.set_value('ENTREZ_GENE_ID','Sample_ID',"")
	methyl_list_df.set_value('ENTREZ_GENE_ID','Tumor',"")
	methyl_list_df.set_value('ENTREZ_GENE_ID','Patient_ID',"")

	# Add to the dataframe the list of methylation values for each gene of interest in each aliquot TCGA
	for gene, dict_value in dict_methyl_list.items():
		for tcga_aliq, beta_list in dict_value.items():
			# get the list of beta_values for gene 'gene' and aliquot 'tcga_aliq' and add it in the proper cell of the dataframe
			methyl_list_df.set_value(tcga_aliq,gene,beta_list)


	# Compute the MEAN for the beta_values:

	# In case the same gene has more than one beta_value for a single sample, compute their median value and set it as the new beta_value.
	# In this way, we will have a single methylation value for each gene in each sample.
	dict_methyl = dict_methyl_list.copy()
	methyl_df = methyl_list_df.copy()

	sum_values = 0
	count_values = 0

	for gene_name, dict_value in dict_methyl.items():
		for tcga_aliq, beta_list in dict_value.items():
			# get the list of beta_values for gene 'entrez_id' and aliquot 'tcga_aliq'
			for v in beta_list:
				if (len(beta_list) != 0):  # if the list of beta_values is not empty
					if not(math.isnan(v)): # if the values considered is not 'nan'
						# consider the current value
						sum_values += v
						count_values += 1
			# if there's more than one beta_value for the same gene in the same sample
			if (count_values != 0):
				# compute the median value
				single_beta_value = float(sum_values / count_values) 
				# set the new single beta_value as the new methylation value for that gene
				# in correspondence of that specific aliquot, rounding it to a float with 8 decimal numbers
				dict_value[tcga_aliq] = round(single_beta_value,8)
				# add this methylation value also in the proper cell of the dataframe
				methyl_df.set_value(tcga_aliq,gene_name,round(single_beta_value,8))
				# reset the variables for the next iteration
				sum_values = 0
				count_values = 0
				single_beta_value = 0

	for i in genesSYM_of_interest: 
		methyl_df[i] = methyl_df[i].apply(lambda y: np.nan if (isinstance(y,list) and len(y)==0) else y)
	

	# Export the dataframe with the single methylation values for each gene of interest and in each TCGA aliquot
	writer = ExcelWriter('./3_TCGA_Data/Methylation/Methylation_Values.xlsx')
	methyl_df.to_excel(writer,'Sheet1')
	writer.save()
	
	return methyl_df