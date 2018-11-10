# coding: utf-8

# Import libraries
import gmql as gl
import pandas as pd
from pandas import ExcelWriter
import pickle
import xlsxwriter


def extract_tfs(cell_lines, gencode_version):

	"""
	The EXTRACT_TFS operation extracts, from ChIP_Seq ENCODE expriments and for assembly GRCh38, the Transcription Factors that bind to promoter regions of genes belonging to specific cell lines, filtering by 'conservative idr thresholded peaks' in order to extract higher quality region data, and removing negative audits in order to keep only high quality data samples. Intermediate results files are exported locally during the execution of the function, while the final set of trasncription factors is returned as a Python dictionary (dict_GeneTF.p), where each target gene (set as key) is associated to the list of TFs binding to its promoters (set as value). 

	:param cell_lines: a list of strings containing the names of the cell lines to analyze (it's possible to consider data from 1 up to 3 cell lines at the same time)
	:param gencode_version: number representing the GENCODE genomic annotations to use (currently, for assembly GRCh38, versions 22, 24 and 27 can be used)
	:return: a Python dictionary
	
	Example::
	
		import genereg as gr
		tfs_dict = gr.TranscriptionFactors.extract_tfs(cell_lines=['K562','MCF7'], gencode_version=22)
	"""
	
	# Check input parameters
	if not ((len(cell_lines) > 0) and (len(cell_lines) < 4)):
		raise ValueError('You have to specify from 1 up to 3 cell lines to investigate')
	
	if gencode_version not in [22, 24, 27]:
		raise ValueError('GRCh38 GENCODE versions available are 22, 24 and 27')

		
	# Execute the query for the extraction of TFs on the remote server, using the PyGMQL Python library
	gl.set_remote_address('http://gmql.eu/gmql-rest/')
	gl.login()
	gl.set_mode('remote')

	# Load the ENCODE datasets to be used in the query
	narrow_dataset = gl.load_from_remote(remote_name='GRCh38_ENCODE_NARROW_NOV_2017', owner='public')
	ann_dataset = gl.load_from_remote(remote_name='GRCh38_ANNOTATION_GENCODE', owner='public')

	# Extract NARROW data of interest
	if len(cell_lines) == 1:
		cell = cell_lines[0]
		narrow = narrow_dataset.meta_select((narrow_dataset['assay'] == 'ChIP-seq') & (narrow_dataset['output_type'] == 'conservative idr thresholded peaks') & (narrow_dataset['biosample_term_name'] == cell) & (narrow_dataset['assembly'] == 'GRCh38') & (narrow_dataset['project'] == 'ENCODE') & (narrow_dataset['file_status'] == 'released') & (~((narrow_dataset['audit_error'] == 'extremely low read depth') | (narrow_dataset['audit_error'] == 'extremely low read length') | (narrow_dataset['audit_warning'] == 'insufficient read depth') | (narrow_dataset['audit_not_compliant'] == 'insufficient read depth') | (narrow_dataset['audit_not_compliant'] == 'insufficient replicate concordance') | (narrow_dataset['audit_not_compliant'] == 'missing input control') | (narrow_dataset['audit_not_compliant'] == 'severe bottlenecking') | (narrow_dataset['audit_not_compliant'] == 'unreplicated experiment'))))
	elif len(cell_lines) == 2:
		cell_1 = cell_lines[0]
		cell_2 = cell_lines[1]
		narrow = narrow_dataset.meta_select((narrow_dataset['assay'] == 'ChIP-seq') & (narrow_dataset['output_type'] == 'conservative idr thresholded peaks') & ((narrow_dataset['biosample_term_name'] == cell_1) | (narrow_dataset['biosample_term_name'] == cell_2)) & (narrow_dataset['assembly'] == 'GRCh38') & (narrow_dataset['project'] == 'ENCODE') & (narrow_dataset['file_status'] == 'released') & (~((narrow_dataset['audit_error'] == 'extremely low read depth') | (narrow_dataset['audit_error'] == 'extremely low read length') | (narrow_dataset['audit_warning'] == 'insufficient read depth') | (narrow_dataset['audit_not_compliant'] == 'insufficient read depth') | (narrow_dataset['audit_not_compliant'] == 'insufficient replicate concordance') | (narrow_dataset['audit_not_compliant'] == 'missing input control') | (narrow_dataset['audit_not_compliant'] == 'severe bottlenecking') | (narrow_dataset['audit_not_compliant'] == 'unreplicated experiment'))))
	elif len(cell_lines) == 3:
		cell_1 = cell_lines[0]
		cell_2 = cell_lines[1]
		cell_3 = cell_lines[2]
		narrow = narrow_dataset.meta_select((narrow_dataset['assay'] == 'ChIP-seq') & (narrow_dataset['output_type'] == 'conservative idr thresholded peaks') & ((narrow_dataset['biosample_term_name'] == cell_1) | (narrow_dataset['biosample_term_name'] == cell_2) | (narrow_dataset['biosample_term_name'] == cell_3)) & (narrow_dataset['assembly'] == 'GRCh38') & (narrow_dataset['project'] == 'ENCODE') & (narrow_dataset['file_status'] == 'released') & (~((narrow_dataset['audit_error'] == 'extremely low read depth') | (narrow_dataset['audit_error'] == 'extremely low read length') | (narrow_dataset['audit_warning'] == 'insufficient read depth') | (narrow_dataset['audit_not_compliant'] == 'insufficient read depth') | (narrow_dataset['audit_not_compliant'] == 'insufficient replicate concordance') | (narrow_dataset['audit_not_compliant'] == 'missing input control') | (narrow_dataset['audit_not_compliant'] == 'severe bottlenecking') | (narrow_dataset['audit_not_compliant'] == 'unreplicated experiment'))))

	# Create the dataset of promoters
	annotations_version = str(gencode_version)
	coding_transcripts_0 = ann_dataset.meta_select((ann_dataset['release_version'] == annotations_version) & (ann_dataset['annotation_type'] == 'transcript'))
	coding_transcripts = coding_transcripts_0.reg_select((coding_transcripts_0.transcript_type == 'protein_coding') & ((coding_transcripts_0.tag == 'basic') | (coding_transcripts_0.tag == 'CCDS')))
	prom_reg = coding_transcripts.reg_project(['gene_id','gene_name','entrez_gene_id'], new_field_dict = {'start': coding_transcripts.start - 2000, 'stop': coding_transcripts.start + 1000})
	prom = prom_reg.group(regs=['gene_name'], regs_aggregates = {'ensembl_id': gl.BAGD('gene_id'), 'entrez_id': gl.BAGD('entrez_gene_id')})

	# Merge all the possible replicas of the same TF, combining them in a single sample
	full_encode = narrow.normal_cover(1, 'ANY', ['experiment_target'])

	# Extract the transcription factors that overlap with at least one promoter region
	res_0 = prom.map(full_encode, refName='prom', expName='TF')
	res_1 = res_0.reg_select(res_0.count_prom_TF > 0)

	# Encode, for each region, the corresponding TF that binds to it as a region attribute
	set_tf = res_1.reg_project(new_field_dict = {'TF': res_1['TF.experiment_target','string']})

	# Merge all the samples into a dataset with a single sample containing all the regions with their binding TFs
	# and remove regions with unknown names for their belonging gene
	merged = set_tf.merge()
	known_genes = merged.reg_select(merged.entrez_id != '')

	# Group the regions by name, setting in the region attribute 'TF' the list of transcription factors that bind to that gene's promoters
	res = known_genes.group(regs=['gene_name'], regs_aggregates = {'ensembl_gene_id': gl.BAGD('ensembl_id'), 'entrez_gene_id': gl.BAGD('entrez_id'), 'TFs': gl.BAGD('TF')})

	# Materialize the results into a GDataframe
	res_Gdf = res.materialize('./(MaterializeResults)')


	# Extract the regions dataframe, where each row corresponds to a region and each column to an attribute
	GeneTF_df = res_Gdf.regs
	# Check the length of the dataframe, that is the number of rows of the dataframe
	length_df = len(GeneTF_df)
	# Set progressive integer numbers as new indexes of the dataframe
	GeneTF_df.index = range(length_df)

	# Convert all columns names into uppercase letters and rename them
	GeneTF_df.columns = map(str.upper, GeneTF_df.columns)
	GeneTF_df.rename(columns={'CHR':'CHROM','START':'LEFT','STOP':'RIGHT','GENE_NAME':'GENE_SYMBOL','TFS':'TFs'}, inplace=True)

	for index, row in GeneTF_df.iterrows():
		tfs_str = row['TFs']
		tfs_list = tfs_str.split(',')
		GeneTF_df.set_value(index,'TFs',tfs_list)


	# Load the list of genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})
	
	# Create a list with the Gene Symbols of the genes of interest
	Symbols = []
	for index, row in EntrezConversion_df.iterrows():
		i = row['GENE_SYMBOL']
		Symbols.append(i)
	N_genes = len(Symbols)


	# Create an empty dictionary with lists as values for each key
	from collections import defaultdict
	dict_GeneTF = defaultdict(list)

	# Set the keys and initialize their values as empty lists
	for v in Symbols:
		dict_GeneTF[v] = []
	dict_length = len(dict_GeneTF)

	# Select from the GeneTF_df only the rows with Gene Symbols of target genes of interest
	for index, row in GeneTF_df.iterrows():   # iterate along the whole dataframe
		# get the current row GENE_SYMBOL
		i = row['GENE_SYMBOL']                  
		for value in Symbols:  # check if the current gene is contained in the list of genes of interest
			if i == value:  # if there's correspondence
				TrFa_list = row.TFs # extract the list of TFs
				for t in TrFa_list:
					# since each gene can have more than one promoter bound by the same TF,
					# only distinct values for each transcription factor should be inserted in the dictionary 
					if t not in dict_GeneTF[i]:
						# add the transcription factor to the list of values corresponding to that gene
						dict_GeneTF[i].append(t)

	# Order alphabetically the list of TFs for each gene of interest
	for k in dict_GeneTF.keys():
		old = dict_GeneTF[k] # get the list of TFs
		sorted_TFs = sorted(old, key=lambda s: s.lower()) # sort the list alphabetically (case-insensitive sorting)
		dict_GeneTF[k] = sorted_TFs


	# For each gene, add to the dictionary its ENTREZ_GENE_ID and the GENE_SET it belongs to.
	# These are added as two lists at the end of the value list of each key (i.e. gene): clearly the list containing
	# the ENTREZ_GENE_ID will always have one string element for each key, while the list containing the gene sets
	# can have one or more elements, depending on the number of sets the corresponding gene belongs to.

	# Get distinct Gene Symbols of genes of interest (considering only once the genes that belongs to multiple sets)
	Symbols_distinct = []
	for value in Symbols:
		if value not in Symbols_distinct:
			Symbols_distinct.append(value)

	for value in Symbols_distinct:
		row = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == value]
		# get the ENTREZ_GENE_ID (in case of gene belonging to multiple sets this list will contain
		# its ENTREZ_GENE_ID as many times as the number of sets it belongs to)
		row_entrez_id = list(row.ENTREZ_GENE_ID) 
		N_eid = len(row_entrez_id)
		if N_eid > 1:
			entrez_id = list(list(row_entrez_id[:N_eid-1])) # consider the ENTREZ_GENE_ID only once
		else: entrez_id = list(row_entrez_id)    
		sets = list(list(row.GENE_SET))
		# add the ENTREZ_GENE_ID and the GENE_SET as elements in the list of values corresponding to the proper gene
		dict_GeneTF[value].append(entrez_id)
		dict_GeneTF[value].append(sets)   

	# So, the general form of this dictionary containing the information we need is the following:
	# dict_GeneTF = {key: value, ...} = {GENE_SYMBOL: [TF1, TF2, TF3, ..., [ENTREZ_GENE_ID], [GENE_SETs]]}

	
	# Store the number of TFs for each gene of interest in a new dictionary
	from collections import defaultdict
	dict_TFs_gene = defaultdict(int)

	# Initialize the dictionary setting the keys and their initial values
	for k in dict_GeneTF.keys():
		dict_TFs_gene[k] = 0
		
	# Set the number of TFs that bind to a gene's promotes as value of the corresping key in the dictionary
	for k in dict_GeneTF.keys():
		transcription_factors = dict_GeneTF[k][:-2]
		number_TFs = len(transcription_factors)
		dict_TFs_gene[k] = number_TFs
	   
	# Generate an histogram showing the previous distribution, that is the number of TFs that bind
	# to each gene's promoters (sorting this number from the highest to the smallest)

	# Convert the dictionary into a pandas dataframe
	TFs_gene_unsorted_df = pd.DataFrame(list(dict_TFs_gene.items()), columns=['GENE_SYMBOL', '#TFs'])

	# Sort the dataframe according to the number of TFs for each gene
	TFs_gene_df = TFs_gene_unsorted_df.sort_values(by='#TFs', ascending=0)

	# Add to the dataframe a column for storing also the Entrez Gene ID of each gene besides the already present Gene Symbols
	TFs_gene_df['ENTREZ_GENE_ID'] = ''

	# Add the correct Entrez Gene ID for each gene
	for index, row in TFs_gene_df.iterrows():
		sym = row['GENE_SYMBOL']
		eid = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == sym, 'ENTREZ_GENE_ID'].iloc[0]
		TFs_gene_df.set_value(index,'ENTREZ_GENE_ID',eid)

	# Export the dataframe into an Excel file
	writer = ExcelWriter('./1_Transcription_Factors/Number of TFs for each gene of interest.xlsx')
	TFs_gene_df.to_excel(writer,'Sheet1',index=False)
	writer.save()

	# Export the dictionary of genes of interest and their TFs, ENTREZ GENE ID and GENE_SETs
	# Save the dictionary into a pickle file
	pickle.dump(dict_GeneTF, open('./1_Transcription_Factors/dict_GeneTF.p', 'wb'))
	 
	# Save the dictionary as a .xlsx file
	workbook = xlsxwriter.Workbook('./1_Transcription_Factors/dict_GeneTF.xlsx')
	worksheet = workbook.add_worksheet()
	# Set the headers of the columns
	worksheet.write(0,0,'GENE_SYMBOL')
	worksheet.write(0,1,'Transcription Factors')

	row = 1
	col = 0
	# Print the dictionary
	for key in dict_GeneTF.keys():
		row += 1
		worksheet.write(row, col, key)
		for item in dict_GeneTF[key]:
			worksheet.write(row, col + 1, ''.join(item))
			if item == dict_GeneTF[key][-2]: # the second to last element is the Entrez Gene ID
				worksheet.write(row, col + 2, 'Entrez Gene ID')
			if item == dict_GeneTF[key][-1]: # the last element is the list of gene sets
				worksheet.write(row, col + 2, 'Gene Set')
			row += 1
	workbook.close()

	# Save the dictionary as a .txt file
	with open ('./1_Transcription_Factors/dict_GeneTF.txt', 'w') as fp:
		for p in dict_GeneTF.items():
			fp.write('%s : %s\n\n' % p)
			
	return dict_GeneTF