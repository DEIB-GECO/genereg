# coding: utf-8

# Import libraries
import pandas as pd
from pandas import ExcelWriter
import pickle
import xlsxwriter


def extract_regulatory_genes():
    
	"""
	The EXTRACT_REGULATORY_GENES operation extracts from the set of Transcription Factors associated to a gene, the list of its candidate regulatory genes, i.e., the genes that encode for those TFs. Intermediate results files are exported locally during the execution of the function, while the final set of trasncription factors is returned as a Python dictionary (dict_RegulGenes.p), where each target gene (set as key) is associated to the list of its candidate regulatory genes (set as value).

	:return: a Python dictionary
	
	Example::
	
		import genereg as gr
		reg_genes_dict = gr.RegulatoryGenes.extract_regulatory_genes()
	"""


	# Starting from the dictionary containing for each gene of interest the TFs that bind to its promoters,
	# extract the names of the genes encoding the TFs in order to identify the candidate regulatory genes of each gene of interest
	dict_GeneTF = pickle.load(open('./1_Transcription_Factors/dict_GeneTF.p', 'rb'))
	TFs_interest = []
	for key, value in dict_GeneTF.items():
		TFs = value[:-2]  # the TFs are all the elements of the value list, except for the last two
		for tf in TFs:
			if tf not in TFs_interest:
				TFs_interest.append(tf)

	# Import the gene-TFs mapping dataframe 
	Mapping_df = pd.read_excel('./0_Genes_Mapping/Genes_Mapping.xlsx',sheetname='Sheet1',header=0,converters={'ENTREZ_GENE_ID':str,'HGNC_ID':str})

	for index, row in Mapping_df.iterrows():
		tfs_str = row['TF']
		if isinstance(tfs_str,str):
			tfs_list = tfs_str.split(', ')
		else:
			tfs_list = []
		Mapping_df.set_value(index,'TF',tfs_list)

	# Extract in a list all the names of the TFs contained in the mapping dataframe
	mapping_df_TFs = []
	for index, row in Mapping_df.iterrows():
		tfs = row['TF']
		if len(tfs) != 0:
			for t in tfs:
				if t not in mapping_df_TFs:
					mapping_df_TFs.append(t)

	# Create a reduced dataframe with all the distinct TFs and their encoding genes, filtering only the TFs of interest previously extracted 
	distinct_TFs = []
	for index, row in Mapping_df.iterrows():
		tfs = row['TF']
		if len(tfs) != 0:
			for t in tfs:
				if t in TFs_interest:
					if t not in distinct_TFs:
						distinct_TFs.append(t)

	from collections import defaultdict
	dict_tf_gene = defaultdict(list)

	for t in distinct_TFs:
		dict_tf_gene[t] = []

	for index, row in Mapping_df.iterrows():
		tf = row['TF']
		gene = row['GENE_SYMBOL']
		for t in tf:
			if t in distinct_TFs:
				dict_tf_gene[t].append(gene)

	TF_Gene_df = pd.DataFrame(list(dict_tf_gene.items()), columns=['TF_NAME', 'GENE_SYMBOL'])
	for index, row in TF_Gene_df.iterrows():
		genes = row['GENE_SYMBOL']
		if len(genes) == 1:
			new_gene = ''.join(genes)
			TF_Gene_df.set_value(index,'GENE_SYMBOL',new_gene)

		
	# Create a new empty dictionary with lists as values for each key (gene)
	from collections import defaultdict
	dict_RegulGenes = defaultdict(list)

	# Set the keys and initialize their values as empty lists
	for v in dict_GeneTF.keys():
		dict_RegulGenes[v] = []

	# Get the TFs of each target gene and extract the names of the genes encoding them from the mapping dataframe
	for key, value in dict_GeneTF.items():
		TFs = value[:-2]  # the TFs are all the elements of the value list, except for the last two
		for tf in TFs:
			# for each TF, search in the mapping dataframe for the name of the encoding gene
			if tf in mapping_df_TFs:
				# get the name (GENE_SYMBOL) of the gene encoding the transcription factor "tf"
				gene_name = TF_Gene_df.loc[TF_Gene_df['TF_NAME'] == tf, 'GENE_SYMBOL'].iloc[0]
				# add the regulatory gene in correspondence of the proper gene in the dictionary
				dict_RegulGenes[key].append(gene_name)
			# in case the transcription factor considered is not mapped in the dataframe,
			# then the name of its encoding gene is unknown ('n/a')
			else: dict_RegulGenes[key].append('n/a')


	# SUMMARY TABLE summarizing for each gene of interest the TFs binding to its promoters and their corresponding encoding genes:

	# Each row of the table is indexed by the Gene Symbols of the genes of interest and progressive integers representing the number of TFs for each gene
	genes_of_interest = []
	for k in dict_GeneTF.keys():   
		genes_of_interest.append(k)

	# Extract the highest number of regulatory genes for a single gene of interest
	highest_n = 0
	for k, v in dict_RegulGenes.items():
		n = len(value)
		if n > highest_n:
			highest_n = n
	top_range = highest_n + 100

	# Define the number of rows in the table for each gene of interest		
	num_lst = []
	for i in list(range(1,top_range)):
		num_lst.append(i)

	# Cartesian product to generate tuples for multi-indexing
	import itertools
	tuples = []
	for i in itertools.product(genes_of_interest,num_lst):
		tuples.append(i)

	# Set the multiple indexes to be used in the dataframe
	index = pd.MultiIndex.from_tuples(tuples, names=['GENE_SYMBOL', '#'])

	# Create the dataframe and initialize the empty cells as empty strings
	info_genes_of_interest = pd.DataFrame('', index = index, columns = ['Transcription Factors','Regulatory Genes','Entrez_Gene_IDs','ENTREZ_GENE_ID','GENE_SET','#TFs','#RegulatoryGenes (distinct)']) 

	# Set the correct Entrez Gene ID for each gene of interest
	for index, row in info_genes_of_interest.iterrows():
		sym = index[0]
		n = index[1]
		if n == 1:
			eid = Mapping_df.loc[Mapping_df['GENE_SYMBOL'] == sym, 'ENTREZ_GENE_ID'].iloc[0]
			info_genes_of_interest.loc[(sym, n),'ENTREZ_GENE_ID'] = eid
			
	# Set the gene sets
	for key, value in dict_GeneTF.items():
		# get the list of gene sets associated to gene 'key'
		# (i.e. the last element of the list related to gene 'key')
		sets = value[-1]
		# set the list of gene sets to the correct cell in the dataframe (in correspondence of index 'key')
		n_path = len(sets)
		if n_path == 1:
			info_genes_of_interest.loc[(key, 1),'GENE_SET'] = sets[0]
		if n_path == 2:
			info_genes_of_interest.loc[(key, 1),'GENE_SET'] = sets[0]
			info_genes_of_interest.loc[(key, 2),'GENE_SET'] = sets[1]
		if n_path == 3:
			info_genes_of_interest.loc[(key, 1),'GENE_SET'] = sets[0]
			info_genes_of_interest.loc[(key, 2),'GENE_SET'] = sets[1]
			info_genes_of_interest.loc[(key, 3),'GENE_SET'] = sets[2]
			
	# Set the TFs
	for key, value in dict_GeneTF.items():
		# get the TFs (i.e. the list of values except for the last two elements)
		tfs = value[:-2]
		# set the list of TFs to the correct cell in the dataframe (in correspondence of index 'key')
		for i in num_lst:
			if i <= len(tfs):
				info_genes_of_interest.loc[(key, i),'Transcription Factors'] = tfs[i-1]
				
	# Set the regulatory genes (both with their Gene Symbols and Entrez Gene IDs)
	for key, value in dict_RegulGenes.items():
		# the set of regulatory genes is the list 'value' corresponding to each key (gene).
		# Set the list of regulatory genes to the correct cell in the dataframe (in correspondence of index 'key')
		for i in num_lst:
			if i <= len(value):
				info_genes_of_interest.loc[(key, i),'Regulatory Genes'] = value[i-1]
				if value[i-1] == 'n/a':
					eid = 'n/a'
				else:
					# get the Entrez Gene ID of the regulatory gene
					eid = Mapping_df.loc[Mapping_df['GENE_SYMBOL'] == value[i-1], 'ENTREZ_GENE_ID'].iloc[0]
				info_genes_of_interest.loc[(key, i),'Entrez_Gene_IDs'] = eid
		
		
	# Remove the empty rows in the dataframe
	for index, row in info_genes_of_interest.iterrows():
		tfs = row['Transcription Factors']
		path = row['GENE_SET']
		if (tfs == '') & (path == ''):
			info_genes_of_interest.drop(index, inplace=True)


	# Extract the distinct candidate regulatory genes for each gene of interest:

	# Remove from the dictionary containing regulatory genes the duplicated genes, if present, in order to have a dictionary with all the distinct candidate regulatory genes for each gene of interest
	for v in dict_GeneTF.keys():
		dict_RegulGenes[v] = []

	for key, value in dict_GeneTF.items():
		TFs = value[:-2]
		for tf in TFs:   
			if tf in mapping_df_TFs:
				gene_name = TF_Gene_df.loc[TF_Gene_df['TF_NAME'] == tf, 'GENE_SYMBOL'].iloc[0]
				if gene_name not in dict_RegulGenes[key]:
					dict_RegulGenes[key].append(gene_name)

	# So, the general form of this second dictionary containing the information about regulatory genes is the following:
	# dict_RegulGenes = {key: value, ...} = {GENE_SYMBOL: [REG_GENE1, REG_GENE2, REG_GENE3, ...]}, where each regulatory gene is identified by its GENE_SYMBOL


	# Export the dictionary of genes of interest and their regulatory genes:

	# Save the dictionary into a pickle file
	pickle.dump(dict_RegulGenes, open('./2_Regulatory_Genes/dict_RegulGenes.p', 'wb'))

	# Only for the sake of clearness, order alphabetically the list of candidate regulatory genes for each gene of interest
	dict_RegulGenes_ord = dict_RegulGenes.copy()
	for k in dict_RegulGenes_ord.keys():
		old = dict_RegulGenes_ord[k]
		sorted_genes = sorted(old)
		dict_RegulGenes_ord[k] = sorted_genes

	# Save the dictionary as a .xlsx file
	workbook = xlsxwriter.Workbook('./2_Regulatory_Genes/dict_RegulGenes.xlsx')
	worksheet = workbook.add_worksheet()
	# Set the headers of the columns
	worksheet.write(0,0,'GENE_SYMBOL')
	worksheet.write(0,1,'ENTREZ_GENE_ID')
	worksheet.write(0,2,'Distinct Regulatory Genes - GENE_SYMBOL')
	worksheet.write(0,3,'Distinct Regulatory Genes - ENTREZ_GENE_ID')

	row = 1
	col = 0
	for key in dict_RegulGenes_ord.keys():
		row += 1
		worksheet.write(row, col, key)
		# get the ENtrez Gene ID of the gene of interest
		eid = Mapping_df.loc[Mapping_df['GENE_SYMBOL'] == key, 'ENTREZ_GENE_ID'].iloc[0]
		worksheet.write(row, col + 1, ''.join(eid))
		
		for item in dict_RegulGenes_ord[key]:
			worksheet.write(row, col + 2, ''.join(item))
			# get the Entrez Gene ID of the regulatory gene
			if item == 'PTRF':
				entrez_id = Mapping_df.loc[Mapping_df['GENE_SYMBOL'] == 'CAVIN1', 'ENTREZ_GENE_ID'].iloc[0]
			else:
				entrez_id = Mapping_df.loc[Mapping_df['GENE_SYMBOL'] == item, 'ENTREZ_GENE_ID'].iloc[0]
			worksheet.write(row, col + 3, ''.join(entrez_id))
			row += 1

	workbook.close()

	# Save the dictionary as a .txt file
	with open ('./2_Regulatory_Genes/dict_RegulGenes.txt', 'w') as fp:
		for p in dict_RegulGenes_ord.items():
			fp.write('%s : %s\n\n' % p)

			
	# Count the number of TFs and distinct regulatory genes for each gene of interest:

	# Store the number of TFs and distinct regulatory genes for each gene of interest in two dictionaries
	from collections import defaultdict
	dict_TFs_genes = defaultdict(int)
	dict_regul_genes = defaultdict(int)

	for k in dict_GeneTF.keys():
		dict_TFs_genes[k] = 0
		
	for k in dict_RegulGenes.keys():
		dict_regul_genes[k] = 0
		
	for k in dict_GeneTF.keys():
		transcription_factors = dict_GeneTF[k][:-2]
		number_TFs = len(transcription_factors)
		dict_TFs_genes[k] = number_TFs

	for k in dict_RegulGenes.keys():
		genes = dict_RegulGenes[k]
		number_genes = len(genes)
		dict_regul_genes[k] = number_genes
    
    
	# Create a table summarizing for each gene of interest the number of TFs binding to its promoters and the number of distinct genes encoding them
	TFs_genes_df = pd.DataFrame(list(dict_TFs_genes.items()), columns=['GENE_SYMBOL', '#TFs'])
	TFs_genes_df.set_index('GENE_SYMBOL', inplace=True)
	regul_genes_df = pd.DataFrame(list(dict_regul_genes.items()), columns=['GENE_SYMBOL', '#RegulatoryGenes (distinct)'])
	regul_genes_df.set_index('GENE_SYMBOL', inplace=True)

	# Join the two dataframes into a single one to have both information together
	TFs_regul_genes_df = TFs_genes_df.join(regul_genes_df)
	TFs_regul_genes_df['GENE_SYMBOL'] = TFs_regul_genes_df.index
	TFs_regul_genes_df.index = range(len(TFs_regul_genes_df))  # set a new progressive index for this table

	# Add to the dataframe a column for storing also the Entrez Gene ID of each gene, besides the already present Gene Symbol
	TFs_regul_genes_df['ENTREZ_GENE_ID'] = ''

	# Add the correct Entrez Gene ID for each gene
	for index, row in TFs_regul_genes_df.iterrows():
		sym = row['GENE_SYMBOL']
		eid = Mapping_df.loc[Mapping_df['GENE_SYMBOL'] == sym, 'ENTREZ_GENE_ID'].iloc[0]
		TFs_regul_genes_df.set_value(index,'ENTREZ_GENE_ID',eid)

	TFs_regul_genes_df_final = TFs_regul_genes_df[['GENE_SYMBOL','ENTREZ_GENE_ID','#TFs','#RegulatoryGenes (distinct)']].copy()
		
	for index, row in TFs_regul_genes_df_final.iterrows():
		gene = row['GENE_SYMBOL']
		n_tfs = row['#TFs']
		n_genes_reg = row['#RegulatoryGenes (distinct)']
		info_genes_of_interest.loc[(gene, 1),'#TFs'] = n_tfs
		info_genes_of_interest.loc[(gene, 1),'#RegulatoryGenes (distinct)'] = n_genes_reg

	# Export the dataframe as a .xlsx file
	writer = ExcelWriter('./2_Regulatory_Genes/Full_TFs-RegulatoryGenes_SUMMARY_Table.xlsx')
	info_genes_of_interest.to_excel(writer,'Sheet1')
	writer.save()
	
	return dict_RegulGenes