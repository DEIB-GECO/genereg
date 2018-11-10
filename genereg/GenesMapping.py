# coding: utf-8

# Import libraries
import pandas as pd
from pandas import ExcelWriter


def genes_mapping():

	"""
	The GENES_MAPPING operation creates creates a mapping table for all the human genes (downloaded from HGNC), providing their current official Gene Symbols, their main IDs, and, if existing, the transcription factors they encode (with their corresponding UniProt IDs). The mapping table is returned as a Pandas dataframe and exported locally in the Excel file 'Genes Mapping.xlsx'.

	:return: a Pandas dataframe

	Example::

		import genereg as gr
		mapping_df = gr.GenesMapping.genes_mapping()
	"""
	
    # Load input data files:
	
	# HGNC
	hgnc_df = pd.read_csv('./0_Genes_Mapping/DATA/HGNC.tsv', sep='\t', names=['HGNC_ID','GENE_SYMBOL','Previous Symbols','Synonyms','ENTREZ_GENE_ID','RefSeq_ID','UniProt_ID','ENSEMBL_GENE_ID'], dtype={'HGNC_ID':str,'ENTREZ_GENE_ID':str},skiprows=1)

	# Convert the elements in the dataframe that has multiple values into a list of values
	for index, row in hgnc_df.iterrows():
		old = row['Previous Symbols']
		syn = row['Synonyms']
		uni = row['UniProt_ID']
		if isinstance(old,str):
			new_old = old.split(', ')
			hgnc_df.set_value(index,'Previous Symbols',new_old)
		if isinstance(syn,str):
			new_syn = syn.split(', ')
			hgnc_df.set_value(index,'Synonyms',new_syn)
		if isinstance(uni,str):
			new_uni = uni.split(', ')
			hgnc_df.set_value(index,'UniProt_ID',new_uni)

	# Update the HGNC_ID by removing the prefix 'HGNC:'
	for index, row in hgnc_df.iterrows():
		old_hgnc = row['HGNC_ID']
		new_hgnc = old_hgnc.split(':')
		hgnc_df.set_value(index,'HGNC_ID',new_hgnc[1])
    
	# Add a new column for storing the TFs
	tf_col = pd.DataFrame(index = hgnc_df.index, columns = ['TF'])
	hgnc_df = hgnc_df.join(tf_col)


	# UNIPROT
	uni_df = pd.read_csv('./0_Genes_Mapping/DATA/UNIPROT.tsv', sep='\t', names=['UniProt_ID','TF_NAME','GENE_SYMBOL'],skiprows=1)

	# Replace the suffix '_HUMAN' with '-human'
	for index, row in uni_df.iterrows():
		name = row['TF_NAME']
		new_name = name.replace('_HUMAN', '-human')
		uni_df.loc[uni_df['TF_NAME'] == name, 'TF_NAME'] = new_name


	# ENCODE
	enc_full_df = pd.read_csv('./0_Genes_Mapping/DATA/ENCODE.tsv', sep='\t', names=['TF_NAME','GENE_SYMBOL','UniProt_ID'],skiprows=2)

	# Add the suffix '-human' to the name of all TFs
	for index, row in enc_full_df.iterrows():
		name = row['TF_NAME']
		new_name = name+'-human'
		enc_full_df.loc[enc_full_df['TF_NAME'] == name, 'TF_NAME'] = new_name

	# Reorder the columns as in UniProt
	enc_df = enc_full_df[['UniProt_ID','TF_NAME','GENE_SYMBOL']].copy()


	# Generate a unique dataframe of TFs:
	
	# Create a list with all the TFs contained in the UniProt dataframe
	transc_factors_UNI = []
	for index, row in uni_df.iterrows():
		tf = row['TF_NAME']
		transc_factors_UNI.append(tf)
		
	# Create a list with all the TFs contained in the ENCODE dataframe
	transc_factors_ENC = []
	for index, row in enc_df.iterrows():
		tf = row['TF_NAME']
		transc_factors_ENC.append(tf)
    
	# Remove from the ENCODE dataframe all the rows that contain a TF that is already present in UniProt
	for index, row in enc_df.iterrows():
		tf = row['TF_NAME']
		if tf in transc_factors_UNI:
			enc_df.drop(index, inplace=True)
        
	# Put the two dataframe together, in order to have a single dataframe with all the transcription factors and the names of the corresponding encoding genes
	tf_df = uni_df.append(enc_df, ignore_index=True)


	# Genes actualization step:
	# for all the genes in the TFs dataframe, check if the Gene Symbols are the current and official ones, otherwise replace the obsolete names or synonyms

	# Convert the GENE_SYMBOL attribute into a list in order to handle those TFs that are codified by multiple genes
	tf_df_list = tf_df.copy()
	for index, row in tf_df_list.iterrows():
		gene = row['GENE_SYMBOL']
		if isinstance(gene,str):
			gene_list = gene.split('; ')
			tf_df_list.set_value(index,'GENE_SYMBOL',gene_list)

	# Create a list with all the Gene Symbols in the 'hgnc_df'
	gene_symbols_HGNC = []
	for index, row in hgnc_df.iterrows():
		name = row['GENE_SYMBOL']
		gene_symbols_HGNC.append(name)

	# Create a list with all the 'Previous Symbols' in the HUGO dataframe
	prev_symbols = []
	for index, row in hgnc_df.iterrows():
		prev = row['Previous Symbols']
		if not(isinstance(prev,float)):
			for i in prev:
				prev_symbols.append(i)
				
	# Create a list with all the 'Synonyms' in the HUGO dataframe
	synonyms = []
	for index, row in hgnc_df.iterrows():
		syn = row['Synonyms']
		if not(isinstance(prev,float)):
			for i in syn:
				synonyms.append(i)

	# Replace obsolete names and gene synonyms
	codifying_genes_to_update = []
	already_updated = []
	for index, row in tf_df_list.iterrows():
		genes = row['GENE_SYMBOL']
		if not(isinstance(genes,float)):
			for g in genes:
				if g in gene_symbols_HGNC:
					already_updated.append(g)
				else:
					if g not in codifying_genes_to_update:
						codifying_genes_to_update.append(g)

	# Create a dictionary to store the current and official name for each of the genes to update
	from collections import defaultdict
	dict_update_symbols = defaultdict(list)

	# The keys of this dictionary are the Gene Symbols to update
	for g in codifying_genes_to_update:
		dict_update_symbols[g] = []

	# Consider obsolete Gene Symbols first
	for index, row in hgnc_df.iterrows():
		sym = row['GENE_SYMBOL']
		prev = row['Previous Symbols']
		for gene in codifying_genes_to_update:
			if not(isinstance(prev,float)) and gene in prev:
				dict_update_symbols[gene].append(sym)
				
	# Consider synonyms by taking into account also UniProt_IDs
	for index, row in hgnc_df.iterrows():
		sym = row['GENE_SYMBOL']
		syn = row['Synonyms']
		uni = row['UniProt_ID']
		for gene in codifying_genes_to_update:
			if len(dict_update_symbols[gene]) == 0:
				if not(isinstance(syn,float)) and gene in syn:
					real_uni = tf_df.loc[tf_df['GENE_SYMBOL'] == gene, 'UniProt_ID'].iloc[0]
					if real_uni in uni:
						dict_update_symbols[gene].append(sym)

	# Check the genes that remain not updated
	count_no_updated_genes = 0 
	for key, value in dict_update_symbols.items():
		if len(value) == 0:
			count_no_updated_genes = count_no_updated_genes + 1 

	updated_genes = []
	for key, value in dict_update_symbols.items():
		if len(value) != 0:
			updated_genes.append(key)


	# Update the gene names according to the dictionary created in the TFs dataframe
	for index, row in tf_df_list.iterrows():
		genes = row['GENE_SYMBOL']
		if not(isinstance(genes,float)):
			new_genes = []
			for g in genes:
				if g in updated_genes:
					new_sym = dict_update_symbols[g][0]
					new_genes.append(new_sym)
				else:
					new_genes.append(g)
			tf_df_list.set_value(index,'GENE_SYMBOL',new_genes)

	# Create a dictionary to store all the encoding genes and the list of their TFs 
	protein_coding_genes = []
	for index, row in tf_df_list.iterrows():
		genes = row['GENE_SYMBOL']
		if not(isinstance(genes,float)):
			for g in genes:
				if g not in protein_coding_genes:
					protein_coding_genes.append(g)


	from collections import defaultdict
	dict_coding_gene_tfs = defaultdict(list)

	for g in protein_coding_genes:
		dict_coding_gene_tfs[g] = []

	for index, row in tf_df_list.iterrows():
		tf = row['TF_NAME']
		genes = row['GENE_SYMBOL']
		if not(isinstance(genes,float)):
			for g in genes:
				if g in protein_coding_genes:
					dict_coding_gene_tfs[g].append(tf)

					
	# Add to the final dataframe the TFs for each gene codifying them
	in_HGNC = list(set(protein_coding_genes) & set(gene_symbols_HGNC))
	not_in_HGNC = list(set(protein_coding_genes) - set(in_HGNC))

	tfs_added = []
	for index, row in hgnc_df.iterrows():
		gene = row['GENE_SYMBOL'] 
		if gene in in_HGNC:
			tfs_list = dict_coding_gene_tfs[gene]
			hgnc_df.set_value(index,'TF',tfs_list)
			tfs_added = tfs_added + tfs_list

	add_df = pd.DataFrame(index=range(len(not_in_HGNC)), columns=['HGNC_ID','GENE_SYMBOL','Previous Symbols','Synonyms','ENTREZ_GENE_ID','RefSeq_ID','UniProt_ID','ENSEMBL_GENE_ID','TF'])
	add_df['GENE_SYMBOL'].iloc[0:len(not_in_HGNC)] = not_in_HGNC
	for index, row in add_df.iterrows():
		gene = row['GENE_SYMBOL']
		tfs_list = dict_coding_gene_tfs[gene]
		add_df.set_value(index,'TF',tfs_list)
		tfs_added = tfs_added + tfs_list
		uni_list = []
		for t in tfs_list:
			uni = tf_df.loc[tf_df['TF_NAME'] == t, 'UniProt_ID'].iloc[0]
			if not(isinstance(uni,float)):
				uni_list.append(uni)
		add_df.set_value(index,'UniProt_ID',uni_list)

	# Put the two dataframe together to create the final mapping table for all genes
	final_df = hgnc_df.append(add_df, ignore_index=True)

	for index, row in final_df.iterrows():
		uni_list = row['UniProt_ID']
		tfs_list = row['TF']
		if not(isinstance(uni_list,float)):
			uni_str = ', '.join(uni_list)
			final_df.set_value(index,'UniProt_ID',uni_str)
		if not(isinstance(tfs_list,float)):
			tfs_str = ', '.join(tfs_list)
			final_df.set_value(index,'TF',tfs_str)    
			  

	# Convert list elements into strings in the mapping dataframe
	for index, row in final_df.iterrows():
		prev = row['Previous Symbols']
		syn = row['Synonyms']
		
		if not(isinstance(prev,float)):
			new_prev = ', '.join(prev)
			final_df.set_value(index,'Previous Symbols',new_prev)
		if not(isinstance(syn,float)):
			new_syn = ', '.join(syn)
			final_df.set_value(index,'Synonyms',new_syn)
			
			
	# Sort the dataframe by gene name
	mapping_table = final_df.sort_values('GENE_SYMBOL').copy()

	# Extract and reorder only columns useful for other users
	mapping_table_filtered = mapping_table[['GENE_SYMBOL','ENTREZ_GENE_ID','ENSEMBL_GENE_ID','HGNC_ID','RefSeq_ID','TF','UniProt_ID']].copy()

	# Export the mapping dataframe as a .xlsx file
	writer = ExcelWriter('./0_Genes_Mapping/Genes_Mapping.xlsx')
	mapping_table_filtered.to_excel(writer,'Sheet1',index=False)
	writer.save()
	
	return mapping_table_filtered