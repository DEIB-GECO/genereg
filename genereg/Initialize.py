# coding: utf-8

# Import libraries
import os
import urllib.request
import pickle
import pandas as pd

def library_init():

	"""
	The LIBRARY_INIT operation performs a preliminary initialization procedure and sets all the initial library settings.
	
	Example::
	
		import genereg as gr
		gr.Initialize.library_init()
	"""
	
	# Build the directories tree for storing files and results
	base_dirs = ['./0_Genes_Mapping','./1_Transcription_Factors','./2_Regulatory_Genes','./3_TCGA_Data','./4_Data_Matrix_Construction','./5_Data_Analysis','./(MaterializeResults)']
	for d in base_dirs:
		if not os.path.exists(d):
			os.makedirs(d)

	if not os.path.exists('./0_Genes_Mapping/DATA'):
		os.makedirs('./0_Genes_Mapping/DATA')

	tcga_dirs = ['./3_TCGA_Data/Methylation','./3_TCGA_Data/Gene_Expression']		
	for t in tcga_dirs:
		if not os.path.exists(t):
			os.makedirs(t)
			
	data_matrixes_dirs = ['./4_Data_Matrix_Construction/Model1','./4_Data_Matrix_Construction/Model2','./4_Data_Matrix_Construction/Model3','./4_Data_Matrix_Construction/Model4','./4_Data_Matrix_Construction/Model5']
	for m in data_matrixes_dirs:
		if not os.path.exists(m):
			os.makedirs(m)


	# Download UniProt and ENCODE lists of Transcription Factors and HGNC human genes list for creating a GENES - TFs mapping table used by some functions of this library.
	# This lists are progressively updated: it is sufficient to re-execute this script for updating them locally in your PC and make the library use the most recently updated versions.
	#
	# 'HGNC.tsv' is a file downloaded from the HGNC website, through the "Custom Downlaods" tool, containing the list of all human genes, each one of them with its current Gene Symbol,
	#  Entrez Gene ID, HGNC ID, RefSeq ID, Ensembl Gene ID and, if present, the UniProt ID related to the TF they encode.
	#  In addition, there are also the symbols that were previously used and that are now obsolete, plus other possible synonyms of these genes (8 columns: HGNC ID, Approved Symbol, Previous Symbols, Synonyms, Entrez Gene ID, RefSeq ID, UniProt ID, Ensembl Gene ID).
	#
	# 'UNIPROT.tsv' is a file downloaded from the UniProt website (http://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+reviewed%3Ayes),
	#  containing the list of all human transcription factors (with  status == 'reviewed'), with their UniProt ID and the names of the corresponding encoding genes (3 columns: UniProt_ID, TF_name and Gene_Symbol).
	#
	# 'ENCODE.tsv' is a file downloaded from the ENCODE website (https://www.encodeproject.org/report/?type=Target&investigated_as=transcription%20factor&organism.scientific_name=Homo%20sapiens&limit=all&field=label&field=gene_name),
	#  containing the list of all human transcription factors and the names of the corresponding encoding genes (2 columns: TF_name and Gene_Symbol).
	hgnc_file = urllib.request.URLopener()
	uniprot_file = urllib.request.URLopener()
	encode_file = urllib.request.URLopener()
	hgnc_file.retrieve('https://raw.githubusercontent.com/Kia23/genereg/master/DATA/HGNC.tsv', './0_Genes_Mapping/DATA/HGNC.tsv')
	uniprot_file.retrieve('https://raw.githubusercontent.com/Kia23/genereg/master/DATA/UNIPROT.tsv', './0_Genes_Mapping/DATA/UNIPROT.tsv')
	encode_file.retrieve('https://raw.githubusercontent.com/Kia23/genereg/master/DATA/ENCODE.tsv', './0_Genes_Mapping/DATA/ENCODE.tsv')


	# Load the list of genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})

	# Save the list of gene sets of interest
	gene_sets = (EntrezConversion_df.GENE_SET.unique()).tolist()
	with open ('./Gene_Sets.txt', 'w') as fp:
		for i in gene_sets:
			fp.write("%s\n" % i)
		
	# Create a dictionary setting the different gene sets of interest as keys and the list of genes of interest belonging to each of them as values
	from collections import defaultdict
	dict_gene_sets = defaultdict(list)
		
	for s in gene_sets:
		dict_gene_sets[s] = []
		
	for i, r in EntrezConversion_df.iterrows():
		n = r['GENE_SYMBOL']
		s = r['GENE_SET']
		dict_gene_sets[s].append(n)
		
	# Export the dictionary
	pickle.dump(dict_gene_sets, open('./dict_gene_sets.p', 'wb'))

				
	# Build the directories tree for storing data analysis results
	gene_sets_file = open('./Gene_Sets.txt', 'r')
	sets = gene_sets_file.read().split('\n')
	sets.remove('')
	gene_sets_file.close()

	gene_sets_dirs = []
	for s in sets:
		elem = './5_Data_Analysis/'+s
		gene_sets_dirs.append(elem)
	for d in gene_sets_dirs:
		if not os.path.exists(d):
			os.makedirs(d)

	analysis_names = ['FeatureSelection','LinearRegression']
	for s in sets:
		analysis_dirs = []
		elem_1 = './5_Data_Analysis/'+s+'/FeatureSelection'
		elem_2 = './5_Data_Analysis/'+s+'/LinearRegression'
		analysis_dirs.append(elem_1)
		analysis_dirs.append(elem_2)
		for a in analysis_dirs:
			if not os.path.exists(a):
				os.makedirs(a)

	model_names = ['M2','M3','M5']
	for s in sets:
		for n in analysis_names:
			models_dirs = []
			elem_1 = './5_Data_Analysis/'+s+'/'+n+'/M2'
			elem_2 = './5_Data_Analysis/'+s+'/'+n+'/M3'
			elem_3 = './5_Data_Analysis/'+s+'/'+n+'/M5'
			models_dirs.append(elem_1)
			models_dirs.append(elem_2)
			models_dirs.append(elem_3)
			for m in models_dirs:
				if not os.path.exists(m):
					os.makedirs(m)			

	for s in sets:
		for model in model_names:
			additional_dirs = []
			elem_1 = './5_Data_Analysis/'+s+'/LinearRegression/'+model+'/Coefficients'
			elem_2 = './5_Data_Analysis/'+s+'/LinearRegression/'+model+'/Confidence Intervals'
			elem_3 = './5_Data_Analysis/'+s+'/LinearRegression/'+model+'/Correlation Matrix'
			additional_dirs.append(elem_1)
			additional_dirs.append(elem_2)
			additional_dirs.append(elem_3)
			for ad in additional_dirs:
				if not os.path.exists(ad):
					os.makedirs(ad)