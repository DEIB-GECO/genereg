# coding: utf-8

# Import libraries
import os
import urllib.request
import pickle

def library_init():

	"""
	The LIBRARY_INIT operation performs a preliminary initialization procedure and sets all the initial library settings.
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
	hgnc_file.retrieve('https://raw.githubusercontent.com/Kia23/genexpreg_cancer/master/DATA/HGNC.tsv', './0_Genes_Mapping/DATA/HGNC.tsv')
	uniprot_file.retrieve('https://raw.githubusercontent.com/Kia23/genexpreg_cancer/master/DATA/UNIPROT.tsv', './0_Genes_Mapping/DATA/UNIPROT.tsv')
	encode_file.retrieve('https://raw.githubusercontent.com/Kia23/genexpreg_cancer/master/DATA/ENCODE.tsv', './0_Genes_Mapping/DATA/ENCODE.tsv')


	# Load the list of genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'PATHWAY':str,'SubClass':str})
	
	# Convert the 'SubClass' attribute into a list
	for index, row in EntrezConversion_df.iterrows():
		subclasses = row['SubClass']
		if isinstance(subclasses,str):
			subclasses_list = subclasses.split(', ')
		else:
			s = ''	
			subclasses_list = s.split(', ')
			subclasses_list.remove('')
		EntrezConversion_df.set_value(index,'SubClass',subclasses_list)

	# Save the list of genomic pathways of interest
	genomic_pathways = (EntrezConversion_df.PATHWAY.unique()).tolist()
	with open ('./Pathways_of_Interest.txt', 'w') as fp:
		for i in genomic_pathways:
			fp.write("%s\n" % i)
		
	# Create a dictionary setting the different genomic pathways of interest as keys and the list of genes of interest belonging to eahc of them as values
	from collections import defaultdict
	dict_pathways = defaultdict(list)
		
	for p in genomic_pathways:
		dict_pathways[p] = []
		
	for i, r in EntrezConversion_df.iterrows():
		n = r['GENE_SYMBOL']
		p = r['PATHWAY']
		dict_pathways[p].append(n)		
		
	# Export the dictionary
	pickle.dump(dict_pathways, open('./dict_pathways.p', 'wb'))

	# Save the list of functional subclasses for each pathway of interest
	for p in genomic_pathways:
		current_pathway_subclasses = []
		for i, r in EntrezConversion_df.iterrows():
			pathway = r['PATHWAY']
			subclasses = r['SubClass']
			if pathway == p:
				for s in subclasses:
					if s not in current_pathway_subclasses:
						current_pathway_subclasses.append(s)
		with open ('./'+p+'_SubClasses.txt', 'w') as fp:
			for i in current_pathway_subclasses:
				fp.write('%s\n' % i)
				
				
	# Build the directories tree for storing data analysis results
	pathway_file = open('./Pathways_of_Interest.txt', 'r')
	pathways = pathway_file.read().split('\n')
	pathways.remove('')
	pathway_file.close()

	pathway_dirs = []
	for p in pathways:
		elem = './5_Data_Analysis/'+p
		pathway_dirs.append(elem)
	for d in pathway_dirs:
		if not os.path.exists(d):
			os.makedirs(d)

	analysis_names = ['FeatureSelection','LinearRegression']
	for p in pathways:
		analysis_dirs = []
		elem_1 = './5_Data_Analysis/'+p+'/FeatureSelection'
		elem_2 = './5_Data_Analysis/'+p+'/LinearRegression'
		analysis_dirs.append(elem_1)
		analysis_dirs.append(elem_2)
		for a in analysis_dirs:
			if not os.path.exists(a):
				os.makedirs(a)

	model_names = ['M2','M3','M5']
	for p in pathways:
		for n in analysis_names:
			models_dirs = []
			elem_1 = './5_Data_Analysis/'+p+'/'+n+'/M2'
			elem_2 = './5_Data_Analysis/'+p+'/'+n+'/M3'
			elem_3 = './5_Data_Analysis/'+p+'/'+n+'/M5'
			models_dirs.append(elem_1)
			models_dirs.append(elem_2)
			models_dirs.append(elem_3)
			for m in models_dirs:
				if not os.path.exists(m):
					os.makedirs(m)			

	for p in pathways:
		for model in model_names:
			additional_dirs = []
			elem_1 = './5_Data_Analysis/'+p+'/LinearRegression/'+model+'/Coefficients'
			elem_2 = './5_Data_Analysis/'+p+'/LinearRegression/'+model+'/Confidence Intervals'
			elem_3 = './5_Data_Analysis/'+p+'/LinearRegression/'+model+'/Correlation Matrix'
			additional_dirs.append(elem_1)
			additional_dirs.append(elem_2)
			additional_dirs.append(elem_3)
			for ad in additional_dirs:
				if not os.path.exists(ad):
					os.makedirs(ad)