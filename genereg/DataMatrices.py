# coding: utf-8

# Import libraries
import pandas as pd
from pandas import ExcelWriter
import numpy as np
import pickle


def create_m1():

	"""
	The CREATE_M1 operation builds the first data matrix for each gene of interest, collecting the current gene expression and methylation values, along with the expression values of all the genes in the same gene set. One data matrix for each target gene is created and exported locally in as many Excel files as the considered genes; while the whole set of M1 matrixes is returned as a Python dictionary (dict_model_v1.p), where each target gene (set as key) is associated to a Pandas dataframe containing M1 data of interest (set as value). 

	:return: a Python dictionary
	
	Example::
	
		import genereg as gr
		m1_dict = gr.DataMatrixes.create_m1()
	"""


	# Load input data:
	
	# Genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})

	# Methylation values for genes of interest
	methyl_df = pd.read_excel('./3_TCGA_Data/Methylation/Methylation_Values.xlsx',sheetname='Sheet1',header=0)

	# Gene expression values for genes of interest
	expr_interest_df = pd.read_excel('./3_TCGA_Data/Gene_Expression/Gene_Expression-InterestGenes.xlsx',sheetname='Sheet1',header=0)


	# Create a list containing the Gene Symbols of the genes of interest
	gene_interest_SYMs = []
	for i, r in EntrezConversion_df.iterrows():
		sym = r['GENE_SYMBOL']
		if sym not in gene_interest_SYMs:
			gene_interest_SYMs.append(sym)

	# Get the TCGA aliquots 
	aliquots = []
	for i, r in methyl_df.iterrows():
		if i != 'ENTREZ_GENE_ID':
			aliquots.append(i)

		
	# Create a dictionary where, for each gene of interest set as key (the model gene), we have a dataframe representing the model (matrix of data) of that gene.
	# This model the expression and methylation values of the model gene in the first and second columns, and the expression of all the genes that belong to the
	# model gene set in the other columns, while the different TCGA aliquots are the indexes of the rows.
	dict_model_v1 = {}

	# Define the variables we need for the computation 
	model_gene_pathways = []  # list of the gene sets the model gene belongs to
	same_pathway_genes = []   # list of the symbols of the genes belonging to the same gene sets as the model gene
	df_columns = []           # list of the model columns names

	# Execute the following code for each gene of interest
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene  # get the Gene Symbol of the current gene
		
		# Get the gene sets of the model gene
		for i, r in EntrezConversion_df.iterrows():
			sym = r['GENE_SYMBOL']
			if sym == model_gene_SYM:
				p = r['GENE_SET']
				model_gene_pathways.append(p)
			
		# Get the genes of interest belonging to the model gene set
		for i, r in EntrezConversion_df.iterrows():
			path = r['GENE_SET']
			if path in model_gene_pathways:
				symbol = r['GENE_SYMBOL']
				if symbol != model_gene_SYM:
					same_pathway_genes.append(symbol)   
		
		# Define the columns of the model gene matrix of data
		df_columns.append('EXPRESSION ('+model_gene_SYM+')')  # the second column contains the expression of the model gene
		df_columns.append('METHYLATION ('+model_gene_SYM+')') # the third column contains the methylation of the model gene
		for g in same_pathway_genes:  
			df_columns.append(g)  # we have a column for each gene in the same gene set of the model gene
		
		# In correspondence of the model gene key in the dictionary,
		# set its model as value, with the proper indexes and column names 
		dict_model_v1[model_gene_SYM] = pd.DataFrame(index = aliquots, columns = df_columns)
		
		# Reset the variables for the next iteration on the next gene of interest
		model_gene_pathways = []
		same_pathway_genes = []
		df_columns = []

	# Fill the models for each gene of interest
	for gene, matrix in dict_model_v1.items():
		
		first_col = 'EXPRESSION ('+gene+')'
		second_col = 'METHYLATION ('+gene+')'
		
		# Add the expression and methylation values of each model gene and for each TCGA aliquot
		for index, row in matrix.iterrows():
			model_expr = expr_interest_df.get_value(index,gene) # get the expression
			model_methyl = methyl_df.get_value(index,gene)  # get the mathylation value
			# set the two values in the correct cell of the matrix
			matrix.set_value(index,first_col,model_expr)
			matrix.set_value(index,second_col,model_methyl)

		# Add the expression values for all the other genes belonging to the same gene set of the model gene
		for index, row in matrix.iterrows():
			for column_name, values in matrix.iteritems(): # iterate along the columns of the dataframe
				# skip the first two columns and add the proper values
				if (column_name != first_col) and (column_name != second_col): 
					expr = expr_interest_df.get_value(index,column_name)
					matrix.set_value(index,column_name,expr)

					
	# Export the dictionary into a pickle file in order to be able to import it back and use it to progressively build the next models for the genes of interest, adding further information.
	pickle.dump(dict_model_v1, open('./4_Data_Matrix_Construction/Model1/dict_model_v1.p', 'wb'))

	# Export the models as .xlsx files
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene
		pathway = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'GENE_SET'].iloc[0]
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'ENTREZ_GENE_ID'].iloc[0]
		file_name = 'Gene_'+gene_ID+'_['+model_gene_SYM+']'+'_('+pathway+')-Model_v1.xlsx'

		writer = ExcelWriter('./4_Data_Matrix_Construction/Model1/'+file_name)
		output_df =  dict_model_v1[model_gene_SYM]
		output_df.to_excel(writer,'Sheet1')
		writer.save()


	# Handle genes belonging to multiple gene sets
	multiple_pathway_genes = []
	n = EntrezConversion_df['GENE_SYMBOL'].value_counts()
	for i, v in n.items():
		if v > 1 :
			multiple_pathway_genes.append(i)
    
	for g in multiple_pathway_genes:
		filtered_df = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g]
		pathways = (filtered_df.GENE_SET.unique()).tolist()
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g, 'ENTREZ_GENE_ID'].iloc[0]
		
		for p in pathways:
			current_pathway_model =  dict_model_v1[g].copy()

			# Extract the genes of interest in the current gene set
			current_pathway_genes = []
			for i, r in EntrezConversion_df.iterrows():
				sym = r['GENE_SYMBOL']
				path = r['GENE_SET']
				if path == p:
					current_pathway_genes.append(sym)
					
			# Extract list of columns in the full model
			all_columns = []
			for column_name, values in current_pathway_model.iteritems():
				if (column_name != 'EXPRESSION ('+g+')') and (column_name != 'METHYLATION ('+g+')'): 
					all_columns.append(column_name)
					
			# Extract the columns to remove form the model
			other_pathway_genes = list(set(all_columns) - set(current_pathway_genes))

			for i in other_pathway_genes:
				if (i != g):
					current_pathway_model.drop(i, axis=1, inplace=True)

			writer = ExcelWriter('./4_Data_Matrix_Construction/Model1/Gene_'+gene_ID+'_['+g+']_('+p+')-Model_v1.xlsx')
			current_pathway_model.to_excel(writer,'Sheet1')
			writer.save()
	
	return dict_model_v1


def create_m2():

	"""
	The CREATE_M2 operation builds the second data matrix for each gene of interest, adding to the first matrix data about the expression of candidate regulatory genes of each gene of interest. One data matrix for each target gene is created and exported locally in as many Excel files as the considered genes; while the whole set of M2 matrixes is returned as a Python dictionary (dict_model_v2.p), where each target gene (set as key) is associated to a Pandas dataframe containing M2 data of interest (set as value). 

	:return: a Python dictionary
	
	Example::
	
		import genereg as gr
		m2_dict = gr.DataMatrixes.create_m2()
	"""

	
	# Load input data:

	# Genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})

	# Models_v1 of genes of interest
	dict_model_v1 = pickle.load(open('./4_Data_Matrix_Construction/Model1/dict_model_v1.p', 'rb'))

	# Distinct regulatory genes for each gene of interest
	dict_RegulGenes = pickle.load(open('./2_Regulatory_Genes/dict_RegulGenes.p', 'rb'))

	# Gene expression values for regulatory genes
	expr_regulatory_df = pd.read_excel('./3_TCGA_Data/Gene_Expression/Gene_Expression-RegulatoryGenes.xlsx',sheetname='Sheet1',header=0)


	# Create a list containing the Gene Symbols of the genes of interest
	gene_interest_SYMs = []
	for i, r in EntrezConversion_df.iterrows():
		sym = r['GENE_SYMBOL']
		if sym not in gene_interest_SYMs:
			gene_interest_SYMs.append(sym)

	# Get the TCGA aliquots 
	aliquots = []
	for i, r in expr_regulatory_df.iterrows():
		if i != 'ENTREZ_GENE_ID':
			aliquots.append(i)


	# Create a dictionary where, for each gene of interest set as key (the model gene), we have a dataframe representing the model (matrix of data) of that gene.
	# This model contains all the information in the first model, plus additional columns with the expression of the regulatory genes for each model gene,
	# while the different TCGA aliquots are the indexes of the rows
	dict_model_v2 = {}

	# Define the variables we need for the computation 
	model_gene_RegulGenes_SYM = []  # list of gene symbols for the regulatory genes of the model gene
	new_columns = []                # list of the new columns names to be added to the model

	# Execute the following code for each gene of interest
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene  # get the Gene Symbol of the current gene
		
		# Get the list of regulatory genes for the model gene
		model_gene_RegulGenes_SYM = dict_RegulGenes[model_gene_SYM]
		
		# Get the first model for the current gene (model_v1)
		model_1_df = dict_model_v1[model_gene_SYM]
		
		# Identify the new columns to be added to the matrix:
		# in this case they are the columns corresponding to regulatory genes of the model gene
		# (be careful not to have duplicated columns, so add only the symbols of the genes
		# that are not already contained in the previous model)
		old_columns = list(model_1_df.columns.values)
		for g in model_gene_RegulGenes_SYM:
			if g not in old_columns:
				new_columns.append(g)
		
		# Create the new part of the model to add
		new_df = pd.DataFrame(index = aliquots, columns = new_columns)

		# Add the expression values for all the new regulatory genes and for each TCGA aliquot
		for index, row in new_df.iterrows():
			for column_name, values in new_df.iteritems(): # iterate along the columns of the dataframe
				expr = expr_regulatory_df.get_value(index,column_name)
				new_df.set_value(index,column_name,expr)  
		
		# Join the two dataframes and create the new model (model_v2)
		model_2_df = model_1_df.join(new_df)    

		# Set the new model in correspondence of the correct model gene key in the new dictionary
		dict_model_v2[model_gene_SYM] = model_2_df

		# Reset the variables for the next iteration on the next gene of interest
		model_gene_RegulGenes_SYM = []
		new_columns = []


	# Check if some genes of interest have their own as candidate regulatory genes. If so, remove that column from the matrix
	for gene in gene_interest_SYMs:
		data_matrix = dict_model_v2[gene]
		matrix_cols = list(data_matrix.columns.values)
		if gene in matrix_cols:
			data_matrix.drop(gene, axis=1, inplace=True)
        

	# Export the dictionary into a pickle file in order to be able to import it back and use it to progressively build the next models for the genes of interest, adding further information
	pickle.dump(dict_model_v2, open('./4_Data_Matrix_Construction/Model2/dict_model_v2.p', 'wb'))

	# Export the models as .xlsx files
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene
		pathway = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'GENE_SET'].iloc[0]
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'ENTREZ_GENE_ID'].iloc[0]
		file_name = 'Gene_'+gene_ID+'_['+model_gene_SYM+']'+'_('+pathway+')-Model_v2.xlsx'

		writer = ExcelWriter('./4_Data_Matrix_Construction/Model2/'+file_name)
		output_df =  dict_model_v2[model_gene_SYM]
		output_df.to_excel(writer,'Sheet1')
		writer.save()

	
	# Handle genes belonging to multiple gene sets
	multiple_pathway_genes = []
	n = EntrezConversion_df['GENE_SYMBOL'].value_counts()
	for i, v in n.items():
		if v > 1 :
			multiple_pathway_genes.append(i)
    
	for g in multiple_pathway_genes:
		filtered_df = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g]
		pathways = (filtered_df.GENE_SET.unique()).tolist()
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g, 'ENTREZ_GENE_ID'].iloc[0]
		
		for p in pathways:
			# Import the 'model_v1' matrix for the current gene
			current_pathway_model = pd.read_excel('./4_Data_Matrix_Construction/Model1/Gene_'+gene_ID+'_['+g+']_('+p+')-Model_v1.xlsx',sheetname='Sheet1',header=0)

			# Get the list of regulatory genes for the model gene
			current_gene_RegulGenes_SYM = dict_RegulGenes[g]
            
			# Create the M2 model for the current gene in the current gene set, identifying the new columns to be added to the matrix
			current_pathway_new_columns = []
			current_pathway_old_columns = list(current_pathway_model.columns.values)
			for gene in current_gene_RegulGenes_SYM:
				if gene not in current_pathway_old_columns:
					current_pathway_new_columns.append(gene)
    
			# Create the new part of the model to add
			current_pathway_new_df = pd.DataFrame(index = aliquots, columns = current_pathway_new_columns)

			# Add the expression values for all the new regulatory genes and for each TCGA aliquot
			for index, row in current_pathway_new_df.iterrows():
				for column_name, values in current_pathway_new_df.iteritems(): # iterate along the columns of the dataframe
					expr = expr_regulatory_df.get_value(index,column_name)
					current_pathway_new_df.set_value(index,column_name,expr)  
		
			# Join the two dataframes and create the new model (model_v2)
			current_pathway_model_2_df = current_pathway_model.join(current_pathway_new_df)   

			# Check if some genes of interest have their own as candidate regulatory genes. If so, remove that column from the matrix
			current_pathway_matrix_cols = list(current_pathway_model_2_df.columns.values)
			if g in current_pathway_matrix_cols:
				current_pathway_model_2_df.drop(g, axis=1, inplace=True)

			writer = ExcelWriter('./4_Data_Matrix_Construction/Model2/Gene_'+gene_ID+'_['+g+']_('+p+')-Model_v2.xlsx')
			current_pathway_model_2_df.to_excel(writer,'Sheet1')
			writer.save()
	
	return dict_model_v2
	

def create_m3():

	"""
	The CREATE_M3 operation builds the third data matrix for the analysis for each gene of interest, adding to the second matrix data about the expression of candidate regulatory genes of genes of interest belonging to the same gene set of the model gene. One data matrix for each target gene is created and exported locally in as many Excel files as the considered genes; while the whole set of M3 matrixes is returned as a Python dictionary (dict_model_v3.p), where each target gene (set as key) is associated to a Pandas dataframe containing M3 data of interest (set as value). 

	:return: a Python dictionary
	
	Example::
	
		import genereg as gr
		m3_dict = gr.DataMatrixes.create_m3()
	"""


	# Load input data:

	# Genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})

	# Models_v2 of genes of interest
	dict_model_v2 = pickle.load(open('./4_Data_Matrix_Construction/Model2/dict_model_v2.p', 'rb'))

	# Distinct regulatory genes for each gene of interest
	dict_RegulGenes = pickle.load(open('./2_Regulatory_Genes/dict_RegulGenes.p', 'rb'))

	# Gene expression values for regulatory genes
	expr_regulatory_df = pd.read_excel('./3_TCGA_Data/Gene_Expression/Gene_Expression-RegulatoryGenes.xlsx',sheetname='Sheet1',header=0)


	# Create a list containing the Gene Symbols of the genes of interest
	gene_interest_SYMs = []
	for i, r in EntrezConversion_df.iterrows():
		sym = r['GENE_SYMBOL']
		if sym not in gene_interest_SYMs:
			gene_interest_SYMs.append(sym)

	# Get the TCGA aliquots 
	aliquots = []
	for i, r in expr_regulatory_df.iterrows():
		if i != 'ENTREZ_GENE_ID':
			aliquots.append(i)


	# Create a dictionary where, for each gene of interest set as key (the model gene), we have a dataframe representing the model (matrix of data) of that gene.
	# This model contains all the information in the second model, plus additional columns with the expression of the regulatory genes for each one of the genes belonging to the model gene set,
	# while the different TCGA aliquots are the indexes of the rows
	dict_model_v3 = {}

	# Define the variables we need for the computation 
	model_gene_pathways = []                # list of the gene sets the model gene belongs to
	same_pathway_genes = []                 # list of the symbols of the genes belonging to the same gene sets as the model gene
	same_pathway_genes_RegulGenes_SYM = []  # list of gene symbols for the regulatory genes of the genes in the same gene set
	new_columns = []                        # list of the new columns names to be added to the model

	# Execute the following code for each gene of interest
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene  # get the Gene Symbol of the current gene
		
		# Get the gene sets of the model gene
		for i, r in EntrezConversion_df.iterrows():
			sym = r['GENE_SYMBOL']
			if sym == model_gene_SYM:
				p = r['GENE_SET']
				model_gene_pathways.append(p)
			
		# Get the genes of interest belonging to the model gene sets
		for i, r in EntrezConversion_df.iterrows():
			path = r['GENE_SET']
			if path in model_gene_pathways:
				symbol = r['GENE_SYMBOL']
				if symbol != model_gene_SYM:
					same_pathway_genes.append(symbol)  
		
		# Get the list of regulatory genes for each one of the genes belonging to the same gene sets of the model gene
		for elem in same_pathway_genes:
			elem_regulatory_genes = dict_RegulGenes[elem]
			same_pathway_genes_RegulGenes_SYM = same_pathway_genes_RegulGenes_SYM + elem_regulatory_genes
		same_pathway_genes_RegulGenes_SYM = list(set(same_pathway_genes_RegulGenes_SYM)) # keep only distinct regulatory genes
		
		# Get the second model for the current gene (model_v2)
		model_2_df = dict_model_v2[model_gene_SYM]
			
		# Identify the new columns to be added to the matrix:
		# in this case they are the columns corresponding to regulatory genes of genes in the
		# same gene sets of our model gene
		# (be careful not to have duplicated columns, so add only the symbols of the genes
		#  that are not already contained in the previous model)
		old_columns = list(model_2_df.columns.values)
		for g in same_pathway_genes_RegulGenes_SYM:
			if g not in old_columns:
				new_columns.append(g)
		
		# Create the new part of the model to add
		new_df = pd.DataFrame(index = aliquots, columns = new_columns)

		# Add the expression values for all the new regulatory genes and for each TCGA aliquot
		for index, row in new_df.iterrows():
			for column_name, values in new_df.iteritems(): # iterate along the columns of the dataframe
				expr = expr_regulatory_df.get_value(index,column_name)
				new_df.set_value(index,column_name,expr)  
			
		# Join the two dataframes and create the new model (model_v3)
		model_3_df = model_2_df.join(new_df)    

		# Set the new model in correspondence of the correct model gene key in the new dictionary
		dict_model_v3[model_gene_SYM] = model_3_df

		# Reset the variables for the next iteration on the next gene of interest
		model_gene_pathways = []
		same_pathway_genes = []
		same_pathway_genes_RegulGenes_SYM = []
		new_columns = []
		
		
	# Remove duplicate columns of the model gene
	for gene in gene_interest_SYMs:
		data_matrix = dict_model_v3[gene]
		matrix_cols = list(data_matrix.columns.values)
		if gene in matrix_cols:
			data_matrix.drop(gene, axis=1, inplace=True)

			
	# Export the dictionary into a pickle file in order to be able to import it back and use it to progressively build the next models for the genes of interest, adding further information
	pickle.dump(dict_model_v3, open('./4_Data_Matrix_Construction/Model3/dict_model_v3.p', 'wb'))

	# Export the models as .xlsx files
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene
		pathway = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'GENE_SET'].iloc[0]
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'ENTREZ_GENE_ID'].iloc[0]
		file_name = 'Gene_'+gene_ID+'_['+model_gene_SYM+']'+'_('+pathway+')-Model_v3.xlsx'

		writer = ExcelWriter('./4_Data_Matrix_Construction/Model3/'+file_name)
		output_df =  dict_model_v3[model_gene_SYM]
		output_df.to_excel(writer,'Sheet1')
		writer.save()

		
	# Handle genes belonging to multiple gene sets
	multiple_pathway_genes = []
	n = EntrezConversion_df['GENE_SYMBOL'].value_counts()
	for i, v in n.items():
		if v > 1 :
			multiple_pathway_genes.append(i)
    
	for g in multiple_pathway_genes:
		filtered_df = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g]
		pathways = (filtered_df.GENE_SET.unique()).tolist()
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g, 'ENTREZ_GENE_ID'].iloc[0]
		
		for p in pathways:
			# Import the 'model_v2' matrix for the current gene
			current_pathway_model = pd.read_excel('./4_Data_Matrix_Construction/Model2/Gene_'+gene_ID+'_['+g+']_('+p+')-Model_v2.xlsx',sheetname='Sheet1',header=0)	
	           
			# Current gene set model
			current_pathway_genes = []
			current_pathway_RegulGenes_SYM = []
			current_pathway_new_columns = []

			# Get the genes of interest belonging to the model gene set
			for i, r in EntrezConversion_df.iterrows():
				path = r['GENE_SET']
				if path == p:
					sym = r['GENE_SYMBOL']
					if sym != g:
						current_pathway_genes.append(sym)  

			# Get the list of regulatory genes for each one of the genes belonging to the same gene sets of the model gene
			for elem in current_pathway_genes:
				elem_regulatory_genes = dict_RegulGenes[elem]
				current_pathway_RegulGenes_SYM = current_pathway_RegulGenes_SYM + elem_regulatory_genes
			current_pathway_RegulGenes_SYM = list(set(current_pathway_RegulGenes_SYM)) # keep only distinct regulatory genes
    
			# Identify the new columns to be added to the matrix
			current_pathway_old_columns = list(current_pathway_model.columns.values)
			for gene in current_pathway_RegulGenes_SYM:
				if gene not in current_pathway_old_columns:
					current_pathway_new_columns.append(gene)

			# Create the new part of the model to add
			current_pathway_new_df = pd.DataFrame(index = aliquots, columns = current_pathway_new_columns)
				
			# Add the expression values for all the new regulatory genes and for each TCGA aliquot
			for index, row in current_pathway_new_df.iterrows():
				for column_name, values in current_pathway_new_df.iteritems():
					expr = expr_regulatory_df.get_value(index,column_name)
					current_pathway_new_df.set_value(index,column_name,expr)  
        
			# Join the two dataframes and create the new model (model_v3)
			current_pathway_model_3_df = current_pathway_model.join(current_pathway_new_df)    

			# Remove duplicate columns of the model gene
			current_pathway_matrix_cols = list(current_pathway_model_3_df.columns.values)
			if g in current_pathway_matrix_cols:
				current_pathway_model_3_df.drop(g, axis=1, inplace=True)

			writer = ExcelWriter('./4_Data_Matrix_Construction/Model3/Gene_'+gene_ID+'_['+g+']_('+p+')-Model_v3.xlsx')
			current_pathway_model_3_df.to_excel(writer,'Sheet1')
			writer.save()	
	
	return dict_model_v3

	
def create_m4():

	"""
	The CREATE_M4 operation builds the fourth data matrix for the analysis for each gene of interest, adding to the third matrix data about the expression of genes of interest belonging to the other gene sets with respect ot the model gene. One data matrix for each target gene is created and exported locally in as many Excel files as the considered genes; while the whole set of M4 matrixes is returned as a Python dictionary (dict_model_v4.p), where each target gene (set as key) is associated to a Pandas dataframe containing M4 data of interest (set as value).

	:return: a Python dictionary
	
	Example::
	
		import genereg as gr
		m4_dict = gr.DataMatrixes.create_m4()
	"""


	# Load input data:
	
	# Genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})

	# Models_v3 of genes of interest
	dict_model_v3 = pickle.load(open('./4_Data_Matrix_Construction/Model3/dict_model_v3.p', 'rb'))

	# Gene expression values for genes of interest
	expr_interest_df = pd.read_excel('./3_TCGA_Data/Gene_Expression/Gene_Expression-InterestGenes.xlsx',sheetname='Sheet1',header=0)


	# Create a list containing the Gene Symbols of the genes of interest
	gene_interest_SYMs = []
	for i, r in EntrezConversion_df.iterrows():
		sym = r['GENE_SYMBOL']
		if sym not in gene_interest_SYMs:
			gene_interest_SYMs.append(sym)

	# Get the TCGA aliquots 
	aliquots = []
	for i, r in expr_interest_df.iterrows():
		if i != 'ENTREZ_GENE_ID':
			aliquots.append(i)
        

	# Create a dictionary where, for each gene of interest set as key (the model gene), we have a dataframe representing the model (matrix of data) of that gene.
	# This model contains all the information of the third model, plus additional columns with the expression of the genes of interest that belong to gene sets different from the ones of the model gene,
	# while the different TCGA aliquots are the indexes of the rows
	dict_model_v4 = {}

	# Define the variables we need for the computation 
	model_gene_pathways = []  # list of the gene sets the model gene belongs to
	other_pathway_genes = []  # list of the symbols of the genes belonging to different gene sets
	new_columns = []          # list of the new columns names to be added to the model

	# Execute the following code for each gene of interest
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene  # get the Gene Symbol of the current gene
		
		# Get the gene sets of the model gene
		for i, r in EntrezConversion_df.iterrows():
			sym = r['GENE_SYMBOL']
			if sym == model_gene_SYM:
				p = r['GENE_SET']
				model_gene_pathways.append(p)
			
		# Get the genes of interest belonging to other gene sets
		for i, r in EntrezConversion_df.iterrows():
			path = r['GENE_SET']
			if (path not in model_gene_pathways) and (path != 'GLUCOSE_METABOLISM'):
				symbol = r['GENE_SYMBOL']
				if symbol not in other_pathway_genes: # consider only once the genes belonging to multiple gene sets
					other_pathway_genes.append(symbol)  
		
		# Get the third model for the current gene (model_v3)
		model_3_df = dict_model_v3[model_gene_SYM]
		
		# Identify the new columns to be added to the matrix:
		# in this case they are the columns corresponding to genes of interest beloging to different
		# gene sets with respect to our model gene
		# (be careful not to have duplicated columns, so add only the symbols of the genes
		# that are not already contained in the previous model)
		old_columns = list(model_3_df.columns.values)
		for g in other_pathway_genes:
			if g not in old_columns:
				new_columns.append(g)
		
		# Create the new part of the model to add
		new_df = pd.DataFrame(index = aliquots, columns = new_columns)

		# Add the expression values for all the these genes of interest belonging to other gene sets and for each TCGA aliquot
		for index, row in new_df.iterrows():
			for column_name, values in new_df.iteritems(): # iterate along the columns of the dataframe
				expr = expr_interest_df.get_value(index,column_name)
				new_df.set_value(index,column_name,expr)  
			
		# Join the two dataframes and create the new model (model_v4)
		model_4_df = model_3_df.join(new_df)    

		# Set the new model in correspondence of the correct model gene key in the new dictionary
		dict_model_v4[model_gene_SYM] = model_4_df

		# Reset the variables for the next iteration on the next gene of interest
		model_gene_pathways = []
		other_pathway_genes = []
		new_columns = []
    
	
	# Export the dictionary into a pickle file in order to be able to import it back and use it to progressively build the next models for the genes of interest, adding further information
	pickle.dump(dict_model_v4, open('./4_Data_Matrix_Construction/Model4/dict_model_v4.p', 'wb'))

	# Export the models as .xlsx files
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene
		pathway = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'GENE_SET'].iloc[0]
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'ENTREZ_GENE_ID'].iloc[0]
		file_name = 'Gene_'+gene_ID+'_['+model_gene_SYM+']'+'_('+pathway+')-Model_v4.xlsx'

		writer = ExcelWriter('./4_Data_Matrix_Construction/Model4/'+file_name)
		output_df =  dict_model_v4[model_gene_SYM]
		output_df.to_excel(writer,'Sheet1')
		writer.save()


	# Handle genes belonging to multiple gene sets
	multiple_pathway_genes = []
	n = EntrezConversion_df['GENE_SYMBOL'].value_counts()
	for i, v in n.items():
		if v > 1 :
			multiple_pathway_genes.append(i)
    
	for g in multiple_pathway_genes:
		filtered_df = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g]
		pathways = (filtered_df.GENE_SET.unique()).tolist()
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g, 'ENTREZ_GENE_ID'].iloc[0]
		
		for p in pathways:
			# Import the 'model_v3' matrix for the current gene
			current_pathway_model = pd.read_excel('./4_Data_Matrix_Construction/Model3/Gene_'+gene_ID+'_['+g+']_('+p+')-Model_v3.xlsx',sheetname='Sheet1',header=0)		
		
			# Current gene set model
			current_pathway_other_genes = []
			current_pathway_new_columns = []

			# Get the genes of interest belonging to other gene sets
			for i, r in EntrezConversion_df.iterrows():
				path = r['GENE_SET']
				if (path != p):
					symbol = r['GENE_SYMBOL']
					if symbol != g:
						current_pathway_other_genes.append(symbol)  

			# Identify the new columns to be added to the matrix
			current_pathway_old_columns = list(current_pathway_model.columns.values)
			for gene in current_pathway_other_genes:
				if gene not in current_pathway_old_columns:
					current_pathway_new_columns.append(gene)

			# Create the new part of the model to add
			current_pathway_new_df = pd.DataFrame(index = aliquots, columns = current_pathway_new_columns)
				
			# Add the expression values for all the these genes of interest belonging to other gene sets and for each TCGA aliquot
			for index, row in current_pathway_new_df.iterrows():
				for column_name, values in current_pathway_new_df.iteritems():
					expr = expr_interest_df.get_value(index,column_name)
					current_pathway_new_df.set_value(index,column_name,expr)  
        
			# Join the two dataframes and create the new model (model_v4)
			current_pathway_model_4_df = current_pathway_model.join(current_pathway_new_df)    

			writer = ExcelWriter('./4_Data_Matrix_Construction/Model4/Gene_'+gene_ID+'_['+g+']_('+p+')-Model_v4.xlsx')
			current_pathway_model_4_df.to_excel(writer,'Sheet1')
			writer.save()	

	return dict_model_v4

	
def create_m5():

	"""
	The CREATE_M5 operation builds the fifth data matrix for the analysis for each gene of interest, adding to the fourth matrix data about the expression of candidate regulatory genes of genes of interest belonging to the other gene sets with respect to the model gene.. One data matrix for each target gene is created and exported locally in as many Excel files as the considered genes; while the whole set of M5 matrixes is returned as a Python dictionary (dict_model_v5.p), where each target gene (set as key) is associated to a Pandas dataframe containing M5 data of interest (set as value).

	:return: a Python dictionary
	
	Example::
	
		import genereg as gr
		m5_dict = gr.DataMatrixes.create_m5()
	"""

	
	# Load input data:

	# Genes of interest
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})

	# Models_v4 of genes of interest
	dict_model_v4 = pickle.load(open('./4_Data_Matrix_Construction/Model4/dict_model_v4.p', 'rb'))

	# Distinct regulatory genes for each gene of interest
	dict_RegulGenes = pickle.load(open('./2_Regulatory_Genes/dict_RegulGenes.p', 'rb'))

	# Gene expression values for regulatory genes
	expr_regulatory_df = pd.read_excel('./3_TCGA_Data/Gene_Expression/Gene_Expression-RegulatoryGenes.xlsx',sheetname='Sheet1',header=0)


	# Create a list containing the Gene Symbols of the genes of interest
	gene_interest_SYMs = []
	for i, r in EntrezConversion_df.iterrows():
		sym = r['GENE_SYMBOL']
		if sym not in gene_interest_SYMs:
			gene_interest_SYMs.append(sym)

	# Get the TCGA aliquots 
	aliquots = []
	for i, r in expr_regulatory_df.iterrows():
		if i != 'ENTREZ_GENE_ID':
			aliquots.append(i)


	# Create a dictionary where, for each gene of interest set as key (the model gene), we have a dataframe representing the model (matrix of data) of that gene.
	# This model contains all the information of the fourth model, plus additional columns with the expression of the regulatory genes for each gene of interest belonging to other gene sets with respect to the model gene,
	# while the different TCGA aliquots are the indexes of the rows
	dict_model_v5 = {}

	# Define the variables we need for the computation 
	model_gene_pathways = []  # list of the gene sets the model gene belongs to
	other_pathway_genes = []  # list of the gene symbol of the genes belonging to different gene sets
	other_pathway_genes_RegulGenes_SYM = []  # list of gene symbols for the regulatory genes of gene in other gene sets
	new_columns = []  # list of the new columns names to be added to the model

	# Execute the following code for each gene of interest
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene  # get the Gene Symbol of the current gene
		
		# Get the gene sets of the model gene
		for i, r in EntrezConversion_df.iterrows():
			sym = r['GENE_SYMBOL']
			if sym == model_gene_SYM:
				p = r['GENE_SET']
				model_gene_pathways.append(p)
			
		# Get the genes of interest belonging to other gene sets
		for i, r in EntrezConversion_df.iterrows():
			path = r['GENE_SET']
			if (path not in model_gene_pathways) and (path != 'GLUCOSE_METABOLISM'):
				symbol = r['GENE_SYMBOL']
				if symbol not in other_pathway_genes: # consider only once the genes belonging to multiple gene sets
					other_pathway_genes.append(symbol) 
		
		# Get the list of regulatory genes for the genes in the other gene sets
		for elem in other_pathway_genes:
			elem_regulatory_genes = dict_RegulGenes[elem]
			other_pathway_genes_RegulGenes_SYM = other_pathway_genes_RegulGenes_SYM + elem_regulatory_genes
		other_pathway_genes_RegulGenes_SYM = list(set(other_pathway_genes_RegulGenes_SYM)) # keep only distinct regulatory genes
				
		# Get the fourth model for the current gene (model_v4)
		model_4_df = dict_model_v4[model_gene_SYM]

		# Identify the new columns to be added to the matrix:
		# in this case they are the columns corresponding to regulatory genes of genes in other gene sets
		# (be careful not to have duplicated columns, so add only the symbols of the genes
		# that are not already contained in the previous model)
		old_columns = list(model_4_df.columns.values)
		for g in other_pathway_genes_RegulGenes_SYM:
			if g not in old_columns:
				new_columns.append(g)
		
		# Create the new part of the model to add
		new_df = pd.DataFrame(index = aliquots, columns = new_columns)

		# Add the expression values for all the new regulatory genes and for each TCGA aliquot
		for index, row in new_df.iterrows():
			for column_name, values in new_df.iteritems(): # iterate along the columns of the dataframe
				expr = expr_regulatory_df.get_value(index,column_name)
				new_df.set_value(index,column_name,expr)  
		
		# Join the two dataframes and create the new model (model_v5)
		model_5_df = model_4_df.join(new_df)    

		# Set the new model in correspondence of the correct model gene key in the new dictionary
		dict_model_v5[model_gene_SYM] = model_5_df

		# Reset the variables for the next iteration on the next gene of interest
		model_gene_pathways = []
		other_pathway_genes = []
		other_pathway_genes_RegulGenes_SYM = []
		new_columns = []
		
		
	# Remove duplicate columns of the model gene
	for gene in gene_interest_SYMs:
		data_matrix = dict_model_v5[gene]
		matrix_cols = list(data_matrix.columns.values)
		if gene in matrix_cols:
			data_matrix.drop(gene, axis=1, inplace=True)


	# Export the dictionary into a pickle file in order to be able to import it back and use it to progressively build the next models for the genes of interest, adding further information
	pickle.dump(dict_model_v5, open('./4_Data_Matrix_Construction/Model5/dict_model_v5.p', 'wb'))

	# Export the models as .xlsx files
	for gene in gene_interest_SYMs:
		
		model_gene_SYM = gene
		pathway = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'GENE_SET'].iloc[0]
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == model_gene_SYM, 'ENTREZ_GENE_ID'].iloc[0]
		file_name = 'Gene_'+gene_ID+'_['+model_gene_SYM+']'+'_('+pathway+')-Model_v5.xlsx'

		writer = ExcelWriter('./4_Data_Matrix_Construction/Model5/'+file_name)
		output_df =  dict_model_v5[model_gene_SYM]
		output_df.to_excel(writer,'Sheet1')
		writer.save()

				
	# Handle genes belonging to multiple gene sets
	multiple_pathway_genes = []
	n = EntrezConversion_df['GENE_SYMBOL'].value_counts()
	for i, v in n.items():
		if v > 1 :
			multiple_pathway_genes.append(i)
    
	for g in multiple_pathway_genes:
		filtered_df = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g]
		pathways = (filtered_df.GENE_SET.unique()).tolist()
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == g, 'ENTREZ_GENE_ID'].iloc[0]
		
		for p in pathways:
			# Import the 'model_v4' matrix for the current gene
			current_pathway_model = pd.read_excel('./4_Data_Matrix_Construction/Model4/Gene_'+gene_ID+'_['+g+']_('+p+')-Model_v4.xlsx',sheetname='Sheet1',header=0)	
	 
			# Current gene set model
			current_pathway_other_genes = []
			current_pathway_other_RegulGenes_SYM = []
			current_pathway_new_columns = []

			# Get the genes of interest belonging to other gene sets
			for i, r in EntrezConversion_df.iterrows():
				path = r['GENE_SET']
				if (path != p):
					sym = r['GENE_SYMBOL']
					if sym != g:
						current_pathway_other_genes.append(sym)  

			# Get the list of regulatory genes for each one of the genes belonging to other gene sets
			for elem in current_pathway_other_genes:
				elem_regulatory_genes = dict_RegulGenes[elem]
				current_pathway_other_RegulGenes_SYM = current_pathway_other_RegulGenes_SYM + elem_regulatory_genes
			current_pathway_other_RegulGenes_SYM = list(set(current_pathway_other_RegulGenes_SYM)) # keep only distinct regulatory genes

			# Identify the new columns to be added to the matrix
			current_pathway_old_columns = list(current_pathway_model.columns.values)
			for gene in current_pathway_other_RegulGenes_SYM :
				if gene not in current_pathway_old_columns:
					current_pathway_new_columns.append(gene)

			# Create the new part of the model to add
			current_pathway_new_df = pd.DataFrame(index = aliquots, columns = current_pathway_new_columns)
				
			# Add the expression values for all the new regulatory genes and for each TCGA aliquot
			for index, row in current_pathway_new_df.iterrows():
				for column_name, values in current_pathway_new_df.iteritems():
					expr = expr_regulatory_df.get_value(index,column_name)
					current_pathway_new_df.set_value(index,column_name,expr)  
        
			# Join the two dataframes and create the new model (model_v5)
			current_pathway_model_5_df = current_pathway_model.join(current_pathway_new_df)    

			# Remove duplicate columns of the model gene
			current_pathway_matrix_cols = list(current_pathway_model_5_df.columns.values)
			if g in current_pathway_matrix_cols:
				current_pathway_model_5_df.drop(g, axis=1, inplace=True)

			writer = ExcelWriter('./4_Data_Matrix_Construction/Model5/Gene_'+gene_ID+'_['+g+']_('+p+')-Model_v5.xlsx')
			current_pathway_model_5_df.to_excel(writer,'Sheet1')
			writer.save()

	return dict_model_v5