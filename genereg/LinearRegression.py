# coding: utf-8

# Import libraries
import pandas as pd
from pandas import ExcelWriter
import numpy as np
import math
import pickle
from sklearn import preprocessing
import statsmodels.api as sm


def linear_regression(gene_set, n_data_matrix):

	"""
	The LINEAR_REGRESSION operation executes the linear regression analysis per gene set and per data matrix, considering as inputs of the model only the features selected during the previous feature selection procedure. Results are exported locally either in Excel or text files.

	:param gene_set: the set of genes of interest to analyze
	:param n_data_matrix: number identifying the data matrix to analyze (only 2,3 and 5 values are permitted)
	
	Example::
	
		import genereg as gr
		gr.LinearRegression.linear_regression(gene_set='DNA_REPAIR', n_data_matrix=2)
		gr.LinearRegression.linear_regression(gene_set='DNA_REPAIR', n_data_matrix=3)
		gr.LinearRegression.linear_regression(gene_set='DNA_REPAIR', n_data_matrix=5)		
	"""

	# Check input parameters	
	if n_data_matrix not in [2, 3, 5]:
		raise ValueError('Data Matrix ERROR! Possible values: {2,3,5}')

	# Define the model to create
	model = str(n_data_matrix)

	# Import the list of genes of interest and extract in a list the Gene Symbols of all the genes belonging to the current gene set
	EntrezConversion_df = pd.read_excel('./Genes_of_Interest.xlsx',sheetname='Sheet1',header=0,converters={'GENE_SYMBOL':str,'ENTREZ_GENE_ID':str,'GENE_SET':str})

	SYMs_current_pathway = []
	for index, row in EntrezConversion_df.iterrows():
		sym = row['GENE_SYMBOL']
		path = row['GENE_SET']
		if path == gene_set:
			SYMs_current_pathway.append(sym)

	# Create a dataframe to store results of linear regression for each gene (i.e. R2 and Adjusted R2 of the linear regression)
	all_r2_df = pd.DataFrame(index=SYMs_current_pathway, columns=['Adj.R2','R2'])  


	for current_gene in SYMs_current_pathway:
		
		# Import the model corresponding to the current gene
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == current_gene, 'ENTREZ_GENE_ID'].iloc[0]
		model_gene_df = pd.read_excel('./4_Data_Matrix_Construction/Model'+model+'/Gene_'+gene_ID+'_['+current_gene+']'+'_('+gene_set+')-Model_v'+model+'.xlsx',sheetname='Sheet1',header=0)
		
		# Set all the unknown values (NaN) to zero
		model_gene_df = model_gene_df.fillna(0)
		
		
		# DATA STANDARDIZATION:
		# Normalize expression values (using the proper scaler from the sklearn library):
		# MinMaxScaler() normalizes values between 0 and 1
		# StandardScaler() performs Z-score normalization
		scaler = preprocessing.StandardScaler() # define the scaler

		# Define the dataframe to normalize
		to_normalize = model_gene_df.copy()
		matrix = to_normalize.values # convert into a numpy array

		# Normalize and convert back to pandas dataframe
		matrix_scaled = scaler.fit_transform(matrix)
		model_gene_df = pd.DataFrame(matrix_scaled, index=to_normalize.index, columns=to_normalize.columns)
		original_model_gene_df = model_gene_df.copy()
		
		
		# Load the list of features selected
		text_file = open('./5_Data_Analysis/'+gene_set+'/FeatureSelection/M'+model+'/Features-Gene_'+gene_ID+'_['+current_gene+'].txt', 'r')
		features_sel = text_file.read().split('\n')
		features_sel.remove('')
		
		# Consider only the features selected and perform a complete linear regression on the data
		if len(features_sel) > 0:

			# Filter the initial matrix of data extracting only the columns selected by feature selection
			cols_to_extract = ['EXPRESSION ('+current_gene+')']+features_sel
			model_gene_df_filtered = model_gene_df[cols_to_extract].copy()

		
			# PERFORM LINEAR REGRESSION

			# Define the features (predictors X) and the target (label y)
			X = model_gene_df_filtered.drop(['EXPRESSION ('+current_gene+')'],1)
			y = model_gene_df_filtered['EXPRESSION ('+current_gene+')']
		
			# Add an intercept to our model 
			X = sm.add_constant(X)

			# Define the linear regression object and fit the model
			lr_model = sm.OLS(y, X).fit()

			# Make predictions
			y_predicted = lr_model.predict(X)

			# Compute the error rate (Root Mean Square Error).
			# RMS Error measures the differences between predicted values and values actually observed
			rms = np.sqrt(np.mean((np.array(y_predicted) - np.array(y)) ** 2))
		
			# Export the summary statistics into a text file
			lr_summary = lr_model.summary()
			lr_summary_file = lr_summary.as_text()
			with open ('./5_Data_Analysis/'+gene_set+'/LinearRegression/M'+model+'/LinReg_Summary-Gene_'+gene_ID+'_['+current_gene+'].txt', 'w') as fp:
				fp.write(lr_summary_file)
				fp.write('\n\nRoot Mean Square Error: RMS = '+str(rms))  

			# Store the R-squared score in the summery dataframe
			all_r2_df.set_value(current_gene, 'Adj.R2', lr_model.rsquared_adj)
			all_r2_df.set_value(current_gene, 'R2', lr_model.rsquared)
		
			# Save the coefficients and the intercept of the model
			coeff = lr_model.params
			coeff_df = pd.DataFrame({'feature':coeff.index, 'coefficient':coeff.values})
			coeff_df = coeff_df.sort_values(by=['coefficient'], ascending=[False])
			coeff_df_ordered = coeff_df[['feature','coefficient']].copy()
			coeff_path = './5_Data_Analysis/'+gene_set+'/LinearRegression/M'+model+'/Coefficients/'
			writer = ExcelWriter(coeff_path+'Coefficients_(M'+model+')-Gene_'+gene_ID+'_['+current_gene+'].xlsx')
			coeff_df_ordered.to_excel(writer,'Sheet1',index=False)
			writer.save()
		   
			
			# Compute and export the confidence intervals for the model coefficients (default: 95%)
			CI_df = lr_model.conf_int()
			CI_df.rename(columns={0: 'min_bound', 1: 'max_bound'}, inplace=True)
		
			# Verify the relevance of each feature by checking if 0 is contained in the confidence interval
			for index, row in CI_df.iterrows():
				min_ci = row['min_bound']
				max_ci = row['max_bound']
				if (0 >= min_ci) and (0 <= max_ci):
					CI_df.set_value(index,'Significant Feature?','NO')
				else:
					CI_df.set_value(index,'Significant Feature?','YES')
		
			# Extract the probabilty that the value of the feature is 0
			p_val = lr_model.pvalues
			p_val_df = pd.DataFrame({'feature':p_val.index, 'p-value':p_val.values})
			for index, row in p_val_df.iterrows():
				feature_gene = row['feature']
				p = row['p-value']
				CI_df.set_value(feature_gene,'P',p)  
		
			ci_path = './5_Data_Analysis/'+gene_set+'/LinearRegression/M'+model+'/ConfidenceIntervals/'
			writer = ExcelWriter(ci_path+'Confidence_Intervals_(M'+model+')-Gene_'+gene_ID+'_['+current_gene+'].xlsx')
			CI_df.to_excel(writer,'Sheet1')
			writer.save()
		
		
			# Compute the correlation matrix between gene data
			if 'METHYLATION ('+current_gene+')' in cols_to_extract:
				model_gene_df_corr = model_gene_df_filtered.copy()
			else:
				cols_for_correlation = ['EXPRESSION ('+current_gene+')','METHYLATION ('+current_gene+')']+features_sel
				model_gene_df_corr = original_model_gene_df[cols_for_correlation].copy()
			corr_matrix = model_gene_df_corr.corr()    
			corr_path = './5_Data_Analysis/'+gene_set+'/LinearRegression/M'+model+'/CorrelationMatrix/'
			writer = ExcelWriter(corr_path+'Correlation_Matrix_(M'+model+')-Gene_'+gene_ID+'_['+current_gene+'].xlsx')
			corr_matrix.to_excel(writer,'Sheet1')
			writer.save()
		
		
	# Export the summary dataframe in an Excel file
	all_r2_df = all_r2_df.sort_values(by=['Adj.R2'], ascending=[False])
	writer = ExcelWriter('./5_Data_Analysis/'+gene_set+'/LinearRegression/M'+model+'/Linear_Regression_R2_SCORES.xlsx')
	all_r2_df.to_excel(writer,'Sheet1')
	writer.save()