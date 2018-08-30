# coding: utf-8

# Import libraries
import pandas as pd
from pandas import ExcelWriter
import numpy as np
import pickle
from sklearn.linear_model import LinearRegression
from mlxtend.feature_selection import SequentialFeatureSelector as SFS


def forward_feature_selection(pathway, n_data_matrix):

	"""
	The FORWARD_FEATURE_SELECTION operation executes the forward feature selection procedure per pathway and per data matrix: it takes as input the name of the pathway to consider and the number of the model to build (i.e., the number of the data matrix to consider) and performs the fetaure selection for all the genes of interest in the selected pathway. Results are exported locally either in Excel or text files.

	:param pathway: genomic pathway to analyze
	:param n_data_matrix: number identifying the data matrix to analyze (only 2,3 and 5 values are permitted)
	"""
	
	
	# Define the model to create
	model = str(n_data_matrix)

	# Import the list of genes of interest and extract in a list the Gene Symbols of all the genes belonging to the current pathway
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

	SYMs_current_pathway = []
	for index, row in EntrezConversion_df.iterrows():
		sym = row['GENE_SYMBOL']
		path = row['PATHWAY']
		if path == pathway:
			SYMs_current_pathway.append(sym)

	# Create a dataframe to store the 5 R2 of the partitios
	r2_df = pd.DataFrame(index=SYMs_current_pathway, columns=[1,2,3,4,5])

	# Create a dataframe to store results of linear regression considering methylation and feature selection for each gene:
	# - total number of features in the current matrix M1 (exluding target expression)
	# - number of discarded features (NaN values)
	# - total number of features selected
	if model == '2':
		summary_results_df = pd.DataFrame(index=SYMs_current_pathway, columns=['TOT Inital N° Features','Discarded Features','N° Features Selected'])
	else:
		summary_results_df = pd.DataFrame(index=SYMs_current_pathway, columns=['TOT Inital N° Features','Discarded Features','Features Available for Selection','N° Features Selected'])

	
	for current_gene in SYMs_current_pathway:
		
		# Import the model corresponding to the current gene
		gene_ID = EntrezConversion_df.loc[EntrezConversion_df['GENE_SYMBOL'] == current_gene, 'ENTREZ_GENE_ID'].iloc[0]
		model_gene_df = pd.read_excel('./4_Data_Matrix_Construction/Model'+model+'/Gene '+gene_ID+' ['+current_gene+']'+' ('+pathway+') - Model_v'+model+'.xlsx',sheetname='Sheet1',header=0)

		tot_n_features = len(model_gene_df.columns) - 1 # all the columns, except for the expression of the model gene (i.e. target)
		# Count the number of columns that contains all NaN values and that will be discarded before the regression
		n_discarded_features = 0
		for col, val in model_gene_df.iteritems():
			s = model_gene_df[col]
			if s.isnull().values.all():
				n_discarded_features = n_discarded_features + 1
			
		# Remove NaN columns
		model_gene_df.dropna(axis=1, how='all', inplace=True)
			
		# Store the first results in the summary dataframe
		summary_results_df.set_value(current_gene, 'TOT Inital N° Features', tot_n_features)
		summary_results_df.set_value(current_gene, 'Discarded Features', n_discarded_features)     
		
		
		if (model == '3') or (model == '5'):
			# Load the list of features selected for the previous model
			previous_model = str(n_data_matrix-1)
			text_file = open('./5_Data_Analysis/'+pathway+'/FeatureSelection/M'+previous_model+'/Features - Gene '+gene_ID+' ['+current_gene+'].txt', 'r')
			prev_features = text_file.read().split('\n')
			prev_features.remove('')
			text_file.close()
		
			# Extract the features of the previous model that have not been selected by the feature selection and that we can remove from the current model before performing regression
			previous_col_names_to_delete = []
			for name, values in previous_model_df.iteritems():
				if ('EXPRESSION' not in name):
					if name not in prev_features:
						previous_col_names_to_delete.append(name)
		
			# Update the model keeping only the set of features we can select from
			current_model_col_names = set(list(model_gene_df.columns.values))
			previous_model_col_names = set(list(previous_model_df.columns.values))
			# if no new columns were added to the current matrix with respect to the previous one and no features were selected for the previous matrix, use again the whole matrix for the feature selection
			if (current_model_col_names.issubset(previous_model_col_names)) & (len(prev_features) == 0):
				model_gene_df = model_gene_df.copy()
			else:
				model_gene_df.drop(list(current_model_col_names & set(previous_col_names_to_delete)), axis=1, inplace=True)
			summary_results_df.set_value(current_gene, 'Features Available for Selection', (len(model_gene_df.columns)-1))
		
		
		# Import the dictionary containing the information needed to split the dataframe in five test sets
		dict_test_split = pickle.load(open('./5_Data_Analysis/dict_test_split.p', 'rb'))
		
		# Split the dataframe into five dataframes that will be used as test sets
		model_gene_df_test1 = model_gene_df.loc[dict_test_split['Test_1']]
		model_gene_df_test2 = model_gene_df.loc[dict_test_split['Test_2']]
		model_gene_df_test3 = model_gene_df.loc[dict_test_split['Test_3']]
		model_gene_df_test4 = model_gene_df.loc[dict_test_split['Test_4']]
		model_gene_df_test5 = model_gene_df.loc[dict_test_split['Test_5']]

		# Define the corresponding five dataframes to be used as training sets
		model_gene_df_train1 = model_gene_df[~model_gene_df.index.isin(model_gene_df_test1.index)]
		model_gene_df_train2 = model_gene_df[~model_gene_df.index.isin(model_gene_df_test2.index)]
		model_gene_df_train3 = model_gene_df[~model_gene_df.index.isin(model_gene_df_test3.index)]
		model_gene_df_train4 = model_gene_df[~model_gene_df.index.isin(model_gene_df_test4.index)]
		model_gene_df_train5 = model_gene_df[~model_gene_df.index.isin(model_gene_df_test5.index)]

		# Now, execute the feature selection five times, each time considering one of the five dataframes as test set.
		# MULTIPLE LINEAR REGRESSION: prediction of gene expression based on the expression of other specific genes 
		
		
		# CASE 1 ---------------------------------------------------------------------------------------------------------
		# Define the parameters:
		case = 1
		model_gene_df_train = model_gene_df_train1.copy()
		model_gene_df_test = model_gene_df_test1.copy()

		# Define the features (predictors X) and the target (label y), together with training and testing sets:
		# features X: model gene methylation and other genes expression values
		# target y: model gene expression value
		X_train = np.array(model_gene_df_train.drop(['EXPRESSION ('+current_gene+')'],1))
		X_test = np.array(model_gene_df_test.drop(['EXPRESSION ('+current_gene+')'],1))
		y_train = np.array(model_gene_df_train['EXPRESSION ('+current_gene+')'])
		y_test = np.array(model_gene_df_test['EXPRESSION ('+current_gene+')'])

		# Reshape y
		y_train = y_train.reshape(-1, 1)
		y_test = y_test.reshape(-1, 1)

		# APPLY FEATURE SELECTION

		# Sequential Feature Selector performs forward fueature selection.
		# The function used is the following:
		#    SFS(estimator, k_features, forward, floating, scoring, cv, n_jobs, ...)
		# where  estimator = scikit-learn classifier or regressor
		#        k_features = number of features to select (int or tuple with a min and max value).
		#                     SFS will consider return any feature combination between min and max
		#                     that scored highest in cross-validation.
		#                     If 'best' is provided, the feature selector will return the feature subset with the best
		#                     cross-validation performance. If 'parsimonious' is provided as an argument,
		#                     the smallest feature subset that is within one standard error of the cross-validation
		#                     performance will be selected.
		#        forward = forward selection if true, backward selection otherwise
		#        floating = allows to implement SFFS or SBFS
		#        scoring = scoring metric
		#                  {accuracy, f1, precision, recall, roc_auc} for classifiers
		#                  {'mean_absolute_error', 'neg_mean_squared_error', 'median_absolute_error', 'r2'} for regressors
		#        cv = cross-validation generator (default: 5)
		#        n_jobs = number of CPUs to use for evaluating different feature subsets in parallel ('-1' means 'all CPUs')

		# Define the linear regression object
		lr = LinearRegression()

		# Count the total number of features
		tot_N_features = int(X_train.shape[1]) 

		# Perform feature selection
		sfs = SFS(lr, k_features='best', forward=True, floating=False, scoring='neg_mean_squared_error', cv=5)

		# Learn model from training data
		sfs = sfs.fit(X_train, y_train)

		# Get all the details of the forward fits:
		# 'get_metric_dict(confidence_interval=0.95)' returns a dictionary, where dictionary keys are the number
		# of iterations (number of feature subsets) and where the value for each key is a second dictionary.
		# The keys of this second dictionary are:
		#    'feature_idx': tuple of the indices of the feature subset
		#    'cv_scores': list with individual CV scores
		#    'avg_score': average of CV scores
		#    'std_dev': standard deviation of the CV score average
		#    'std_err': standard error of the CV score average
		#    'ci_bound': confidence interval bound of the CV score average (around the computed cross-validation scores)
		# and they each have a different value in each iteration.
		# So, the general struture of this dictionary is the following:
		# {Iteration_1 : {feature_idx: tuple_of_values, cv_scores: list_of_values, avg_score: value,...},
		#  Iteration_2 : {feature_idx: tuple_of_values, cv_scores: list_of_values, avg_score: value,...},
		#  Iteration_3 : {feature_idx: tuple_of_values, cv_scores: list_of_values, avg_score: value,...}, ...}
		result_dict = sfs.get_metric_dict()
		
		# Compute the mean of cross-validation scores
		mean_cv_scores = []
		for i in np.arange(1,tot_N_features+1): # values are generated within the interval [start, stop), including start but excluding stop
			# since cv_scores are negative numbers in the previous dictionary, I have to add a '-' to compute the mean
			mean_cv_scores.append(-np.mean(result_dict[i]['cv_scores']))
		
		# Get the number of features selected, which corresponds to the number of features selected
		# in correspondence of the minimum of the mean cross-validation scores
		idx_1 = np.argmin(mean_cv_scores)+1

		# Get the features indexes for the best forward fit and convert them to list
		feature_idx_1 = result_dict[idx_1]['feature_idx']
		selected_features_indexes = list(feature_idx_1)

		# Extract the names of these features
		X_df = model_gene_df.drop(['EXPRESSION ('+current_gene+')'],1)
		X_df_columns = list(X_df.columns)

		columns_selected_1 = []
		for index in selected_features_indexes:
			columns_selected_1.append(X_df_columns[index])

		# FIT THE NEW MODEL

		# Define the new training and test set, according to the features selected
		X_train_sfs = X_train[:, feature_idx_1]
		X_test_sfs = X_test[:, feature_idx_1]

		# Define the linear regression object and train the model using training sets
		lr.fit(X_train_sfs, y_train)

		# Make predictions
		y_predicted = lr.predict(X_test_sfs)

		# Set the expected results
		y_expected = y_test

		# Compute the values of R-squared
		train_R2_1 = lr.score(X_train_sfs, y_train)
		test_R2_1 = lr.score(X_test_sfs, y_test)
		
		r2_df.set_value(current_gene,case,test_R2_1)

			
		# CASE 2 ---------------------------------------------------------------------------------------------------------
		# Define the parameters:
		case = 2
		model_gene_df_train = model_gene_df_train2.copy()
		model_gene_df_test = model_gene_df_test2.copy()

		# Define the features (predictors X) and the target (label y), together with training and testing sets:
		X_train = np.array(model_gene_df_train.drop(['EXPRESSION ('+current_gene+')'],1))
		X_test = np.array(model_gene_df_test.drop(['EXPRESSION ('+current_gene+')'],1))
		y_train = np.array(model_gene_df_train['EXPRESSION ('+current_gene+')'])
		y_test = np.array(model_gene_df_test['EXPRESSION ('+current_gene+')'])

		# Reshape y
		y_train = y_train.reshape(-1, 1)
		y_test = y_test.reshape(-1, 1)

		# APPLY FEATURE SELECTION

		lr = LinearRegression()
		tot_N_features = int(X_train.shape[1]) 
		sfs = SFS(lr, k_features='best', forward=True, floating=False, scoring='neg_mean_squared_error', cv=5)
		sfs = sfs.fit(X_train, y_train)

		# Get all the details of the forward fits:
		result_dict = sfs.get_metric_dict()
		
		# Compute the mean of cross-validation scores
		mean_cv_scores = []
		for i in np.arange(1,tot_N_features+1): 
			mean_cv_scores.append(-np.mean(result_dict[i]['cv_scores']))

		# Get the number of features selected
		idx_2 = np.argmin(mean_cv_scores)+1

		# Get the features indexes for the best forward fit and convert them to list
		feature_idx_2 = result_dict[idx_2]['feature_idx']
		selected_features_indexes = list(feature_idx_2)

		# Extract the names of these features
		X_df = model_gene_df.drop(['EXPRESSION ('+current_gene+')'],1)
		X_df_columns = list(X_df.columns)

		columns_selected_2 = []
		for index in selected_features_indexes:
			columns_selected_2.append(X_df_columns[index])
			 
		# FIT THE NEW MODEL

		X_train_sfs = X_train[:, feature_idx_2]
		X_test_sfs = X_test[:, feature_idx_2]

		lr.fit(X_train_sfs, y_train)
		y_predicted = lr.predict(X_test_sfs)
		y_expected = y_test

		# Compute the values of R-squared
		train_R2_2 = lr.score(X_train_sfs, y_train)
		test_R2_2 = lr.score(X_test_sfs, y_test)
		
		r2_df.set_value(current_gene,case,test_R2_2)
		
		
		# CASE 3 ---------------------------------------------------------------------------------------------------------
		# Define the parameters:
		case = 3
		model_gene_df_train = model_gene_df_train3.copy()
		model_gene_df_test = model_gene_df_test3.copy()

		# Define the features (predictors X) and the target (label y), together with training and testing sets:
		X_train = np.array(model_gene_df_train.drop(['EXPRESSION ('+current_gene+')'],1))
		X_test = np.array(model_gene_df_test.drop(['EXPRESSION ('+current_gene+')'],1))
		y_train = np.array(model_gene_df_train['EXPRESSION ('+current_gene+')'])
		y_test = np.array(model_gene_df_test['EXPRESSION ('+current_gene+')'])

		# Reshape y
		y_train = y_train.reshape(-1, 1)
		y_test = y_test.reshape(-1, 1)

		# APPLY FEATURE SELECTION

		lr = LinearRegression()
		tot_N_features = int(X_train.shape[1]) 
		sfs = SFS(lr, k_features='best', forward=True, floating=False, scoring='neg_mean_squared_error', cv=5)
		sfs = sfs.fit(X_train, y_train)

		# Get all the details of the forward fits:
		result_dict = sfs.get_metric_dict()
		
		# Compute the mean of cross-validation scores
		mean_cv_scores = []
		for i in np.arange(1,tot_N_features+1): 
			mean_cv_scores.append(-np.mean(result_dict[i]['cv_scores']))

		# Get the number of features selected
		idx_3 = np.argmin(mean_cv_scores)+1

		# Get the features indexes for the best forward fit and convert them to list
		feature_idx_3 = result_dict[idx_3]['feature_idx']
		selected_features_indexes = list(feature_idx_3)

		# Extract the names of these features
		X_df = model_gene_df.drop(['EXPRESSION ('+current_gene+')'],1)
		X_df_columns = list(X_df.columns)

		columns_selected_3 = []
		for index in selected_features_indexes:
			columns_selected_3.append(X_df_columns[index])
				
		# FIT THE NEW MODEL

		X_train_sfs = X_train[:, feature_idx_3]
		X_test_sfs = X_test[:, feature_idx_3]

		lr.fit(X_train_sfs, y_train)
		y_predicted = lr.predict(X_test_sfs)
		y_expected = y_test

		# Compute the values of R-squared
		train_R2_3 = lr.score(X_train_sfs, y_train)
		test_R2_3 = lr.score(X_test_sfs, y_test)
		
		r2_df.set_value(current_gene,case,test_R2_3)

		
		# CASE 4 ---------------------------------------------------------------------------------------------------------
		# Define the parameters:
		case = 4
		model_gene_df_train = model_gene_df_train4.copy()
		model_gene_df_test = model_gene_df_test4.copy()

		# Define the features (predictors X) and the target (label y), together with training and testing sets:
		X_train = np.array(model_gene_df_train.drop(['EXPRESSION ('+current_gene+')'],1))
		X_test = np.array(model_gene_df_test.drop(['EXPRESSION ('+current_gene+')'],1))
		y_train = np.array(model_gene_df_train['EXPRESSION ('+current_gene+')'])
		y_test = np.array(model_gene_df_test['EXPRESSION ('+current_gene+')'])

		# Reshape y
		y_train = y_train.reshape(-1, 1)
		y_test = y_test.reshape(-1, 1)

		# APPLY FEATURE SELECTION

		lr = LinearRegression()
		tot_N_features = int(X_train.shape[1]) 
		sfs = SFS(lr, k_features='best', forward=True, floating=False, scoring='neg_mean_squared_error', cv=5)
		sfs = sfs.fit(X_train, y_train)

		# Get all the details of the forward fits:
		result_dict = sfs.get_metric_dict()

		# Compute the mean of cross-validation scores
		mean_cv_scores = []
		for i in np.arange(1,tot_N_features+1): 
			mean_cv_scores.append(-np.mean(result_dict[i]['cv_scores']))

		# Get the number of features selected
		idx_4 = np.argmin(mean_cv_scores)+1

		# Get the features indexes for the best forward fit and convert them to list
		feature_idx_4 = result_dict[idx_4]['feature_idx']
		selected_features_indexes = list(feature_idx_4)

		# Extract the names of these features
		X_df = model_gene_df.drop(['EXPRESSION ('+current_gene+')'],1)
		X_df_columns = list(X_df.columns)

		columns_selected_4 = []
		for index in selected_features_indexes:
			columns_selected_4.append(X_df_columns[index])
		
		# FIT THE NEW MODEL

		X_train_sfs = X_train[:, feature_idx_4]
		X_test_sfs = X_test[:, feature_idx_4]

		lr.fit(X_train_sfs, y_train)
		y_predicted = lr.predict(X_test_sfs)
		y_expected = y_test

		# Compute the values of R-squared
		train_R2_4 = lr.score(X_train_sfs, y_train)
		test_R2_4 = lr.score(X_test_sfs, y_test)
		
		r2_df.set_value(current_gene,case,test_R2_4)

		
		# CASE 5 ---------------------------------------------------------------------------------------------------------
		# Define the parameters:
		case = 5
		model_gene_df_train = model_gene_df_train5.copy()
		model_gene_df_test = model_gene_df_test5.copy()

		# Define the features (predictors X) and the target (label y), together with training and testing sets:
		X_train = np.array(model_gene_df_train.drop(['EXPRESSION ('+current_gene+')'],1))
		X_test = np.array(model_gene_df_test.drop(['EXPRESSION ('+current_gene+')'],1))
		y_train = np.array(model_gene_df_train['EXPRESSION ('+current_gene+')'])
		y_test = np.array(model_gene_df_test['EXPRESSION ('+current_gene+')'])

		# Reshape y
		y_train = y_train.reshape(-1, 1)
		y_test = y_test.reshape(-1, 1)

		# APPLY FEATURE SELECTION

		lr = LinearRegression()
		tot_N_features = int(X_train.shape[1]) 
		sfs = SFS(lr, k_features='best', forward=True, floating=False, scoring='neg_mean_squared_error', cv=5)
		sfs = sfs.fit(X_train, y_train)

		# Get all the details of the forward fits:
		result_dict = sfs.get_metric_dict()
		
		# Compute the mean of cross-validation scores
		mean_cv_scores = []
		for i in np.arange(1,tot_N_features+1): 
			mean_cv_scores.append(-np.mean(result_dict[i]['cv_scores']))

		# Get the number of features selected
		idx_5 = np.argmin(mean_cv_scores)+1

		# Get the features indexes for the best forward fit and convert them to list
		feature_idx_5 = result_dict[idx_5]['feature_idx']
		selected_features_indexes = list(feature_idx_5)

		# Extract the names of these features
		X_df = model_gene_df.drop(['EXPRESSION ('+current_gene+')'],1)
		X_df_columns = list(X_df.columns)

		columns_selected_5 = []
		for index in selected_features_indexes:
			columns_selected_5.append(X_df_columns[index])
				
		# FIT THE NEW MODEL

		X_train_sfs = X_train[:, feature_idx_5]
		X_test_sfs = X_test[:, feature_idx_5]

		lr.fit(X_train_sfs, y_train)
		y_predicted = lr.predict(X_test_sfs)
		y_expected = y_test

		# Compute the values of R-squared
		train_R2_5 = lr.score(X_train_sfs, y_train)
		test_R2_5 = lr.score(X_test_sfs, y_test)
		
		r2_df.set_value(current_gene,case,test_R2_5)

		
		# Compute the mean of the five training and test R2 scores
		train_R2 = (train_R2_1+train_R2_2+train_R2_3+train_R2_4+train_R2_5)/5
		test_R2 = (test_R2_1+test_R2_2+test_R2_3+test_R2_4+test_R2_5)/5

		# Take the names of the features selected in the five cases, create their intersection and export the list of columns selected in a .txt file
		features_intersection = list(set(columns_selected_1) & set(columns_selected_2) & set(columns_selected_3) & set(columns_selected_4) & set(columns_selected_5))
		with open ('./5_Data_Analysis/'+pathway+'/FeatureSelection/M'+model+'/Features - Gene '+gene_ID+' ['+current_gene+'].txt', 'w') as fp:
			features_intersection_sorted = sorted(features_intersection)
			for i in features_intersection_sorted:
				fp.write('%s\n' % i)
		summary_results_df.set_value(current_gene, 'N° Features Selected', len(features_intersection))  

	# Export the summary dataframe in an Excel file
	writer = ExcelWriter('./5_Data_Analysis/'+pathway+'/FeatureSelection/M'+model+'/Feature Selection SUMMARY.xlsx')
	summary_results_df.to_excel(writer,'Sheet1')
	writer.save()