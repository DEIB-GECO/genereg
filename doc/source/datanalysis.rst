Data Analysis
============================================
In this phase all the data arranged in the previous phases are processed according to a suitable machine learning algorithm. In order to overcome potential issues related to overfitting and high computational complexity, the procedure is broken down into sub-processes with increasing complexity and a linear regression model is built for each target gene and for each one of the data matrixes previously defined.

Each model is built subsequently executing the following two steps:

	* a **FEATURE SELECTION** process that allows you to select the best subset of features for the model gene, removing inputs that are known to be non-significant for the regulation of the model gene expression beforehand. In particular, for each gene of interest and for each input matrix, the feature selection process is performed five times: in order to reduce the bias, the set of TCGA data samples is randomly split into five, possibly equal, groups of samples which are used to create five different testing sets (this partition is indeed performed only once at the beginning of the data analysis, and then the same five subsets of samples are used for processing all the genes. Therefore, feature selection is performed five times, one for each generated testing set (using the remaining samples as training set), according to a *k-fold cross-validation* process, setting *k=5*. The intersection of the five sets of extracted features is computed to obtain the final selected features for the current gene in the current matrix. Two different types of feature selection are available:
	
		* **foreward feature selection**, for which three different methods can be implemented by the library, according to how features selected in the previous steps are treated:
		
			* the default incremental method (**type='ffs_default'**), which implements the forewad feature selection retrieving in each step (i.e. for each data matrix) only the features extracted as most relevant in the previous step and re-evaluating them with all the new features in the current matrix in a new feature selection process;
			
			* an incremental method without the re-evaluation of previous features (**type='ffs_no_reval'**), which implements the forewad feature selection considering at each step (i.e. for each data matrix) only the new features added to the current matrix w.r.t. the previous one, maintaining the features selected in the previous step as fixed and without re-evaluating them in the new feature selection process (this means that the set of features selected for each data matrix will surely contain all the features selected for the previous one, along with one or more new relevant features);
			
			* a complete method (**type='all'**)which implements the foreward feature selection considering all the features in the current matrix , regardless of the features selected in the previous step. This method is useful when you want to analyze all the potential regulatory features all together, without any previous pre-selection process (e.g. this is extremely useful if applied to M5, which contains the whole set of candidate regulatory features). However, you have to pay attention to this type of selection process, because the computational complexity will inevitably be affected.

		* **Lasso**, where the features selected are the ones with a coefficient different from zero, according to the Lasso algorithm.

|

The following images illustrates the *k-fold cross-validation* process implemented during the *foreward feature selection* and the details of how the default incremental method works, respectively.

.. image:: images/fselection.png

|

.. image:: images/fselection2.png

|

and

	* a **REGRESSION** process that fits a linear model on each individual target gene and for each of its data matrixes, starting from the set of selected features as input and normalizing (i.e. *Z-score normalization*) all the data.

Here it is an example of the linear regression process for the sample gene TKT and the selection of its relevant features:

.. image:: images/lreg.png

|

However, since it is possible to introduce an alteration in the results due to the order in which features are included in the models, the regression procedure implemented here is logically based only on the biological knowledge and it therefore analyzes the full M2 matrix (which in turn contains M1) at first, including all at once in the model all the features whose exchanged order could cause different results in the analysis. Subsequently, data analysis on M3 and M5 is performed, finally executing only three regression processes for each target gene, assessing the impact of the model gene set first and then the influence of the other gene sets considered.

The set of fuctions used to perform the analysis is the following:

``feature_selection(gene_set, n_data_matrix, type)``

	The FEATURE_SELECTION operation executes the feature selection procedure per gene set and per data matrix, according to the specified type of selection: it takes as input the name of the gene set to consider and the number of the model to build (i.e., the number of the data matrix to consider) and performs the specified feature selection for all the genes of interest in the selected set. Results are exported locally either in Excel or text files.
	
	**Parameters:**
	
	* *gene_set*: the set of genes of interest to analyze
	
	* *n_data_matrix*: number identifying the data matrix to analyze (only 2,3 and 5 values are permitted)
	
	* *type*: the type of feature selection to perform (possibile values are {'ffs_default','ffs_no_reval','lasso','all'})
	
	Example::

		import genereg as gr
		gr.FeatureSelection.feature_selection(gene_set='DNA_REPAIR', n_data_matrix=2, type=ffs_default)
		gr.FeatureSelection.feature_selection(gene_set='DNA_REPAIR', n_data_matrix=3, type=ffs_default)
		gr.FeatureSelection.feature_selection(gene_set='DNA_REPAIR', n_data_matrix=5, type=ffs_default)

|

``linear_regression(gene_set, n_data_matrix)``

	The LINEAR_REGRESSION operation executes the linear regression analysis per gene set and per data matrix, considering as inputs of the model only the features selected during the previous feature selection procedure. Results are exported locally either in Excel or text files.
	
	**Parameters:**
	
	* *gene_set*: the set of genes of interest to analyze
	
	* *n_data_matrix*: number identifying the data matrix to analyze (only 2,3 and 5 values are permitted)
	
	Example::

		import genereg as gr
		gr.LinearRegression.linear_regression(gene_set='DNA_REPAIR', n_data_matrix=2)
		gr.LinearRegression.linear_regression(gene_set='DNA_REPAIR', n_data_matrix=3)
		gr.LinearRegression.linear_regression(gene_set='DNA_REPAIR', n_data_matrix=5)

|

``summarize_reg(gene_set, n_data_matrix)``

	The SUMMARIZE_REG operation summarizes all the data analysis results, by collecting them in convenient tables that exported locally in Excel files.
	
	**Parameters:**
	
	* *gene_set*: the set of genes of interest to summarize
	
	* *n_data_matrix*: number identifying the data matrix to summarize (only 2,3 and 5 values are permitted)
	
	Example::

		import genereg as gr
		gr.SummaryResults.summarize_reg(gene_set='DNA_REPAIR', n_data_matrix=2)
		gr.SummaryResults.summarize_reg(gene_set='DNA_REPAIR', n_data_matrix=3)
		gr.SummaryResults.summarize_reg(gene_set='DNA_REPAIR', n_data_matrix=5)

|

``summarize_r2(gene_set)``

	The SUMMARIZE_R2 operation summarizes R2 and Adjusted R2 scores for each target gene in each regression model, storing them locally in a single Excel file.
	
	**Parameters:**
	
	* *gene_set*: the set of genes of interest to summarize
	
	Example::

		import genereg as gr
		gr.SummaryResults.summarize_r2(gene_set='DNA_REPAIR')

|

``best_genes(gene_set)``

	The BEST_GENES operation collects the target genes with the best linear fit (Adjusted R2 >= 0.6) in the three regression models, storing them locally in a single Excel file.
	
	**Parameters:**
	
	* *gene_set*: the set of genes of interest to summarize
	
	Example::

		import genereg as gr
		gr.SummaryResults.best_genes(gene_set='DNA_REPAIR')
	