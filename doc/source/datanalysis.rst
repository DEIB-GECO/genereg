Data Analysis
============================================
In this phase all the data arranged in the previous phases are processed according to a suitable machine learning algorithm. In order to overcome potential issues related to overfitting and high computational complexity, the procedure is broken down into sub-processes with increasing complexity and a linear regression model is built for each target gene and for each one of the data matrixes previously defined.

Each model is built subsequently executing the following two steps:

	* a **FEATURE SELECTION** process that allows selecting the best subset of features for the target gene, removing inputs that are known to be non-significant for the regulation of the target gene expression beforehand. In particular, for each target gene of interest and for each input matrix, the feature selection process is performed five times: in order to reduce the bias, the set of TCGA data samples is randomly split into five, possibly equal, groups of samples which are used to create five different testing sets (this partition is indeed performed only once at the beginning of the data analysis, and then the same five subsets of samples are used for processing all the genes). Therefore, feature selection is performed five times, one for each generated testing set (using the remaining samples as training set), according to a *k-fold cross-validation* process, setting *k=5*. The intersection of the five sets of extracted features is computed to obtain the final selected features for the target gene in the considered data matrix. Two different types of feature selections are available:
	
		* **forward feature selection**, for which three different methods are available in the library, according to how features selected in the previous steps are treated:
		
			* the default incremental method (**type='ffs_default'**), which implements the forward feature selection considering in each step (i.e. for each data matrix) only the features extracted as most relevant in the previous step and re-evaluating them with all the new features in the current matrix in a new feature selection process;
			
			* an incremental method without the re-evaluation of previous features (**type='ffs_no_reval'**), which implements the forward feature selection considering at each step (i.e. for each data matrix) only the new features added to the current matrix w.r.t. the previous one, maintaining the features selected in the previous step as fixed and without re-evaluating them in the new feature selection process (this means that the set of features selected for each data matrix will surely contain all the features selected for the previous data matrix, along with one or more new relevant features);
			
			* a complete method (**type='all'**) which implements the forward feature selection considering all the features in the current matrix, regardless of the features selected in the previous step. This method allows analyzing all the potential regulatory features all together, without any previous pre-selection process (e.g. in particular this can be useful if applied to M5, which contains the whole set of candidate regulatory features). However, pay attention to this type of selection process, because the computational complexity will inevitably increase, in case greatly.

		* **Lasso**, where the features selected are the ones with a coefficient different from zero, according to the Lasso algorithm.

|

The following images illustrate the *k-fold cross-validation* process implemented during the *forward feature selection* and the details of how the default incremental method works, respectively.

.. image:: images/fselection.png

|

.. image:: images/fselection2.png

|

and

	* a **REGRESSION** process that fits a linear model on each individual target gene and for each of its data matrixes, starting from the set of selected features as input and normalizing (i.e. *Z-score normalization*) all the data.

Here it is an example of the linear regression process for the sample gene TKT and the selection of its relevant features:

.. image:: images/lreg.png

|

However, an alteration in the results can be introduced due to the order in which features are included in the models (i.e. considered data matrixes). Thus, the regression procedure implemented is logically based only on the biological meaning of the considered features; it therefore analyzes the full M2 matrix (which contains M1) at first, including all at once the features whose exchanged order could cause different results in the analysis. Subsequently, data analysis on M3 and M5 is performed, finally executing only three regression processes for each target gene, assessing the impact of the genes in the target gene set first and then the influence of the genes in the other gene sets considered.

The set of functions used to perform the analysis is the following:

``feature_selection(gene_set, n_data_matrix, type)``

	The FEATURE_SELECTION operation executes the feature selection procedure per gene set and per data matrix, according to the specified type of selection: it takes as input the name of the gene set to consider and the number of the model to build (i.e., the number of the data matrix to consider) and performs the specified feature selection for all the genes of interest in the selected set. Results are exported locally in Excel and text files.
	
	**Parameters:**
	
	* *gene_set*: the name of the set of genes of interest to analyze
	
	* *n_data_matrix*: number identifying the data matrix to analyze (only 2, 3 and 5 values are permitted)
	
	* *type*: the type of feature selection to perform (possible values are {'ffs_default', 'ffs_no_reval', 'lasso', 'all'})
	
	**INPUT FILES:** Genes_of_Interest.xlsx from *./*, data matrixes from *./4_Data_Matrix_Construction/*
	
	**OUTPUT_FILES:** `Features - Gene ID [SYM].txt <https://raw.githubusercontent.com/Kia23/genereg/master/DATA/sample_files/Features%20-%20Gene%20672%20%5BBRCA1%5D.txt>`_, `Feature Selection SUMMARY.xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Feature%20Selection%20SUMMARY.xlsx>`_ (click on the files to download an example)
	
	Example::

		import genereg as gr
		gr.FeatureSelection.feature_selection(gene_set='DNA_REPAIR', n_data_matrix=2, type=ffs_default)
		gr.FeatureSelection.feature_selection(gene_set='DNA_REPAIR', n_data_matrix=3, type=ffs_default)
		gr.FeatureSelection.feature_selection(gene_set='DNA_REPAIR', n_data_matrix=5, type=ffs_default)

|

``linear_regression(gene_set, n_data_matrix)``

	The LINEAR_REGRESSION operation executes the linear regression analysis per gene set and per data matrix, considering as inputs of the model only the features selected during the previous feature selection procedure. Results are exported locally in Excel and text files.
	
	**Parameters:**
	
	* *gene_set*: the name of the set of genes of interest to analyze
	
	* *n_data_matrix*: number identifying the data matrix to analyze (only 2, 3 and 5 values are permitted)
	
	**INPUT FILES:** Genes_of_Interest.xlsx from *./*, data matrixes from *./4_Data_Matrix_Construction/*, features selected from *./5_Data_Analysis/.../FeatureSelection/M.../*
	
	**OUTPUT_FILES:** `LinReg Summary - Gene ID [SYM].txt <https://raw.githubusercontent.com/Kia23/genereg/master/DATA/sample_files/LinReg%20Summary%20-%20Gene%20672%20%5BBRCA1%5D.txt>`_, `Coefficients (model) - Gene Gene ID [SYM].xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Coefficients%20(M3)%20-%20Gene%20672%20%5BBRCA1%5D.xlsx>`_, `Confidence Intervals (model) - Gene Gene ID [SYM].xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Confidence%20Intervals%20(M3)%20-%20Gene%20672%20%5BBRCA1%5D.xlsx>`_, `Correlation Matrix (model) - Gene Gene ID [SYM].xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Correlation%20Matrix%20(M3)%20-%20Gene%20672%20%5BBRCA1%5D.xlsx>`_, `Linear Regression R2 SCORES.xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Linear%20Regression%20R2%20SCORES.xlsx>`_ (click on the files to download an example)
	
	Example::

		import genereg as gr
		gr.LinearRegression.linear_regression(gene_set='DNA_REPAIR', n_data_matrix=2)
		gr.LinearRegression.linear_regression(gene_set='DNA_REPAIR', n_data_matrix=3)
		gr.LinearRegression.linear_regression(gene_set='DNA_REPAIR', n_data_matrix=5)

|

``summarize_reg(gene_set, n_data_matrix)``

	The SUMMARIZE_REG operation summarizes all the data analysis results, by collecting them in convenient tables that are exported locally in Excel files.
	
	**Parameters:**
	
	* *gene_set*: the name of the set of genes of interest to summarize
	
	* *n_data_matrix*: number identifying the data matrix to summarize (only 2, 3 and 5 values are permitted)
	
	**INPUT FILES:** Genes_of_Interest.xlsx from *./*, dict_RegulGenes.p from *./2_Regulatory_Genes/*, ata matrixes from *./4_Data_Matrix_Construction/*, Feature Selection SUMMARY.xlsx from *./5_Data_Analysis/.../FeatureSelection/M.../*, Linear Regression R2 SCORES.xlsx from *./5_Data_Analysis/.../LinearRegression/M.../*
	
	**OUTPUT_FILES:** `Feature Selection and Linear Regression.xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Feature%20Selection%20and%20Linear%20Regression.xlsx>`_, `Relevant Features - Gene ID [SYM].xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Relevant%20Features%20-%20Gene%20672%20%5BBRCA1%5D.xlsx>`_, `Order of Features Selected.xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Order%20of%20Features%20Selected.xlsx>`_ (click on the files to download an example)
	
	Example::

		import genereg as gr
		gr.SummaryResults.summarize_reg(gene_set='DNA_REPAIR', n_data_matrix=2)
		gr.SummaryResults.summarize_reg(gene_set='DNA_REPAIR', n_data_matrix=3)
		gr.SummaryResults.summarize_reg(gene_set='DNA_REPAIR', n_data_matrix=5)

|

``summarize_r2(gene_set)``

	The SUMMARIZE_R2 operation summarizes R2 and Adjusted R2 scores for each target gene in each regression model, storing them locally in a single Excel file.
	
	**Parameters:**
	
	* *gene_set*: the name of the set of genes of interest to summarize
	
	**INPUT FILES:** Genes_of_Interest.xlsx from *./*, Feature Selection and Linear Regression.xlsx from *./5_Data_Analysis/.../*
	
	**OUTPUT_FILES:** `R2 and Adj.R2 Scores.xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/R2%20and%20Adj.R2%20Scores.xlsx>`_ (click on the file to download an example)
	
	Example::

		import genereg as gr
		gr.SummaryResults.summarize_r2(gene_set='DNA_REPAIR')

|

``best_genes(gene_set)``

	The BEST_GENES operation collects the target genes with the best linear fit (Adjusted R2 >= 0.6) in the three regression models built, storing them locally in a single Excel file.
	
	**Parameters:**
	
	* *gene_set*: the name of the set of genes of interest to summarize
	
	**INPUT FILES:** Genes_of_Interest.xlsx from *./*, R2 and Adj.R2 Scores.xlsx from *./5_Data_Analysis/.../*
	
	**OUTPUT_FILES:** `Best Genes.xlsx <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Best%20Genes.xlsx>`_ (click on the file to download an example)
	
	Example::

		import genereg as gr
		gr.SummaryResults.best_genes(gene_set='DNA_REPAIR')
	