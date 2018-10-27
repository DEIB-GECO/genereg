ENCODE Data Extraction
============================================
The first data extraction phase consists in selecting from ENCODE all the Transcription Factors (TFs) having binding sites located in the promoter regions of genes of interest. Then, starting from these TFs, candidate regulatory genes for each gene of interest are identified:

.. image:: images/TFs.png


The query for extracting TFs of interest is implemented according to the `PyGMQL <https://pygmql.readthedocs.io/en/latest/index.html>`_ syntax and data are retrieved from public datasets available on the `GMQL <http://gmql.eu/gmql-rest/>`_
system.

|

-------------------------------------
Extraction of Transcription Factors
-------------------------------------

``extract_tfs(cell_lines, gencode_version)``

	The EXTRACT_TFS operation extracts, from ChIP_Seq ENCODE experiments and for assembly GRCh38, the Transcription Factors of specific cell lines that bind to promoter regions of genes in the specified version of the GENCODE genomic annotations, filtering by 'conservative idr thresholded peaks' in order to extract higher quality region data, and removing negative audits in order to keep only high quality data samples. Intermediate result files are exported locally during the execution of the function, while the final set of transcription factors is returned as a Python dictionary (dict_GeneTF.p), where each target gene (set as key) is associated with the list of TFs binding to its promoters (set as value).
	
	**Parameters:**
	
	* *cell_lines*: a list of strings containing the names of the cell lines to analyze (it is possible to consider data from 1 up to 3 cell lines at the same time)
	
	* *gencode_version*: number representing the GENCODE genomic annotations to use (currently, for assembly GRCh38, versions 22, 24 and 27 can be used)
	
	**Return:** a Python dictionary
	
	Example::

		import genereg as gr
		tfs_dict = gr.TranscriptionFactors.extract_tfs(cell_lines=['K562','MCF7'], gencode_version=22)

|

-------------------------------------
Identification of Regulatory Genes
-------------------------------------

``extract_regulatory_genes()``

	The EXTRACT_REGULATORY_GENES operation extracts from the set of Transcription Factors associated with a gene, the list of its candidate regulatory genes, i.e., the genes that encode for those TFs. Intermediate results files are exported locally during the execution of the function, while the final set of transcription factors is returned as a Python dictionary (dict_RegulGenes.p), where each target gene (set as key) is associated with the list of its candidate regulatory genes (set as value).
	
	**Return:** a Python dictionary
	
	Example::

		import genereg as gr
		reg_genes_dict = gr.RegulatoryGenes.extract_regulatory_genes()

|

The following image explains how the extraction process for each target gene works and how its TFs and regulatory genes are identified:

.. image:: images/encode.png

