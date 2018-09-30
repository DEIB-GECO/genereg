Initialization
============================================

The library automatically creates directories for storing all the intermediate and final results in the current workspace of the user (i.e. the directory Python it is executed in), in order to be able to use them during the execution, when needed.
For this reason, it is advised to execute Python direclty within the directory you want the results to be saved in (let us call it "*library workspace*") and it is highly recommended not to move any file until the whole process has been fully executed, that is all the functions of this library have been run and the whole analysis procedure has been perfomed.

------------------
Genes of Interest
------------------
The first thing you need to do is defining your genes of interest for the pathology you want to analyze: these are the target genes whose regulation systems will be assesed during the analysis. Genes belonging to one or multiple gene sets are supported. So, you can analyze genes from different gene sets of your interest, by simply specifying the name of each set, as explained below.

Here it is how this list of genes has to be structured:
    
	* create an Excel file (.xlsx) with 3 columns ('GENE_SYMBOL', 'ENTREZ_GENE_ID', 'GENE_SET') and one row for each target gene
	
		* GENE_SYMBOL is the official gene name

		* ENTREZ_GENE_ID is the numerical ID associated with the gene

		* GENE_SET is the name of the gene set the corresponding gene belongs to 

	* save this file as "*Genes_of_Interest.xlsx*" in the *library workspace*

An example of this file can be downloaded `here <https://github.com/Kia23/genereg/raw/master/DATA/sample_files/Genes_of_Interest.xlsx>`_. In this case, three different gene sets are analyzed, representing three genomic pathways relevant for the *Ovarian Serous Cystadenocarcinoma*. Each target gene belongs to a specific pathway according to the biological functions it is involved in.

|

-----------------------
Library Initialization
-----------------------
You have to execute the following function in order to initialize the library and to correctly set all the initial library settings:

``library_init()``

	The LIBRARY_INIT operation performs a preliminary initialization procedure and sets all the initial library settings.
	
	Example::

		import genereg as gr
		gr.Initialize.library_init()

This function builds the complete directory tree for storing files and results in the *library workspace*, it downloads and saves the complete lists of Transcription Factors (from both UniProt and ENCODE) and human genes (from HGNC) to be used for the **Gene - TFs Mapping** phase, and it creates some useful soon-to-be-used files, starting from the "*Genes_of_Interest.xlsx*" table.

|

**You are now ready to begin!**