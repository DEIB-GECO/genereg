Initialization
============================================

The library automatically creates directories for storing all the intermediate and final results in the current workspace of the user (i.e. the directory Python is executed in), in order to be able to use them during the execution, when needed.
For this reason, it's advised to execute Python direclty within the directory you want the results to be saved in (let's call it "*library workspace*") and it's highly recommended not to move any file until the whole process has been fully executed, that is all the functions of this library have been run and the whole analysis procedure has been perfomed.

------------------
Genes of Interest
------------------
The first thing you need to do is defining your set of genes of interest for the cencer type you want to analyze: these are the target genes whose regulation systems will be assesd during the analysis.

Here is how this set of genes has to be structured:
    
	* create an Excel file (.xlsx) with 4 columns ('GENE_SYMBOL', 'ENTREZ_GENE_ID', 'PATHWAY', 'SubClass') and one row for each target gene
	
		* GENE_SYMBOL is the official gene name

		* ENTREZ_GENE_ID is the numerical ID associated to the gene

		* PATHWAY is the genomic pathway the gene belongs to, usually defined according to the functions the gene is involved in

		* SubClass is the list of the functional classes the gene is assigned to within its own genomic pathway, defined according to the specific functions the gene participates in within its own patwhay

	* save this file as "*Genes_of_Interest.xlsx*" in the *library workspace*

An example of this file, containing target genes belonging to three relevant pathways for the *Ovarian Serous Cystadenocarcinoma*, can be downloaded `here <https://github.com/Kia23/genexpreg-cancer/raw/master/DATA/sample_files/Genes_of_Interest.xlsx>`_.

|

-----------------------
Library Initialization
-----------------------
You have to execute the following function in order to initialize the library and to correctly set all the initial library settings:

..  automodule:: Initialize
    :members:

This function builds the complete directory tree for storing files and results in the *library workspace*, it downloads and saves the complete lists of Transcription Factors (from both UniProt and ENCODE) and human genes (from HGNC) to be used for the **Gene - TFs Mapping** phase, and it creates some useful soon-to-be-used files, starting from the "*Genes_of_Interest.xlsx*" table.

|

**You are now ready to begin!**