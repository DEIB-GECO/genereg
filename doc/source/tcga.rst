TCGA Data Extraction
============================================
The second data extraction phase consists in extracting the methylation level (expressed as the **mean** of its *beta_values*) and the gene expression value of each target gene and for each data sample under analysis.
These types of data are currently available in GMQL for the TCGA dataset, which is the one that is used in this method.

The queries for extracting methylation and gene expression data of interest are implemented acconding to the `PyGMQL <https://pygmql.readthedocs.io/en/latest/index.html>`_ syntax and data are retrieved from public datasets available on the `GMQL <http://gmql.eu/gmql-rest/>`_
system.

Only common data samples having both methylation and expression values in TCGA are selected. Each data sample corresponds to a specific patient, which is here identified by a unique string, called *Sample Barcode*, with the following structure: **TCGA-xx-xxxx-xxx**.


-------------------------------------
Extraction of Methylation values
-------------------------------------

..  automodule:: genereg.Methylation
    :members:

|

Here it is an example of extraction: in this case the methylation sites of interest are the ones falling within slightly wider areas than the target genes promoters, going from 4000 bases upstream to 1000 bases downstream from the TSSs of the genes of interest:

.. image:: images/methyl.png

This interval can be useful, for example, for assessing the impact of a methylated promoter on the target gene, which participates in considerably reducing the gene expression.

Here it is a sample excerpt of the final methylation values table, containing TCGA data samples as rows and target genes as columns:

.. image:: images/methyltab.png

|

-------------------------------------
Extraction of Gene Expression values
-------------------------------------

..  automodule:: genereg.GeneExpression
    :members:

