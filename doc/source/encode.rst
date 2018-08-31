ENCODE Data Extraction
============================================
The first data extraction phase consists in selecting from ENCODE all the Transcription Factors (TFs) having binding sites located in the promoter regions of genes of interest. Then, starting from these TFs, candidate regulatory genes for each gene of interest are identified:

.. image:: images/TFs.png


The query for extracting TFs of interest is implemented acconding to the `PyGMQL <https://pygmql.readthedocs.io/en/latest/index.html>`_ syntax and data are retrieved from public datasets available on the `GMQL <http://gmql.eu/gmql-rest/>`_
system.

|

-------------------------------------
Extraction of Transcription Factors
-------------------------------------

..  automodule:: TranscriptionFactors
    :members:

|

-------------------------------------
Identification of Regulatory Genes
-------------------------------------

..  automodule:: RegulatoryGenes
    :members:

The following image explains the extraction process for each target gene works and how its TFs and regulatory genes are identified:

.. image:: images/encode.png

