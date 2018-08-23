# genexpreg_cancer

[![Documentation Status](https://readthedocs.org/projects/genexpreg-cancer/badge/?version=latest)](https://genexpreg-cancer.readthedocs.io/en/latest/?badge=latest)

This library provides a method for analyzing the regulation systems of target genes belonging to genomic pathways relevant for a specific TCGA cancer type. Adopting a linear regression approach, it builds a predictive model for the regulation of target genes expression, in order to identify the relevant features that explain the regulation of each gene of interest within patients affected by the tumor under analysis.
The regulation system of each target gene is analyzed singularly and independently, by assessing the effect of some specific factors on its expression:
* **DNA methylation** (more specifically, the mean promotorial methylation level of the target gene);
* **expression of target genes** belonging to the **same pathway** of the model gene;
* **expression of target genes** belonging to the **other pathways** with respect to the model gene;
* **expression of candidate regulatory genes** (i.e. those genes encoding for transcription factors having binding sites located in the promoter regions of genes of interest).

Thus, the matter is understanding the relationships between the activity of each target gene and the genes belonging either to the same pathway or to the other relevant pathways, and the relationships between all such target genes and their candidate regulatory genes: this may lead to identify potential common regulators along each pathway, or frequent regulators with a key role in the regulation systems of genes of interest, eventually predicting their potential oncogenic role. Whenever a correlation existsm an assessment of the potential influence that the gene methylation mey have on its expression is also made.
This approach does not analyze all the existing correlation among target genes and their regulatory features, but it identifies those associations that more contributes to the target genes expression regulation.
The results are focused by design to the **best-predicting sets of features**, leaving out potential regulators with important biological functions, but with an extremely low predictive power with respect to the expression of the target gene. 

The full documentation can be found at the following link: http://genexpreg-cancer.readthedocs.io
