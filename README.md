# Project Description

Project 1 provides analysis of a public microarray dataset. We reproduce results from comparisons of C3 and C4 tumor subtypes from Marisa et al. (2013). We normalized, filtered, and analyzed the data from 135 samples to discover relationships among them.

# Contributors

<<<<<<< HEAD
List contributor names and github user names, or email addresses if desired
=======
Names: Konrad Thorner, Aishwarya Deengar, Jia Liu, Morgan Rozman

Github: kthorner, AishwaryaD1, jialiu0103, morganroz

Email: kthorner@bu.edu, , jiliu@bu.edu, mrozman@bu.edu
 
# Repository Contents

### project1_jia.R

Normalization of microarray data. Compute standard quality control metrics on normalized data and visualize the distribution of samples using Principal Component Analysis.

### konrad_analyst.R

Filters the gene expression matrix, performs hierarchial clustering, constructs a heatmap, and find differential gene expression between clusters with t-tests.

### project1_aishwarya.Rmd
Scripts for mapping the probeset IDs to gene symbols, selecting the top 1000 up- and down-regulated genes and reporting the top 10 of these up- and down-regulated genes, fisher test to compute hypergeometric statistics and p-values comparing overlap for each gene set and each gene list, adjust the p-values for multiple hypotheses using the Benjamini-Hochberg (FDR) procedure and append this adjusted p-value column to the data frame and sorting the genesets according to their p values and reporting the table.
