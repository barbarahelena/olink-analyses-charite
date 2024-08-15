# Analysis Olink data

## Set-up mouse experiment
There were two conditions at 1 timepoint: 
- Control + CI = Carotid Injury + PBS, this is the control, N = 5
- Pneumonia + CI = Carotid Injury + streptococcus pneumonia, this is the infection group, N = 5

## Olink panel
The Olink mouse exploratory panel (ww.olink.com/mouse-exploratory) was performed using EDTA plasma from the mice (40 ul samples). This panel includes 92 proteins of which the list can be found in the data folder. All samples passed quality control (QC).

## Explanation on data processing
In the data processing script (`230421_olink.R`), first we open the sample data frame and protein expression data frame. From the protein expression table, we extract the proteins and their metadata (such as olink ID and Uniprot ID).

For determining significant group differences, I used 2-sided t-tests for the 92 proteins. I applied FDR (Benjamini-Hochberg) multiple testing correction on the p-values to generate a q-value (= adjusted p-value). Thereafter, we calculate the group means, sd and n per protein.
The table with the results from the t-tests and the means per group was saved as `proteins_ttest_diff.csv` in the results folder.

The expression data, so the protein expression values for all samples and the 92 proteins, were saved in `proteins_expression_data.csv`. For the heatmaps and boxplots I used scaled expression data, so the expression that was measured scaled to a mean of 0 and sd of 1, to make the proteins easier to visualize in for instance a heatmap. For the PCA plot i used unscaled data, as scaling is already performed in the PCA function itself.

## PCA
For the principal component analysis (PCA) in `230425_pca.R`, I used the `prcomp` function with scaling and centering of the data. The first two principal components are plotted. I made a version with labels to make it easier to assess which samples are more alike.

## Heatmaps
Heatmaps were drawn using the ComplexHeatmap R package (v2.12.1), see script `230425_heatmaps.R`. I made different versions: different selections of proteins (all proteins, or sig based on p-value or q-value (= FDR-adjusted)); showing a scale for significance and without. 

## Boxplots
Boxplots were drawn with ggplot2 while comparing the two groups with `stat_compare_means` from the ggpubr package (v0.6.0) using t-tests. See script `230425_boxplots.R`.

## Combination plots
The three plots above were combined in one plot in the `230425_combineplots.R` script.

## R packages
To improve reproducibility, we created a lockfile in August 2024 and reran all the analyses to ensure that everything still works. The renv lockfile, activate.R and Rprofile can be found in this repository. Please refer to the renv documentation for instructions on how to use these to recreate the R environment with all packages needed.


