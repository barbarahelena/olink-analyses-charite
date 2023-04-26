# Analysis Olink data

## Set-up mouse experiment
There were two conditions at 1 timepoint: 
- CI+PBS = Carotid Injury + PBS, this is the control, N = 5
- CI+S.pn. = Carotid Injury + streptococcus pneumonia, this is the infection group, N = 5

## Olink panel
The Olink mouse exploratory panel (ww.olink.com/mouse-exploratory) was performed using EDTA plasma from the mice (40 ul samples). This panel includes 92 proteins of which the list can be found in the data folder. All samples passed quality control (QC). However, some of the proteins had values under the limit of detection (LOD). The LOD is defined by Olink as 3 times the standard deviation over background. Olink indicated that they see that values under lod do generally not increase false-positive findings, and they recommend using all data including those falling below their detection limits. Therefore, no filtering was applied to the proteins in this analysis.

<img src="https://github.com/barbarahelena/olink-analyses-charite/blob/main/results/proteins_belowLOD.png" width = '300'>

## Explanation on data processing
In the data processing script (`230421_olink.R`), first we open the sample data frame and protein expression data frame. From the protein expression table, we extract the proteins and their metadata (such as olink ID and Uniprot ID). These data were saved into `proteins_metadata.csv` in the data folder.
In this script, we also select the proteins with less than 5 samples of interest (PBS or Strep) below the LOD.

For determining significant group differences, I used Welch's 2-sample t-tests for the 71 proteins. I applied FDR (aka Benjamini-Hochberg) multiple testing correction on the p-values to generate a q-value (= adjusted p-value). Thereafter, we calculate the group means, sd and n per protein.
The table with the results from the t-tests and the means per group was saved as `proteins_ttest_diff.csv`.

The expression data, so the protein expression values for all samples and the 71 proteins, were saved in `proteins_expression_data.csv`. For the heatmaps and boxplots I used scaled expression data, so the expression that was measured scaled to a mean of 0 and sd of 1, to make the proteins easier to visualize in for instance a heatmap. For the PCA plot i used unscaled data, as scaling is already performed in the PCA function itself.

## PCA
For the principal component analysis (PCA) in `230425_pca.R`, I used the `prcomp` function with scaling and centering of the data. The first two principal components are plotted. I made a version with labels to make it easier to assess which samples are more alike.

## Heatmaps
Heatmaps were drawn using the ComplexHeatmap R package (v2.12.1), see script `230425_heatmaps.R`. I made different versions: different selections of proteins (all proteins, or sig based on p-value or q-value (= FDR-adjusted)); showing a scale for significance and without. 

## Boxplots
Boxplots were drawn with ggplot2 while comparing the two groups with `stat_compare_means` from the ggpubr package (v0.6.0) using t-tests. See script `230425_boxplots.R`.

## Combination plots
The three plots above were combined in one plot in the `230425_combineplots.R` script.
