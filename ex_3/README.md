# Exercise 3

## Differential Expression Analysis

We computed the differential expression using the results() function, with the level of 0.05. The threshold for log-fold change was set to 1 (so the change is at least twice as high/low).

Then we also sortered the results by the adjusted p-value to see the most prominent different expressions.

## MA-plot

We plotted two versions of the MA-plot. The log fold usually includes a lot of noise from genes that are not supported by high count of reads. The first plot does not account for this noise, while the second plot does. We used alpha=0.05.

MA plot with noise:

![With noise](MA_plot_noise.png)

MA plot without noise:

![Without noise](MA_plot.png)

### Plotting

In the next task, we plotted the normalized counts for the GJB2 gene (a gene with the 9th lowest adjusted p-value - done in the first task by sorting the adjusted p-values). The function `plotCounts()` automatically normalizes counts by the estimated size factors, so no prior normalization is necessary. Here is the resulting plot:

![Normalized counts of expression of the GJB2 gene](GJB2_normalized.png)


