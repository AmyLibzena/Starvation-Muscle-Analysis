## Clear environment
rm(list = ls())

## Load necessary libraries
library(edgeR)

## Set input and output file paths
input <- "input/StarvationMuscle_count.tsv"  # Input file with both WT and ob/ob data
output_dir <- "output/"  # Directory for saving output files

## Load the data
count <- read.table(input, sep = "\t", header = TRUE, row.names = 1)
count <- as.matrix(count)

## Define groups (WT and ob/ob) and time points
genotypes <- factor(rep(c("WT", "OB"), each = 5 * 8))  # 5 replicates each for WT and ob/ob at 8 time points
time_points <- factor(rep(1:8, each = 5, times = 2))  # 8 time points, 5 replicates per time point for both genotypes

## Create a design matrix for pairwise comparison between WT and ob/ob at each time point
design <- model.matrix(~ 0 + genotypes:time_points)  # Interaction between genotypes and time points
colnames(design) <- make.names(levels(genotypes:time_points))  # Make column names syntactically valid

## Create DGEList object with counts and group information
d <- DGEList(counts = count, group = interaction(genotypes, time_points))

## Calculate normalization factors
d <- calcNormFactors(d)

## Estimate dispersion parameters
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

## Fit the GLM
fit <- glmFit(d, design)

## Actual time points for output file suffix
time_point_suffix <- c(0, 2, 4, 6, 8, 12, 16, 24)

## Perform pairwise comparison for each time point
for (tp in 1:8) {
  # Manually define the contrast names for each time point
  if (tp == 1) {
    contrast <- makeContrasts(WT.1 - OB.1, levels = design)
  } else if (tp == 2) {
    contrast <- makeContrasts(WT.2 - OB.2, levels = design)
  } else if (tp == 3) {
    contrast <- makeContrasts(WT.3 - OB.3, levels = design)
  } else if (tp == 4) {
    contrast <- makeContrasts(WT.4 - OB.4, levels = design)
  } else if (tp == 5) {
    contrast <- makeContrasts(WT.5 - OB.5, levels = design)
  } else if (tp == 6) {
    contrast <- makeContrasts(WT.6 - OB.6, levels = design)
  } else if (tp == 7) {
    contrast <- makeContrasts(WT.7 - OB.7, levels = design)
  } else if (tp == 8) {
    contrast <- makeContrasts(WT.8 - OB.8, levels = design)
  }
  
  # Perform likelihood ratio test for this time point
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Get the top tags (differentially expressed genes)
  table <- as.data.frame(topTags(lrt, n = nrow(count)))
  
  # Write the results to a file with a suffix indicating the actual time point (0, 2, 4, etc.)
  output_file <- paste(output_dir, "result_StarvationMuscle_ANOVA_", time_point_suffix[tp], "h.txt", sep = "")
  write.table(table, file = output_file, col.names = TRUE, row.names = TRUE, sep = "\t")
}
