## Clear environment
rm(list = ls())

## Load necessary libraries
library(edgeR)

## Set input and output file paths
input <- "input/StarvationMuscle_count_OB.tsv"                                    # Input file
output <- "output/result_StarvationMuscle_ANOVA_OB.txt"                           # Output file

## Load the data
count <- read.table(input, sep = "\t", header = TRUE, row.names = 1)
count <- as.matrix(count)

## Define groups (8 time points, 5 replicates each)
group_all <- factor(rep(1:8, each = 5))  # 8 time points, 5 replicates per time point

## Create design matrix for one-way ANOVA (using time points as groups)
design <- model.matrix(~ group_all)

## Create DGEList object with counts and group information
d <- DGEList(counts = count, group = group_all)

## Calculate normalization factors
d <- calcNormFactors(d)

## Estimate dispersion parameters
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

## Fit the GLM
fit <- glmFit(d, design)

## Perform likelihood ratio test (one-way ANOVA-like test)
lrt <- glmLRT(fit, coef = 2:length(unique(group_all)))  # Compare across all time points

## Get the top tags (differentially expressed genes)
table <- as.data.frame(topTags(lrt, n = nrow(count)))

## Write the results to a file
write.table(table, file = paste(output, sep = ""), col.names = TRUE, row.names = TRUE, sep = "\t")
