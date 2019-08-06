# R functions for transforming, standardizing, and merging marker gene datasets to preserve relative ratios of taxa.

__merge.dataset__ takes an unlimited number of input datasets (phyloseq-class) and 
merges them into a single phyloseq-class object. First, datasets are transformed into relative abundances 
to preserve relative ratios. Datasets are then merged into a single object using the merge_phyloseq function
within the phyloseq package. If count data are desired, relative abundances can be transformed into standardized
count data based on either the geometric mean read count or the median read count of all datasets merged. 
This is done by calculating the respective value, multiplying that statistic by the number of datasets merged, 
and then dividing that value by the total number of samples in the merged dataset.