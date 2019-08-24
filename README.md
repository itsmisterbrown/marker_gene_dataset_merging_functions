# R functions for transforming, standardizing, and merging marker gene datasets to preserve relative ratios of taxa.

__merge.dataset__ takes an unlimited number of input datasets (phyloseq-class) and 
merges them into a single phyloseq-class object. First, datasets are transformed into relative abundances 
to preserve relative ratios. Datasets are then merged into a single object using the merge_phyloseq function
within the phyloseq package. If count data are desired, relative abundances can be transformed into standardized
count data based on either the geometric mean read count or the median read count of all datasets merged. 

for standardizing to the geometric mean read count per sample, use the parameter __standardization="geomeans"__  
</br>
for standardizing to the median read count per sample, use the parameter __standardization="median"__  
</br>
for a relative abundance-transformed dataset (default), use the parameter __standardization="relative"__ or leave blank
