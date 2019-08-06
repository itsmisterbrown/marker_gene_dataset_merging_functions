#marker gene dataset merging functions

#calculate the total read count across the dataset
sumRuns <- function(x){
  y <- phyloseq::sample_sums(x)
  sum(y)
}

#calculate the geometric mean of an input vector of read counts
geomeans.r <- function(df, p=NULL, na.rm=FALSE){
  if (is.vector(df)){
    df <- matrix(df, nrow=1)
  }
  if (is.null(p)){
    exp(rowMeans(log(df), na.rm=na.rm))
  } else {
    exp((rowSums(log(df) %*% diag(p)))/sum(p))
  }
}

#merge several phyloseq-class objects. This function takes an unlimited number of input datasets and 
#merges them into a single phyloseq-class object. First, datasets are transformed into relative abundances 
#to preserve relative ratios. Datasets are then merged into a single object using the merge_phyloseq function
#within the phyloseq package. If count data are desired, relative abundances can be transformed into standardized
#count data based on either the geometric mean read count or the median read count of all datasets merged. 
#This is done by calculating the respective value, multiplying that statistic by the number of datasets merged, 
#and then dividing that value by the total number of samples in the merged dataset.

merge.dataset <- function(..., standardization="relative", return.all=FALSE){
  #create list of datasets to merge
  phy.list <- list(...)
  dataset.count <- ...length()
  
  #relative abundance transform to preserve ratios
  phy.list.rel <- lapply(phy.list, phyloseq::transform_sample_counts, fun = function(x) x/sum(x))
  #merge the relative abundance transformed datasets
  phy.obj.merged <- do.call(phyloseq::merge_phyloseq, phy.list.rel)
  #normalize to geometric mean read counts, if desired
  if (standardization=="geomeans"){
    runSums <- sapply(phy.list, FUN = sumRuns)
    run.geomeans <- geomeans.r(runSums)
    sample.geomeans <- ((run.geomeans * dataset.count)/phyloseq::nsamples(phy.obj.merged))
    phy.obj.merged <- phyloseq::transform_sample_counts(phy.obj.merged, fun = function(x) round(x * sample.geomeans))
    run.median <- NULL
    sample.median <- NULL
  } else if (standardization=="median"){
    runSums <- sapply(phy.list, FUN = sumRuns)
    run.median <- median(runSums)
    sample.median <- ((run.median * dataset.count)/phyloseq::nsamples(phy.obj.merged))
    phy.obj.merged <- phyloseq::transform_sample_counts(phy.obj.merged, fun = function(x) round(x * sample.median))
    run.geomeans <- NULL
    sample.geomeans <- NULL
  } else {
    run.median <- NULL
    run.geomeans <- NULL
    sample.geomeans <- NULL
    sample.median <- NULL
  }
  
  # Build return list
  l.return = list()
  if (return.all==FALSE){
    return(phy.obj.merged)
  } else {
    l.return[['merged.phyloseq']] <- phy.obj.merged
    l.return[['geometric.mean.dataset.read.count']] <- run.geomeans
    l.return[['median.dataset.read.count']] <- run.median
    l.return[['geometric.mean.sample.read.count']] <- sample.geomeans
    l.return[['median.sample.read.count']] <- sample.median
    
  }
  
  return(l.return)
}
