# This file implements data loading functions

# NYT dataset code is based on code from
# https://github.com/RoheLab/vsp-paper
# -----------------------------------------------------------------


source("functions/util.R")



# Generates a simulated sparse Gaussian dataset
get_sparse_normal = function(p0, n=1000, k=25, all_sparse=F) {
  # p0: probability of zero
  # n: number of rows
  # k: number of columns
  # all_sparse: T if sparsity applies to all columns; T if all but last column
  
  A = matrix(0, n, k)
  for (j in 1:k)
  {
    if (j==k & !all_sparse)
    {
      # No sparsity for last column
      A[,j] = rnorm(n)
    } else {
      A[,j] = rbinom(n,1,1-p0) * rnorm(n)
    }
  }
  
  # Skew positive and order by column variance
  A = skew_positive(A)
  A_var = apply(A, 2, var)
  A = A[,sort(A_var, decreasing=T, index.return=T)$ix]
  
  return(A)
}



# Downloads/loads NYT articles dataset
get_NYT = function(file, from_saved=T, save_data=F) {
  # file: file name to load/save data
  # from_saved: T if load from saved data; F if download
  # save_data: T if save data; F otherwise
  
  if (from_saved)
  {
    load(file=file)
  } else {
    # Download from https://archive.ics.uci.edu/ml/datasets/Bag+of+Words
    el = read_delim("https://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/docword.nytimes.txt.gz",
                    skip=3, delim=" ", col_names=F)
    el = as.tbl(el) %>% as.matrix
    A = spMatrix(nrow=max(el[,1]),ncol=max(el[,2]),i=el[,1],j=el[,2],x=el[,3]) %>% as("dgCMatrix")
    
    # Add words as column names
    words = read_delim("https://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/vocab.nytimes.txt",
                       delim="\n", col_names=F)$X1
    colnames(A) = words
    
    # Save data
    if (save_data) {save(A,file=file)}
  }
  return(A)
}