library(ape)
library(phangorn)
library(doParallel)

for(i in dir('functions', pattern = 'R$')){
      print(paste('load', i))
  source(paste0('functions/', i))
}

fasta_files <- paste0('remove_pandinas/', grep('.FASTA$', dir('remove_pandinas'), value = T))
get_nj_tree <- function(fasta_file, max_0_br = 10){
  f_temp <- ape::read.dna(fasta_file, format = 'fasta')
  di_temp <- phangorn::dist.hamming(phangorn::phyDat(f_temp))
  nj_temp <- ape::nj(di_temp)
  if(sum(nj_temp$edge.length == 0) < max_0_br){
    nj_poly <- nj_collapse(f_temp)
    return(list(fasta_file, nj_poly))
  }
}

###########################
# Estimate collapsed trees for a maximum of 20 branches with length of 0
cl <- makeCluster(15)
registerDoParallel(cl)
poly_trees_20 <- foreach(i = 1:length(fasta_files), .combine = c) %dopar% list(get_nj_tree(fasta_files[i], max_0_br = 20))
stopCluster(cl)
rem_20 <- sapply(1:length(poly_trees_20), function(x) !is.null(poly_trees_20[[x]]))
poly_trees_20 <- poly_trees_20[rem_20]

for(i in 1:length(poly_trees_20)){
  write.tree(poly_trees_20[[i]][[2]], file = 'collapsed_20_brs.trees', tree.names = gsub('^([a-z]|_|)+/', '', poly_trees_20[[i]][[1]]), append = TRUE)
}

###########################
# Estimate collapsed trees for a maximum of 10 branches with length of 0
cl <- makeCluster(15)
registerDoParallel(cl)
poly_trees_10 <- foreach(i = 1:length(fasta_files), .combine = c) %dopar% list(get_nj_tree(fasta_files[i], max_0_br = 10))
stopCluster(cl)
rem_10 <- sapply(1:length(poly_trees_10), function(x) !is.null(poly_trees_10[[x]]))
poly_trees_10 <- poly_trees_10[rem_10]

for(i in 1:length(poly_trees_10)){
  write.tree(poly_trees_10[[i]][[2]], file = 'collapsed_10_brs.trees', tree.names = gsub('^([a-z]|_|)+/', '', poly_trees_10[[i]][[1]]), append = TRUE)
}

############################
# Estimate collapsed trees for a maximum of 5 branches with length of 0
cl <- makeCluster(15)
registerDoParallel(cl)
poly_trees_5 <- foreach(i = 1:length(fasta_files), .combine = c) %dopar% list(get_nj_tree(fasta_files[i], max_0_br = 5))
stopCluster(cl)
rem_5 <- sapply(1:length(poly_trees_5), function(x) !is.null(poly_trees_5[[x]]))
poly_trees_5 <- poly_trees_5[rem_5]

for(i in 1:length(poly_trees_5)){
  write.tree(poly_trees_5[[i]][[2]], file = 'collapsed_5_brs.trees', tree.names = gsub('^([a-z]|_|)+/', '', poly_trees_5[[i]][[1]]), append = TRUE)
}




