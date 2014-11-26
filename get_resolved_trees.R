library(ape)
library(phangorn)
library(doParallel)

for(i in dir('functions', pattern = 'R$')){
      print(paste('load', i))
  source(paste0('functions/', i))
}

fasta_files <- paste0('~/Desktop/phyt_infestans/remove_pandinas/', grep('.FASTA$', dir('~/Desktop/phyt_infestans/remove_pandinas'), value = T))
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
# Estimate collapsed trees for a maximum of 21 branches with length of 0 (20%)
cl <- makeCluster(15)
registerDoParallel(cl)
poly_trees_21 <- foreach(i = 1:length(fasta_files), .combine = c) %dopar% list(get_nj_tree(fasta_files[i], max_0_br = 21))
stopCluster(cl)
rem_21 <- sapply(1:length(poly_trees_21), function(x) !is.null(poly_trees_21[[x]]))
poly_trees_21 <- poly_trees_21[rem_21]

for(i in 1:length(poly_trees_21)){
  write.tree(poly_trees_21[[i]][[2]], file = 'collapsed_21_brs.trees', tree.names = gsub('^([a-z]|_|)+/', '', poly_trees_21[[i]][[1]]), append = TRUE)
}

###########################
# Estimate collapsed trees for a maximum of 52 branches with length of 0 (50%)
cl <- makeCluster(15)
registerDoParallel(cl)
poly_trees_52 <- foreach(i = 1:length(fasta_files), .combine = c) %dopar% list(get_nj_tree(fasta_files[i], max_0_br = 52))
stopCluster(cl)
rem_52 <- sapply(1:length(poly_trees_52), function(x) !is.null(poly_trees_52[[x]]))
poly_trees_52 <- poly_trees_52[rem_52]

for(i in 1:length(poly_trees_52)){
  write.tree(poly_trees_52[[i]][[2]], file = 'collapsed_52_brs.trees', tree.names = gsub('^([a-z]|_|)+/', '', poly_trees_52[[i]][[1]]), append = TRUE)
}
