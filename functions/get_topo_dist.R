get_topo_dist <- function(tree_file_name, out_file_name = 'out_topo_dist.txt', n_clusters = 2){

  require(ape)
  require(doParallel)
  require(bigmemory)
  
  trees_tr <- read.tree(tree_file_name)
print('read trees')

  topo_mat <- big.matrix(ncol = length(trees_tr), nrow = length(trees_tr), dimnames = list(names(trees_tr), names(trees_tr)))

print('made matrix')

  library(doParallel)

  cat(paste('I will register', n_clusters, 'clusters\n'))
  cl <- makeCluster(n_clusters)

  registerDoParallel(cl)

for(tree_max in 2:length(trees_tr)){
cat('calculating distance ', tree_max-1, ' of ', length(trees_tr), '\n')

topo_mat[tree_max, 1:(tree_max - 1)] <- foreach(x = 1:(tree_max - 1), .combine = c) %dopar%  ape::dist.topo(trees_tr[[tree_max]], trees_tr[[x]])
}

  stopCluster(cl)
  
  write.big.matrix(topo_mat, filename = out_file_name, row.names = T, col.names = T)
  cat('Done. The PH85 tree distances are saved in', out_file_name, '\n')
}