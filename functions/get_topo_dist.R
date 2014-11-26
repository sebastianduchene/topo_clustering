get_topo_dist <- function(tree_file_name, out_file_name = 'out_topo_dist.txt', n_clusters = 2){

  require(ape)
  require(doParallel)

  trees_tr <- read.tree(tree_file_name)

  topo_mat <- matrix(NA, length(trees_tr), length(trees_tr))

  colnames(topo_mat) <- names(trees_tr)
  rownames(topo_mat) <- names(trees_tr)

  for(i in 1:nrow(topo_mat)){
      topo_mat[i, ] <- paste(rownames(topo_mat)[i], colnames(topo_mat))
  }

  vector_names <- topo_mat[lower.tri(topo_mat)]

  get_topo_dist <- function(tree_list, tree_names){
	      tree_index_1 <- which(names(tree_list) == tree_names[1])
	      tree_index_2 <- which(names(tree_list) == tree_names[2])

	      if(length(tree_index_1) > 1 || length(tree_index_2) > 1){
	        stop('There are duplicate tree names in the list')
	      }

	      return(ape::dist.topo(tree_list[[tree_index_1]], tree_list[[tree_index_2]]))
  }

  library(doParallel)
  cat(paste('I will register', n_clusters, 'clusters\n'))
  cl <- makeCluster(n_clusters)
  registerDoParallel(cl)
  cat('Clusters registered. Open a terminal window and type top -u to see the number of R processes running\n')
  cat('I will calculate', length(vector_names), ' tree distances. Please be patient\n')
  
  topo_dists_line <- foreach(i = 1:length(vector_names), .combine = cbind) %dopar% get_topo_dist(trees_tr, strsplit(vector_names[i], ' ')[[1]])
  colnames(topo_dists_line) <- vector_names
  stopCluster(cl)

  trees_run <- which(lower.tri(topo_mat), arr.ind = T)

  for(i in 1:nrow(trees_run)){
    topo_mat[trees_run[i, 1], trees_run[i, 2]] <- topo_dists_line[i]     
  # print(topo_mat[trees_run[i, 1], 1])
  }

  write.table(topo_mat, file = out_file_name, row.names = T, col.names = T)
  cat('Done. The PH85 tree distances are saved in', out_file_name, '\n')
}