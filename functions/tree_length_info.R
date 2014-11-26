library(ape)
library(cluster)

# example:
# tree_length_info('tree_file.txt', 'mds_dists.txt')
# The parameter are:
# optim_trees: trees with the same topology and branch lengths optimised. They should have names.
# mds_sbsdmin: MDS of the sbsdmin distances
# k_set: The number of clock clusters/ pacemakers
# tl_class: The number of tree length classes to use


tree_length_info <- function(optim_trees, mds_sbsdmin, k_set = 10, tl_class = 10, clustering_info_name = 'clus_info.txt', tree_info_name = 'tree_info.txt'){
  require(ape)
  require(cluster)		 

  all_trees <- read.tree(optim_trees)
  mds_mat <- read.table(mds_sbsdmin, as.is = T)

  n_p_groups <- ceiling(length(all_trees) / tl_class)

  tr_lens_sorted <- sort(sapply(all_trees, function(x) sum(x$edge.length)), dec = T)

  tree_sets <- seq(from = 1, to = length(tr_lens_sorted), by = n_p_groups)
  if(!(length(tr_lens_sorted) %in% tree_sets)) tree_sets <- c(tree_sets, length(all_trees))

  tree_set_list <- list()
  for(i in 1:(length(tree_sets) -1)){
    tree_set_list[[i]] <- (tree_sets[i] + 1):(tree_sets[i+1])
  }
  tree_set_list[[1]] <- c(1, tree_set_list[[1]])

  tree_set_mat <- matrix(NA, length(all_trees), 2)
  for(k in 1:length(tree_set_list)){
    tree_set_mat[tree_set_list[[k]], 1] <- tree_set_list[[k]]
    tree_set_mat[tree_set_list[[k]], 2] <- k
  }

  rownames(tree_set_mat) <- names(tr_lens_sorted)

  clust_data <- pam(mds_mat, k = k_set)
  clus_k <- as.matrix(clust_data$clustering)

  med_dist <- as.matrix(dist(clust_data$medoids))

  tree_set_mat_binded <- cbind(tree_set_mat, clus_k[match(rownames(tree_set_mat), rownames(clus_k)), ], tr_lens_sorted)
  colnames(tree_set_mat_binded) <- c('rank', 'tl_class', 'clock_cluster', 'tl')
  tree_set_mat_binded <- as.data.frame(tree_set_mat_binded)

  # get number of clusters within tree 'clubs'

  tree_quant_mat <- as.data.frame(matrix(NA, length(unique(tree_set_mat_binded$tl_class)), 5))
  tree_quant_mat[, 1] <- 1:nrow(tree_quant_mat)

  for(m in 1:nrow(tree_quant_mat)){
      tree_quant_mat[m, 2] <-   length(unique(tree_set_mat_binded$clock_cluster[tree_set_mat_binded$tl_class == m]))
      tree_quant_mat[m, 3] <- mean(tree_set_mat_binded$tl[tree_set_mat_binded$tl_class == m])
      class_clusts <- unique(tree_set_mat_binded$clock_cluster[tree_set_mat_binded$tl_class == m])
      tree_quant_mat[m, 4] <- mean(as.dist(med_dist[class_clusts, class_clusts], diag = F))
      tree_quant_mat[m, 5] <- paste(class_clusts, collapse = ',')
  }

  colnames(tree_quant_mat) <- c('tl_class', 'n_clusters', 'mean_tl', 'inter_clus_d', 'cluster_id')

  write.table(tree_quant_mat, file= clustering_info_name)
  write.table(tree_set_mat_binded, file = tree_info_name)

  return(list(tree_quant_mat, tree_set_mat_binded))
}









