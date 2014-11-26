funs <- dir('functions', pattern = 'R$')
for(i in funs){
  source(paste0('functions/', i))
}

#opt_topo <- optim_clusters_topo('poly_topo_dists.txt', n_clusters = 5, kmax = 10)

# recursive function to convert the matrix to numeric. Limited by the stack size
mat_to_num <- function(mat){
  if(ncol(mat) == 2){
    return(cbind(as.numeric(mat[, 1]), as.numeric(mat[, 2])))
    }else if(ncol(mat) > 2){
      return(cbind(as.numeric(mat[, 1]), mat_to_num(mat[, -1])))
      }else{
        return(matrix(as.numeric(mat), nrow(mat), ncol(mat)))
      }
}


read_dist_matrix <- function(f_name){
  dist_mat <- as.matrix(read.table(f_name, head = T, as.is = T))
  dist_mat[upper.tri(dist_mat)] <- NA
  dist_mat <- matrix(as.numeric(dist_mat), nrow(dist_mat), ncol(dist_mat), dimnames = list(rownames(dist_mat), colnames(dist_mat)))
  return(dist_mat)
}

#dist20 <- read_dist_matrix('collapsed_20_brs_tdist.txt')
#mds20 <- optim_clusters_topo('collapsed_20_brs_tdist.txt', n_clusters = 10, kmax = 100, b_reps = 50, out_cluster_id = 'opt_topo_20.txt', out_clus_info = 'clusinfo_20.txt', plot_clustering = T)

#dist10 <- read_dist_matrix('collapsed_10_brs_tdist.txt')
#mds10 <- optim_clusters_topo('collapsed_10_brs_tdist.txt', n_clusters = 10, kmax = 100, b_reps = 50, out_cluster_id = 'opt_topo_10.txt', out_clus_info = 'clusinfo_10.txt', plot_clustering = T)

dist5 <- read_dist_matrix('collapsed_5_brs_tdist.txt')
mds10 <- optim_clusters_topo('collapsed_5_brs_tdist.txt', n_clusters = 10, kmax = 100, b_reps = 50, out_cluster_id = 'opt_topo_5.txt', out_clus_info = 'clusinfo_5.txt', plot_clustering = T)
