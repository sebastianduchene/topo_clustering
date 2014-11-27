

funs <- dir('functions', pattern = 'R$')
for(i in funs){
  source(paste0('functions/', i))
}

#opt_topo <- optim_clusters_topo('poly_topo_dists.txt', n_clusters = 5, kmax = 10)

# recursive function to convert the matrix to numeric. Limited by the stack size TOO DEEP RECURSION FOR LARGE OBJECTS. USE
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


#dist21 <- read_dist_matrix('collapsed_21_brs_tdist.txt')
mds52 <- optim_clusters_topo('collapsed_52_brs_tdist.txt', n_clusters = 10, kmax = 100, b_reps = 50, out_cluster_id = 'opt_topo_52.txt', out_clus_info = 'clusinfo_52.txt', plot_clustering = T)



