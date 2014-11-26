library(ape)


#c5_trees <- read.tree('nj_trees.trees')
get_sbsd_dist <-  function(tree_file_name, out_file_name = 'out_sbsd_dist.txt',  n_clusters = 2){
  
  require(doParallel)
  require(ape)
  trees <- read.tree(tree_file_name)

  bsd.dist <- function(tree1 , tree2){
	       list.tr <- list()
	       list.tr[[1]] <- tree1
	       list.tr[[2]] <- tree2
        # In every case the tree that is rescaled is tree2, which is the shortest
	lens <- c(sum(tree1$edge.length), sum(tree2$edge.length))
	tree1 <- list.tr[lens==max(lens)][[1]]
	tree2 <- list.tr[lens==min(lens)][[1]]

        # This is the objective function for the optimisation implementation with optim()
        tree.dist.opt <- function(x){
            tree3 <- tree2
            tree3$edge.length <- tree2$edge.length*x
            return(ape::dist.topo(tree1, tree3, method="score"))
        }

        # Optimisation of the tree distance with a starting value of 0
        opt.dist <- optim(0, fn=tree.dist.opt, method = "Brent", lower = 0, upper = 50)

        min.bdi <- opt.dist$value
        scaling <- opt.dist$par

        # Scaling for tree2
        tree2.scaled <- tree2
        tree2.scaled$edge.length <- tree2$edge.length * scaling

        # The trees are scaled so that the mean branch length is 0.05
        # This is an arbitrary value, but is useful so that the total tree lengths are
        # in similar scales
	root.scaling <- 0.05 / mean(c(tree1$edge.length[tree1$edge.length > 0.00001] , tree2.scaled$edge.length[tree2.scaled$edge.length > 0.00001]))

	tree1.root.scaled <- tree1
	tree2.root.scaled <- tree2.scaled

	tree1.root.scaled$edge.length <- tree1$edge.length * root.scaling
	tree2.root.scaled$edge.length <- tree2.scaled$edge.length * root.scaling

	min.bdi.root.scaled <- ape::dist.topo(tree1.root.scaled, tree2.root.scaled, method="score")
        res.vect <- c(min.bdi.root.scaled, scaling, min.bdi)
        names(res.vect) <- c("min.bdi.scaled", "scaling.factor", "min.bdi")
        return(res.vect)
  }

  sbsd_mat <- matrix(NA, length(trees), length(trees))

  colnames(sbsd_mat) <- names(trees)
  rownames(sbsd_mat) <- names(trees)

  for(i in 1:nrow(sbsd_mat)){
      sbsd_mat[i, ] <- paste(rownames(sbsd_mat)[i], colnames(sbsd_mat))
  }

  vector_names <- sbsd_mat[lower.tri(sbsd_mat)]

  get_sbsd_dist <- function(tree_list, tree_names){
	      tree_index_1 <- which(names(tree_list) == tree_names[1])
	      tree_index_2 <- which(names(tree_list) == tree_names[2])

	      if(length(tree_index_1) > 1 || length(tree_index_2) > 1){
	        stop('There are duplicate tree names in the list')	      			      
	      }
	      
	      return(bsd.dist(tree_list[[tree_index_1]], tree_list[[tree_index_2]])[1])
  }

  cl <- makeCluster(n_clusters)
  registerDoParallel(cl)

  sbsd_dists <- foreach(i = 1:length(vector_names), .combine = cbind) %dopar% get_sbsd_dist(trees, strsplit(vector_names[i], ' ')[[1]])
  colnames(sbsd_dists) <- vector_names
  stopCluster(cl)

  trees_run <- which(lower.tri(sbsd_mat), arr.ind = T)

  for(i in 1:nrow(trees_run)){
      sbsd_mat[trees_run[i, 1], trees_run[i, 2]] <- sbsd_dists[i]
  }
  write.table(sbsd_mat, file = out_file_name, row.names = T, col.names = T)
  return(sbsd_mat)

}

#sbsd_1 <- sbsdmin_mat(c1_trees, 5)
#write.table(sbsd_1, file = 'sbsd_c1.txt', row.names = T, col.names = T)
#sbsd_5 <- get_sbsd_dist('twenty.trees', n_clusters = 4)
#write.table(sbsd_5, file = 'sbsd.txt', row.names = T, col.names = T)
