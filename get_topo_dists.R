library(doParallel)

for(i in dir('functions', pattern = 'R$')){
      cat('loading', i, '\n')
      source(paste0('functions/', i))
}


cat('starting get topo\n')

get_topo_dist(tree_file_name = 'collapsed_20_brs.trees', out_file_name = 'collapsed_20_brs_tdist.txt', n_clusters = 10)

get_topo_dist(tree_file_name = 'collapsed_10_brs.trees', out_file_name = 'collapsed_10_brs_tdist.txt', n_clusters = 10)

get_topo_dist(tree_file_name = 'collapsed_5_brs.trees', out_file_name = 'collapsed_5_brs_tdist.txt', n_clusters = 10)
