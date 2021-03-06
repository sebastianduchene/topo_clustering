library(ape)

#fasta_files <- grep('fasta$', dir('../gdata_7/'), value = T)

nj_collapse <- function(dna_dat){

   
   boot_trees <- lapply(1:100, function(x) ape::nj(ape::dist.dna(dna_dat[, sample(1:ncol(dna_dat), ncol(dna_dat), replace = T)])))
   class(boot_trees) <- 'multiPhylo'


#   return(dna_nj_tr)
    return(ape::consensus(boot_trees, p = 0.51))
}

#f_1 <- read.dna(paste0('../gdata_7/', fasta_files[1]), format = 'fasta')

#test_1 <- nj_collapse(f_1)