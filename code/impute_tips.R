library(tidyverse)
library(phytools)
library(parallel)
source("code/pe_samp_funs.R")

# load full empirical tree
fern_tre <- read.tree("output/nitta_node-named.tre") 

# load taxonomy table
tax_tbl <- read.csv("output/fern_taxonomy_ncbi.csv") 

# load pruned trees
rand_fls <- list.files("output/pruned_trees", "rand_pruned",     all.files = T, full.names = T)
late_fls <- list.files("output/pruned_trees", "lateral_pruned",  all.files = T, full.names = T)
marg_fls <- list.files("output/pruned_trees", "marginal_pruned", all.files = T, full.names = T)

dir.create("output/imputed_trees/")

run_rand_rand <- T ; run_late_rand <- T ; run_marg_rand <- T # pruned and randomly imputed trees
run_rand_tax <- T  ; run_late_tax <- T  ; run_marg_tax <- T  # pruned and taxonomically imputed trees

n_cores <- 3

# add tips location at random ####
# random pruned
rand_frc <- readr::parse_number(rand_fls)
if(run_rand_rand){
  for(j in 1:length(rand_fls)){
    targ_tr_ls <- readRDS(rand_fls[[j]])
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) 
    length(targ_tr_ls) == length(tips2add_ls)
    
    impt_tr_rmd <- list()
    for(i in 1:length(targ_tr_ls)){
      impt_tr_rmd[[i]] <- add.random(tree = targ_tr_ls[[i]], tips = tips2add_ls[[i]]) 
    }
    impt_tr_rmd %>% saveRDS(paste0("output/imputed_trees/rand_pr_rand_imp_", rand_frc[j], ".rds"))
  }
} 

# lateral pruned
late_frc <- readr::parse_number(late_fls)
if(run_late_rand) {
  for(j in 1:length(late_fls)){
    targ_tr_ls <- readRDS(late_fls[[j]])
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) 
    length(targ_tr_ls) == length(tips2add_ls)
    
    impt_tr_rmd <- list()
    for(i in 1:length(targ_tr_ls)){
      impt_tr_rmd[[i]] <- add.random(tree = targ_tr_ls[[i]], tips = tips2add_ls[[i]]) 
    }
    impt_tr_rmd %>% saveRDS(paste0("output/imputed_trees/late_pr_rand_imp_", late_frc[j], ".rds"))
  }
}

# marginal pruned
marg_frc <- readr::parse_number(marg_fls)
if(run_marg_rand){
  for(j in 1:length(marg_fls)){
    targ_tr_ls <- readRDS(marg_fls[[j]])
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) 
    length(targ_tr_ls) == length(tips2add_ls)
    
    impt_tr_rmd <- list()
    for(i in 1:length(targ_tr_ls)){
      impt_tr_rmd[[i]] <- add.random(tree = targ_tr_ls[[i]], tips = tips2add_ls[[i]]) 
    }
    impt_tr_rmd %>% saveRDS(paste0("output/imputed_trees/marg_pr_rand_imp_", marg_frc[j], ".rds"))
  }
}

# add tips to based on taxonomy ####
# random pruned
rand_frc <- readr::parse_number(rand_fls)
if(run_rand_tax){
  for(i in 1:length(rand_fls)){ 
    trees_ls <- readRDS(rand_fls[[i]]) # load set of trees
    trees_ls <- mclapply(trees_ls, name_nodes, m_tax = tax_tbl, mc.cores = n_cores) # name nodes 
    
    targ_tr_ls <- trees_ls # set trees to impute
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) # get missing tips per tree
    
    for(k in 1:length(targ_tr_ls)){
      gens2add <- strsplit(tips2add_ls[[k]], '_', fixed = F) %>% # get the genera to add in the target tree
        sapply("[", 1) %>% unique()
      gens_tree <- strsplit(targ_tr_ls[[k]]$tip.label, '_', fixed = F) %>% # get the genera in the target tree
        sapply("[", 1) %>% unique() 
        
      gen_diff <- setdiff(gens2add, gens_tree) # get the missing genera in the target tree
      
      tps_gen_impt <- c()
      if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be inputed first
        for(j in 1:length(gen_diff)){
          nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(family) %>% # get node to impute genus at family level
            unique() %>% grepl(pattern = ., x = targ_tr_ls[[k]]$node.label) %>% which
          
          ifelse(test = !is_empty(nd_where), yes = print(paste0("Imputed at family level")), 
                 no = { nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(order) %>% # get node to impute genus at order
                   unique() %>% grepl(pattern = ., x = targ_tr_ls[[k]]$node.label) %>% 
                   which ; ifelse(!is_empty(nd_where), yes = print(paste0("Imputed at order level")),
                                  no = {nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(class) %>%  # get node to impute genus at class
                                    unique() %>% grepl(pattern = ., x = targ_tr_ls[[k]]$node.label) %>% 
                                    which ; ifelse(test = !is_empty(nd_where), 
                                                   yes = {nd_where <- 1 ; print(paste0("Imputed at root node"))},
                                                   no = print(paste0("Imputed at class level")))
                                  }
                   )
                 }
          )
          
          tps_2add <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(tip_names)
          targ_tr_ls[[k]] <- bind.tip(tree = targ_tr_ls[[k]], tip.label = tps_2add[1], where = c(Ntip(trees_ls[[k]])+nd_where)) # adds the first tip/genus
          tps_gen_impt[j] <- tps_2add[1]
        }
        impt_tips <- setdiff(tips2add_ls[[k]], tps_gen_impt)
        targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = impt_tips, where = "random")

      } else { # otherwise impute species to genus already in the tree
        targ_tr_ls[[k]] <- add.species.to.genus(targ_tr_ls[[k]], tips2add_ls[[k]], where = "random")
      }
      ifelse(targ_tr_ls[[k]]$tip.label %in% fern_tre$tip.label %>% sum == length(fern_tre$tip.label), 
             print(paste0("All tips imputed")), print(paste0("WARNING ! unequal number of tips")))
      print(paste0("-> tree: ", k))
    }
    print(paste0("## file: ", i))
    targ_tr_ls %>% saveRDS(paste0("output/imputed_trees/rand_pr_tax_imp_", rand_frc[i], ".rds"))
  }
} else {
  #targ_tr_ls <- readRDS(paste0("output/imputed_trees/rand_pr_tax_imp_", rand_frc[i], ".rds"))
}

sapply(targ_tr_ls, Ntip)
sapply(targ_tr_ls, is.rooted) %>% table()
sapply(targ_tr_ls, is.ultrametric) %>% table()

# lateral pruned
late_frc <- readr::parse_number(late_fls)
if(run_late_tax){
  for(i in 1:length(late_fls)){ 
    trees_ls <- readRDS(late_fls[[i]]) # load set of trees
    trees_ls <- mclapply(trees_ls, name_nodes, m_tax = tax_tbl, mc.cores = n_cores) # name nodes 
    
    targ_tr_ls <- trees_ls # set trees to impute
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) # get missing tips per tree
    
    for(k in 1:length(targ_tr_ls)){
      gens2add <- strsplit(tips2add_ls[[k]], '_', fixed = F) %>% # get the genera to add in the target tree
        sapply("[", 1) %>% unique()
      gens_tree <- strsplit(targ_tr_ls[[k]]$tip.label, '_', fixed = F) %>% # get the genera in the target tree
        sapply("[", 1) %>% unique() 
      
      gen_diff <- setdiff(gens2add, gens_tree) # get the missing genera in the target tree
      
      tps_gen_impt <- c()
      if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be inputed first
        for(j in 1:length(gen_diff)){
          nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(family) %>% # get node to impute genus at family level
            unique() %>% grepl(pattern = ., x = targ_tr_ls[[k]]$node.label) %>% which
          
          ifelse(test = !is_empty(nd_where), yes = print(paste0("Imputed at family level")), 
                 no = { nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(order) %>% # get node to impute genus at order
                   unique() %>% grepl(pattern = ., x = targ_tr_ls[[k]]$node.label) %>% 
                   which ; ifelse(!is_empty(nd_where), yes = print(paste0("Imputed at order level")),
                                  no = {nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(class) %>%  # get node to impute genus at class
                                    unique() %>% grepl(pattern = ., x = targ_tr_ls[[k]]$node.label) %>% 
                                    which ; ifelse(test = !is_empty(nd_where), 
                                                   yes = {nd_where <- 1 ; print(paste0("Imputed at root node"))},
                                                   no = print(paste0("Imputed at class level")))
                                  }
                   )
                 }
          )
          
          tps_2add <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(tip_names)
          targ_tr_ls[[k]] <- bind.tip(tree = targ_tr_ls[[k]], tip.label = tps_2add[1], where = c(Ntip(trees_ls[[k]])+nd_where)) # adds the first tip/genus
          tps_gen_impt[j] <- tps_2add[1]
        }
        impt_tips <- setdiff(tips2add_ls[[k]], tps_gen_impt)
        targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = impt_tips, where = "random")
        
      } else { # otherwise impute species to genus already in the tree
        targ_tr_ls[[k]] <- add.species.to.genus(targ_tr_ls[[k]], tips2add_ls[[k]], where = "random")
      }
      ifelse(targ_tr_ls[[k]]$tip.label %in% fern_tre$tip.label %>% sum == length(fern_tre$tip.label), 
             print(paste0("All tips imputed")), print(paste0("WARNING ! unequal number of tips")))
      print(paste0("-> tree: ", k))
    }
    print(paste0("## file: ", i))
    targ_tr_ls %>% saveRDS(paste0("output/imputed_trees/late_pr_tax_imp_", late_frc[i], ".rds"))
  }
} else {
  #targ_tr_ls <- readRDS(paste0("output/imputed_trees/rand_pr_tax_imp_", rand_frc[i], ".rds"))
}

sapply(targ_tr_ls, Ntip)
sapply(targ_tr_ls, is.rooted) #%>% table()
sapply(targ_tr_ls, is.ultrametric)# %>% table()

# marginal pruned
marg_frc <- readr::parse_number(marg_fls)


sapply(targ_tr_ls, Ntip)
sapply(targ_tr_ls, is.rooted) %>% table()
sapply(targ_tr_ls, is.ultrametric)# %>% table()

# vis comparative trees
library(ggtree) ; library(ggpubr)

t1 <- fern_tre %>% 
  ggtree(layout = "roundrect", ladderize = F) + 
  geom_nodepoint(size = 1.5, alpha = .55)
t2 <- targ_tr_ls[[50]] %>% 
  ggtree(layout = "roundrect", ladderize = F) +
  geom_nodepoint(size = 1.5, alpha = .55)

ggarrange(t1, t2, ncol = 2, labels = c("A", "B"))




