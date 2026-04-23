library(tidyverse)
library(phytools) 
library(parallel)

source("code/pe_samp_funs.R")

# load full empirical tree
fern_tre <- read.tree("output/nitta_node-named.tre") 

# load taxonomy table
tax_tbl <- read.csv("output/fern_taxonomy_ncbi.csv") %>% 
  mutate(tip_names = sub(pattern = "(\\_.*?)\\_", "\\1-", tip_names))

# pruning conditional controls
run_rand_rand <- T ; run_late_rand <- T ; run_marg_rand <- T # pruned and randomly imputed trees
run_rand_tax <- T  ; run_late_tax <- T  ; run_marg_tax <- T  # pruned and taxonomically imputed trees
meth2impt <- c("root", "random") 
n_cores <- 3 # 8

names2load <- c("fern", "young", "old", "cnt", "exp")

# imputation ####
for(m in 1:length(names2load)){
  to_load <- names2load[m]
  out_fold <- to_load

# load pruned trees
rand_fls <- list.files(paste0("output/pruned_trees/", to_load), "rand_pruned",     all.files = T, full.names = T)
late_fls <- list.files(paste0("output/pruned_trees/", to_load), "lateral_pruned",  all.files = T, full.names = T)
marg_fls <- list.files(paste0("output/pruned_trees/", to_load), "marginal_pruned", all.files = T, full.names = T)

rand_frc <- readr::parse_number(rand_fls)
late_frc <- readr::parse_number(late_fls)
marg_frc <- readr::parse_number(marg_fls)

dir.create(paste0("output/imputed_trees/", to_load), recursive = T)

## add tips at random ####
# random pruned
if(run_rand_rand){
  for(j in 1:length(rand_fls)){
    targ_tr_ls <- readRDS(rand_fls[[j]])
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) 
    length(targ_tr_ls) == length(tips2add_ls)
    
    impt_tr_rmd <- list()
    for(i in 1:length(targ_tr_ls)){
      impt_tr_rmd[[i]] <- add.random(tree = targ_tr_ls[[i]], tips = tips2add_ls[[i]]) 
    }
    impt_tr_rmd %>% 
      saveRDS(paste0("output/imputed_trees/", out_fold,"/rand_pr_rand_imp_", rand_frc[j], ".rds"))
  }
} 

# lateral pruned
if(run_late_rand) {
  for(j in 1:length(late_fls)){
    targ_tr_ls <- readRDS(late_fls[[j]])
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) 
    length(targ_tr_ls) == length(tips2add_ls)
    
    impt_tr_rmd <- list()
    for(i in 1:length(targ_tr_ls)){
      impt_tr_rmd[[i]] <- add.random(tree = targ_tr_ls[[i]], tips = tips2add_ls[[i]]) 
    }
    impt_tr_rmd %>% 
      saveRDS(paste0("output/imputed_trees/", out_fold, "/late_pr_rand_imp_", late_frc[j], ".rds"))
  }
}

# marginal pruned
if(run_marg_rand){
  for(j in 1:length(marg_fls)){
    targ_tr_ls <- readRDS(marg_fls[[j]])
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) 
    length(targ_tr_ls) == length(tips2add_ls)
    
    impt_tr_rmd <- list()
    for(i in 1:length(targ_tr_ls)){
      impt_tr_rmd[[i]] <- add.random(tree = targ_tr_ls[[i]], tips = tips2add_ls[[i]]) 
    }
    impt_tr_rmd %>% 
      saveRDS(paste0("output/imputed_trees/", out_fold,"/marg_pr_rand_imp_", marg_frc[j], ".rds"))
  }
}

## add tips to tree based on taxonomy ####
for(f in 1:length(meth2impt)){
  where2impt <- meth2impt[f]

  # random pruned
if(run_rand_tax){ 
  for(i in 1:length(rand_fls)){ 
    trees_ls <- readRDS(rand_fls[[i]]) # load set of trees
     for(h in 1:length(trees_ls)){trees_ls[[h]]$node.label <- c()}
     trees_ls <- mclapply(trees_ls, name_nodes, m_tax = tax_tbl, mc.cores = n_cores) # name nodes 
    
    targ_tr_ls <- trees_ls # set trees to impute
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) # get missing tips per tree
    
    for(k in 1:length(targ_tr_ls)){
      gens2add <- strsplit(tips2add_ls[[k]], '_', fixed = F) %>% # get the genera from the tips to be added in the target tree
        sapply("[", 1) %>% unique() %>% sort()
      gens_tree <- strsplit(targ_tr_ls[[k]]$tip.label, '_', fixed = F) %>% # get the genera in the target tree
        sapply("[", 1) %>% unique() %>% sort
        
      gen_diff <- setdiff(gens2add, gens_tree) %>% sort # get the missing genera in the target tree that need to be imputed first
      
      tps_gen_impt <- c() ; impt_tps_vec <- c()
      if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be imputed first
        for(j in 1:length(gen_diff)){
          nd_where <- node2where(gen = gen_diff[j], tree = targ_tr_ls[[k]]) # get the nodes to impute genus
          tps_2add <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(tip_names) # get tips to impute per genus
          
          # impute the first tip per genus to anchor the clade
          if(length(tps_2add) > 1) {  # in case of more than one tip
            targ_tr_ls[[k]] <- TreeTools::AddTip(tree = targ_tr_ls[[k]], 
                                                 where = c(Ntip(targ_tr_ls[[k]]) + nd_where), label = tps_2add[1]) %>% 
              castor::extend_tree_to_height() %>% .$tree
            for(l in 2:length(tps_2add)){ # then impute the other tips in the same genus
              tem_tr <- add.species.to.genus(tree = targ_tr_ls[[k]], 
                                             species = tps_2add[l], where = where2impt)
              if(is.ultrametric(tem_tr)) {targ_tr_ls[[k]] <- tem_tr} else {targ_tr_ls[[k]] <- force.ultrametric(tem_tr, method = "nnls", message = F)}
            }
          } else { # in case of one tip per genus
            targ_tr_ls[[k]] <- TreeTools::AddTip(tree = targ_tr_ls[[k]], 
                                                 where = c(Ntip(targ_tr_ls[[k]]) + nd_where), label = tps_2add) %>% 
              castor::extend_tree_to_height() %>% .$tree
          }
          impt_tps_vec[[k]] <- c(impt_tps_vec[[k]], tps_2add)
        }
        
        tips2add_ls[[k]] <- setdiff(tips2add_ls[[k]], impt_tps_vec[[k]]) %>% sort
        for(l in 1:length(tips2add_ls[[k]])){
          tem_tr <- add.species.to.genus(tree = targ_tr_ls[[k]], 
                                         species = tips2add_ls[[k]][l], where = where2impt)
          if(is.ultrametric(tem_tr)) {targ_tr_ls[[k]] <- tem_tr} else {targ_tr_ls[[k]] <- force.ultrametric(tem_tr, method = "nnls", message = F)}
        }

      } else { # otherwise impute species to genus already in the tree
        for(l in 1:length(tips2add_ls[[k]])){
          tem_tr <- add.species.to.genus(tree = targ_tr_ls[[k]], 
                                         species = tips2add_ls[[k]][l], where = where2impt)
          if(is.ultrametric(tem_tr)) {targ_tr_ls[[k]] <- tem_tr} else {targ_tr_ls[[k]] <- force.ultrametric(tem_tr, method = "nnls", message = F)}
          
          }
      }
      ifelse(targ_tr_ls[[k]]$tip.label %in% fern_tre$tip.label %>% sum == length(fern_tre$tip.label), 
             print(paste0("*** All tips imputed ***")), print(paste0("WARNING ! unequal number of tips")))
      print(paste0("-> tree: ", k))
    }
    print(paste0("## file: ", i))
    targ_tr_ls %>% 
      saveRDS(paste0("output/imputed_trees/", out_fold,"/rand_pr_tax_imp_", where2impt, "_", rand_frc[i], ".rds"))
  }
} else {
  for(i in 1:length(rand_fls)){
    targ_tr_ls <- 
      readRDS(paste0("output/imputed_trees/", out_fold,"/rand_pr_tax_imp_", where2impt, "_",  rand_frc[i], ".rds"))
    sapply(targ_tr_ls, Ntip) %>% table() %>% print
  }
}

sapply(targ_tr_ls, Ntip)
sapply(targ_tr_ls, is.rooted) #%>% table()
sapply(targ_tr_ls, is.ultrametric) #%>% table()

# lateral pruned
if(run_late_tax){ 
  for(i in 1:length(late_fls)){ 
    trees_ls <- readRDS(late_fls[[i]]) # load set of trees
    for(h in 1:length(trees_ls)){trees_ls[[h]]$node.label <- c()}
    trees_ls <- mclapply(trees_ls, name_nodes, m_tax = tax_tbl, mc.cores = n_cores) # name nodes 
    
    targ_tr_ls <- trees_ls # set trees to impute
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) # get missing tips per tree
    
    for(k in 1:length(targ_tr_ls)){
      gens2add <- strsplit(tips2add_ls[[k]], '_', fixed = F) %>% # get the genera from the tips to be added in the target tree
        sapply("[", 1) %>% unique() %>% sort()
      gens_tree <- strsplit(targ_tr_ls[[k]]$tip.label, '_', fixed = F) %>% # get the genera in the target tree
        sapply("[", 1) %>% unique() %>% sort
      
      gen_diff <- setdiff(gens2add, gens_tree) %>% sort # get the missing genera in the target tree that need to be imputed first
      
      tps_gen_impt <- c() ; impt_tps_vec <- c()
      if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be imputed first
        for(j in 1:length(gen_diff)){
          nd_where <- node2where(gen = gen_diff[j], tree = targ_tr_ls[[k]]) # get the nodes to impute genus
          tps_2add <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(tip_names) # get tips to impute per genus
          
          # impute the first tip per genus to anchor the clade
          if(length(tps_2add) > 1) {  # in case of more than one tip
            targ_tr_ls[[k]] <- TreeTools::AddTip(tree = targ_tr_ls[[k]], 
                                      where = c(Ntip(targ_tr_ls[[k]]) + nd_where), label = tps_2add[1]) %>% 
              castor::extend_tree_to_height() %>% .$tree
            for(l in 2:length(tps_2add)){ # then impute the other tips in the same genus
              tem_tr <- add.species.to.genus(tree = targ_tr_ls[[k]], 
                                             species = tps_2add[l], where = where2impt)
              if(is.ultrametric(tem_tr)) {targ_tr_ls[[k]] <- tem_tr} else {targ_tr_ls[[k]] <- force.ultrametric(tem_tr, method = "nnls", message = F)}
            }
          } else { # in case of one tip per genus
            targ_tr_ls[[k]] <- TreeTools::AddTip(tree = targ_tr_ls[[k]], 
                                                 where = c(Ntip(targ_tr_ls[[k]]) + nd_where), label = tps_2add) %>% 
              castor::extend_tree_to_height() %>% .$tree
          }
          impt_tps_vec[[k]] <- c(impt_tps_vec[[k]], tps_2add)
        }
        
        tips2add_ls[[k]] <- setdiff(tips2add_ls[[k]], impt_tps_vec[[k]]) %>% sort
        for(l in 1:length(tips2add_ls[[k]])){
          tem_tr <- add.species.to.genus(tree = targ_tr_ls[[k]], 
                                         species = tips2add_ls[[k]][l], where = where2impt)
          if(is.ultrametric(tem_tr)) {targ_tr_ls[[k]] <- tem_tr} else {targ_tr_ls[[k]] <- force.ultrametric(tem_tr, method = "nnls", message = F)}
        }
        
      } else { # otherwise impute species to genus already in the tree
        for(l in 1:length(tips2add_ls[[k]])){
          tem_tr <- add.species.to.genus(tree = targ_tr_ls[[k]], 
                                         species = tips2add_ls[[k]][l], where = where2impt)
          if(is.ultrametric(tem_tr)) {targ_tr_ls[[k]] <- tem_tr} else {targ_tr_ls[[k]] <- force.ultrametric(tem_tr, method = "nnls", message = F)}
          
        }
      }
      ifelse(targ_tr_ls[[k]]$tip.label %in% fern_tre$tip.label %>% sum == length(fern_tre$tip.label), 
             print(paste0("*** All tips imputed ***")), print(paste0("WARNING ! unequal number of tips")))
      print(paste0("-> tree: ", k))
    }
    print(paste0("## file: ", i))
    targ_tr_ls %>% 
      saveRDS(paste0("output/imputed_trees/", out_fold,"/late_pr_tax_imp_", where2impt, "_", late_frc[i], ".rds"))
  } 
} else {
  for(i in 1:length(late_fls)){
    targ_tr_ls <- 
      readRDS(paste0("output/imputed_trees/", out_fold,"/late_pr_tax_imp_", where2impt, "_",  late_frc[i], ".rds"))
    sapply(targ_tr_ls, Ntip) %>% table() %>% print
  }
}

sapply(targ_tr_ls, Ntip)
sapply(targ_tr_ls, is.rooted) #%>% table()
sapply(targ_tr_ls, is.ultrametric)# %>% table()

# marginal pruned
if(run_marg_tax){ 
  for(i in 1:length(marg_fls)){ 
    trees_ls <- readRDS(marg_fls[[i]]) # load set of trees
    for(h in 1:length(trees_ls)){trees_ls[[h]]$node.label <- c()}
    trees_ls <- mclapply(trees_ls, name_nodes, m_tax = tax_tbl, mc.cores = n_cores) # name nodes 
    
    targ_tr_ls <- trees_ls # set trees to impute
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) # get missing tips per tree
    
    for(k in 1:length(targ_tr_ls)){
      gens2add <- strsplit(tips2add_ls[[k]], '_', fixed = F) %>% # get the genera from the tips to be added in the target tree
        sapply("[", 1) %>% unique() %>% sort()
      gens_tree <- strsplit(targ_tr_ls[[k]]$tip.label, '_', fixed = F) %>% # get the genera in the target tree
        sapply("[", 1) %>% unique() %>% sort
      
      gen_diff <- setdiff(gens2add, gens_tree) %>% sort # get the missing genera in the target tree that need to be imputed first
      
      tps_gen_impt <- c() ; impt_tps_vec <- c()
      if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be imputed first
        for(j in 1:length(gen_diff)){
          nd_where <- node2where(gen = gen_diff[j], tree = targ_tr_ls[[k]]) # get the nodes to impute genus
          tps_2add <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(tip_names) # get tips to impute per genus
          
          # impute the first tip per genus to anchor the clade
          if(length(tps_2add) > 1) {  # in case of more than one tip
            targ_tr_ls[[k]] <- TreeTools::AddTip(tree = targ_tr_ls[[k]], 
                                      where = c(Ntip(targ_tr_ls[[k]]) + nd_where), label = tps_2add[1]) %>% 
              castor::extend_tree_to_height() %>% .$tree
            for(l in 2:length(tps_2add)){ # then impute the other tips in the same genus
              tem_tr <- add.species.to.genus(tree = targ_tr_ls[[k]], 
                                             species = tps_2add[l], where = where2impt)
              if(is.ultrametric(tem_tr)) {targ_tr_ls[[k]] <- tem_tr} else {targ_tr_ls[[k]] <- force.ultrametric(tem_tr, method = "nnls", message = F)}
            }
          } else { # in case of one tip per genus
            targ_tr_ls[[k]] <- TreeTools::AddTip(tree = targ_tr_ls[[k]], 
                                                 where = c(Ntip(targ_tr_ls[[k]]) + nd_where), label = tps_2add) %>% 
              castor::extend_tree_to_height() %>% .$tree
          }
          impt_tps_vec[[k]] <- c(impt_tps_vec[[k]], tps_2add)
        }
        
        tips2add_ls[[k]] <- setdiff(tips2add_ls[[k]], impt_tps_vec[[k]]) %>% sort
        for(l in 1:length(tips2add_ls[[k]])){
          tem_tr <- add.species.to.genus(tree = targ_tr_ls[[k]], 
                                         species = tips2add_ls[[k]][l], where = where2impt)
          if(is.ultrametric(tem_tr)) {targ_tr_ls[[k]] <- tem_tr} else {targ_tr_ls[[k]] <- force.ultrametric(tem_tr, method = "nnls", message = F)}
        }
        
      } else { # otherwise impute species to genus already in the tree
        for(l in 1:length(tips2add_ls[[k]])){
          tem_tr <- add.species.to.genus(tree = targ_tr_ls[[k]], 
                                         species = tips2add_ls[[k]][l], where = where2impt)
          if(is.ultrametric(tem_tr)) {targ_tr_ls[[k]] <- tem_tr} else {targ_tr_ls[[k]] <- force.ultrametric(tem_tr, method = "nnls", message = F)}
          
        }
      }
      ifelse(targ_tr_ls[[k]]$tip.label %in% fern_tre$tip.label %>% sum == length(fern_tre$tip.label), 
             print(paste0("All tips imputed")), print(paste0("WARNING ! unequal number of tips")))
      print(paste0("-> tree: ", k))
    }
    print(paste0("## file: ", i))
    targ_tr_ls %>% 
      saveRDS(paste0("output/imputed_trees/", out_fold,"/marg_pr_tax_imp_", where2impt, "_", marg_frc[i], ".rds"))
  }
} else {
  for(i in 1:length(marg_fls)){
    targ_tr_ls <- 
      readRDS(paste0("output/imputed_trees/", out_fold,"/marg_pr_tax_imp_", where2impt, "_",  marg_frc[i], ".rds"))
    sapply(targ_tr_ls, Ntip) %>% table() %>% print
  }
}

sapply(targ_tr_ls, Ntip)
sapply(targ_tr_ls, is.rooted)# %>% table()
sapply(targ_tr_ls, is.ultrametric)# %>% table()

print(paste0("END ", where2impt, " imputation"))
  }
}

# check tips were imputed ok ####
k = 6
tips2add_ls[[k]] %>% sort()
clade_name <- "Equisetum" #"Dryopteris" #"Polypodiaceae"

tips2add_ls[[k]] %in% targ_tr_ls[[k]]$tip.label %>% table()

tt <- grep(pattern = clade_name, x = tips2add_ls[[k]], value = T) # the tips that should be imputed

tt %in% targ_tr_ls[[k]]$tip.label

# mrca clade in the new tree
grep(pattern = clade_name, x = targ_tr_ls[[k]]$tip.label) # see the nodes/tips

mrca_node <- targ_tr_ls[[k]]$tip.label[grep(pattern = clade_name, x = targ_tr_ls[[k]]$tip.label)] %>% 
  getMRCA(phy = targ_tr_ls[[k]], tip = .)

mrca_node

clad_extr <- extract.clade(phy = targ_tr_ls[[k]], node = mrca_node)
clad_extr %>% plot(no.margin = T, cex = .5)
nodelabels(text = clad_extr$node.label, col = "red", frame = "none", cex = .6)

# extract clade by node label
n_node <- grep(pattern = clade_name, x = targ_tr_ls[[k]]$node.label)
extract.clade(phy = targ_tr_ls[[k]], node = Ntip(targ_tr_ls[[k]])+n_node) %>% 
  plot(no.margin = T, cex = .75)

# sub tree in original and imputed
par(mfrow = c(1, 2))
fern_tre$tip.label[grep(pattern = clade_name, x = fern_tre$tip.label)] %>% 
  getMRCA(phy = fern_tre, tip = .) %>% extract.clade(phy = fern_tre, node = .) %>% 
  plot(no.margin = T, cex = .5, main = "\nFull tree")
clad_extr %>% plot(no.margin = T, cex = .5,  main = "\nImputed tree")
par(mfrow = c(1, 1))

# original/complete tree
par(mfrow = c(1, 3))
plot(fern_tre, no.margin = T, cex = .0005, main = "\nFull tree")
#nodelabels(text = fern_tre$node.label, col = "red", frame = "none")

# prunned an re-named tree before imputing
plot(trees_ls[[k]], no.margin = T, cex = .0005, main = "\nPruned tree")
#nodelabels(text = trees_ls[[k]]$node.label, col = "red", frame = "none")

# prunned an re-named tree before imputing
plot(targ_tr_ls[[k]], no.margin = T, cex = .0005, main = "\nImputed tree")
#nodelabels(text = trees_ls[[k]]$node.label, col = "red", frame = "none")
par(mfrow = c(1, 1))
