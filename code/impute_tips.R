library(tidyverse)
library(phytools) 
library(parallel)

source("code/pe_samp_funs.R")

# load full empirical tree
fern_tre <- read.tree("output/nitta_node-named.tre") 

# load taxonomy table
tax_tbl <- read.csv("output/fern_taxonomy_ncbi.csv") %>% 
  mutate(tip_names = sub(pattern = "(\\_.*?)\\_", "\\1-", tip_names))

# load pruned trees
rand_fls <- list.files("output/pruned_trees", "rand_pruned",     all.files = T, full.names = T)
late_fls <- list.files("output/pruned_trees", "lateral_pruned",  all.files = T, full.names = T)
marg_fls <- list.files("output/pruned_trees", "marginal_pruned", all.files = T, full.names = T)

rand_frc <- readr::parse_number(rand_fls)
late_frc <- readr::parse_number(late_fls)
marg_frc <- readr::parse_number(marg_fls)

dir.create("output/imputed_trees/")

run_rand_rand <- F ; run_late_rand <- F ; run_marg_rand <- F # pruned and randomly imputed trees
run_rand_tax <- T  ; run_late_tax <- T  ; run_marg_tax <- T  # pruned and taxonomically imputed trees
meth2impt <- c("root", "random") 
n_cores <- 3

# add tips location at random ####
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
    impt_tr_rmd %>% saveRDS(paste0("output/imputed_trees/rand_pr_rand_imp_", rand_frc[j], ".rds"))
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
    impt_tr_rmd %>% saveRDS(paste0("output/imputed_trees/late_pr_rand_imp_", late_frc[j], ".rds"))
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
    impt_tr_rmd %>% saveRDS(paste0("output/imputed_trees/marg_pr_rand_imp_", marg_frc[j], ".rds"))
  }
}

# add tips to based on taxonomy ####
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
             print(paste0("All tips imputed")), print(paste0("WARNING ! unequal number of tips")))
      print(paste0("-> tree: ", k))
    }
    print(paste0("## file: ", i))
    targ_tr_ls %>% 
      saveRDS(paste0("output/imputed_trees/rand_pr_tax_imp_", where2impt, "_", rand_frc[i], ".rds"))
  }
} else {
  for(i in 1:length(rand_fls)){
    targ_tr_ls <- readRDS(paste0("output/imputed_trees/rand_pr_tax_imp_", where2impt, "_",  rand_frc[i], ".rds"))
    sapply(targ_tr_ls, Ntip) %>% table() %>%  print
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
             print(paste0("All tips imputed")), print(paste0("WARNING ! unequal number of tips")))
      print(paste0("-> tree: ", k))
    }
    print(paste0("## file: ", i))
    targ_tr_ls %>% 
      saveRDS(paste0("output/imputed_trees/late_pr_tax_imp_", where2impt, "_", late_frc[i], ".rds"))
  }
} else {
  for(i in 1:length(late_fls)){
    targ_tr_ls <- readRDS(paste0("output/imputed_trees/late_pr_tax_imp_", where2impt, "_",  late_frc[i], ".rds"))
    sapply(targ_tr_ls, Ntip) %>% table() %>%  print
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
      saveRDS(paste0("output/imputed_trees/marg_pr_tax_imp_", where2impt, "_", marg_frc[i], ".rds"))
  }
} else {
  for(i in 1:length(marg_fls)){
    targ_tr_ls <- readRDS(paste0("output/imputed_trees/marg_pr_tax_imp_", where2impt, "_",  marg_frc[i], ".rds"))
    sapply(targ_tr_ls, Ntip) %>% table() %>%  print
  }
}

sapply(targ_tr_ls, Ntip)
sapply(targ_tr_ls, is.rooted)# %>% table()
sapply(targ_tr_ls, is.ultrametric)# %>% table()

print(paste0("END ", where2impt, " imputation"))
}

# vis comparative trees ####
library(ggtree) ; library(patchwork)
# pruned trees
rand_prun_tree <- readRDS(rand_fls[[3]])[[10]] 
late_prun_tree <- readRDS(late_fls[[3]])[[10]]
marg_prun_tree <- readRDS(marg_fls[[4]])[[10]]

# taxonomic at node imputation
rand_tax_imp_tree <- readRDS(paste0("output/imputed_trees/rand_pr_tax_imp_root_",  rand_frc[3], ".rds"))[[10]]
late_tax_imp_tree <- readRDS(paste0("output/imputed_trees/late_pr_tax_imp_root_",  late_frc[3], ".rds"))[[10]]
marg_tax_imp_tree <- readRDS(paste0("output/imputed_trees/marg_pr_tax_imp_root_",  marg_frc[4], ".rds"))[[10]]

# taxonomic at random within
rand_tax_imp_rand_tree <- readRDS(paste0("output/imputed_trees/rand_pr_tax_imp_random_",  rand_frc[3], ".rds"))[[10]]
late_tax_imp_rand_tree <- readRDS(paste0("output/imputed_trees/late_pr_tax_imp_random_",  late_frc[3], ".rds"))[[10]]
marg_tax_imp_rand_tree <- readRDS(paste0("output/imputed_trees/marg_pr_tax_imp_random_",  marg_frc[4], ".rds"))[[10]]

typ_tre <- "rectangular" # "fan" #"roundrect" # 

# full tree
fern_full <- fern_tre %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) + 
  geom_nodepoint(size = .25, alpha = .55) +
  labs(title = "Full tree")

# Randomly pruned
rand_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, rand_prun_tree$tip.label), 
                    keept = intersect(fern_tre$tip.label, rand_prun_tree$tip.label))

t1.1 <- fern_tre %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .15, alpha = .55)  +
  labs(title = "Random pruning (30%)") #+ geom_tiplab(aes(subset = label %in% grp1$pruned), color = 'firebrick', size = 1.2)
 
ran_prn <- groupOTU(t1.1, rand_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "red"))

# Randomly pruned and taxonomically imputed (at node)
rand_prn_ls <- list(imputed = setdiff(fern_tre$tip.label, rand_prun_tree$tip.label), 
                    keept = intersect(fern_tre$tip.label, rand_prun_tree$tip.label))
t2.1 <- rand_tax_imp_tree %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .15, alpha = .55)  +
  labs(title = "Random pruning (30%)", subtitle = "Taxonomic imputed (node)") #+ geom_tiplab(aes(subset = label %in% grp2$imputed), color = 'firebrick', size = 1.2) #+ scale_x_reverse()

ran_prn_nd_imp <- groupOTU(t2.1, rand_prn_ls, 'tips') + aes(color = tips) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = c("red", "black")) 

# Randomly pruned and taxonomically imputed (randomly within)
t3.1 <- rand_tax_imp_rand_tree %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .15, alpha = .55)  +
  labs(title = "Random pruning (30%)", subtitle = "Taxonomic imputed (binary)") #+ geom_tiplab(aes(subset = label %in% grp4$imputed), color = 'firebrick', size = 1.2)

ran_prn_bin <- groupOTU(t3.1, rand_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "black"))

# Lateraly pruned 
late_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, late_prun_tree$tip.label), 
                    keept = intersect(fern_tre$tip.label, late_prun_tree$tip.label))

t1.2 <- fern_tre %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .15, alpha = .55)  +
  labs(title = "Lateral pruning (30%)") #+ geom_tiplab(aes(subset = label %in% grp1$pruned), color = 'firebrick', size = 1.2)

late_prn <- groupOTU(t1.2, late_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "red"))

# Lateraly pruned and taxonomically imputed (at node)
late_prn_ls <- list(imputed = setdiff(fern_tre$tip.label, late_prun_tree$tip.label), 
                    keept = intersect(fern_tre$tip.label, late_prun_tree$tip.label))

t2.2 <- late_tax_imp_tree %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .15, alpha = .55)  +
  labs(title = "Lateral pruning (30%)", subtitle = "Taxonomic imputed (node)")

late_prn_nd_imp <- groupOTU(t2.2, late_prn_ls, 'tips') + aes(color = tips) + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "black"))

# Lateraly pruned and taxonomically imputed (randomly within)
t3.2 <- late_tax_imp_rand_tree %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .15, alpha = .55)  +
  labs(title = "Lateral pruning (30%)", subtitle = "Taxonomic imputed (binary)") #+ geom_tiplab(aes(subset = label %in% grp4$imputed), color = 'firebrick', size = 1.2)

late_prn_bin <- groupOTU(t3.2, late_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "black"))

if(typ_tre == "fan"){
  pl_lay <- "
  AABB
  AABB
  CCDD
  CCDD
  "
} else {
  pl_lay <- "
  AA#BBCCDD
  AA#BBCCDD
  AA#BBCCDD
  "
}

# random pruning
free(fern_full) + free(ran_prn) + free(ran_prn_nd_imp) + free(ran_prn_bin) + 
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/comp_trees_rand_prn.pdf"), 
       width = 30, height = 15, units = "cm", dpi = 200)
# lateral pruning 
free(fern_full) + free(late_prn) + free(late_prn_nd_imp) + free(late_prn_bin) +
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/comp_trees_late_prn.pdf"), 
       width = 30, height = 15, units = "cm", dpi = 200)

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
