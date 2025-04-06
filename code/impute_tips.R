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
where2impt <- "root" # "random"
n_cores <- 2

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
# random pruned
if(run_rand_tax){ 
  for(i in 1:length(rand_fls)){ 
    trees_ls <- readRDS(rand_fls[[i]]) # load set of trees
      for(l in 1:length(trees_ls)){trees_ls[[l]]$node.label <- c()}
    trees_ls <- mclapply(trees_ls, name_nodes, m_tax = tax_tbl, mc.cores = n_cores) # name nodes 
    
    targ_tr_ls <- trees_ls # set trees to impute
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) # get missing tips per tree
    
    for(k in 1:length(targ_tr_ls)){
      gens2add <- strsplit(tips2add_ls[[k]], '_', fixed = F) %>% # get the genera from the tips to be added in the target tree
        sapply("[", 1) %>% unique() %>% sort()
      gens_tree <- strsplit(targ_tr_ls[[k]]$tip.label, '_', fixed = F) %>% # get the genera in the target tree
        sapply("[", 1) %>% unique() %>% sort
        
      gen_diff <- setdiff(gens2add, gens_tree) %>% sort # get the missing genera in the target tree that need to be imputed first
      
      tps_gen_impt <- c()
      if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be imputed first
        for(j in 1:length(gen_diff)){
          nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(family) %>% # get node to impute genus at family level
            unique() %>% grepl(pattern = paste0("^", . , "$"), x = targ_tr_ls[[k]]$node.label) %>% which
          
          ifelse(test = !is_empty(nd_where), yes = print(paste0("Imputed at family level")), 
                 no = { nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(order) %>% # get node to impute genus at order
                   unique() %>% grepl(pattern = paste0("^", . , "$"), x = targ_tr_ls[[k]]$node.label) %>% 
                   which ; ifelse(!is_empty(nd_where), yes = print(paste0("Imputed at order level")),
                                  no = {nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(class) %>%  # get node to impute genus at class
                                    unique() %>% grepl(pattern = paste0("^", . , "$"), x = targ_tr_ls[[k]]$node.label) %>% 
                                    which ; ifelse(test = is_empty(nd_where), 
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
        #targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = impt_tips, where = "root")
        for(l in 1:length(impt_tips)){
          targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = impt_tips[l], where = where2impt)
          }

      } else { # otherwise impute species to genus already in the tree
        #targ_tr_ls[[k]] <- add.species.to.genus(targ_tr_ls[[k]], tips2add_ls[[k]], where = "random")
        for(l in 1:length(tips2add_ls[[k]])){
          targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = tips2add_ls[[k]][l], where = where2impt)
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
      for(l in 1:length(trees_ls)){trees_ls[[l]]$node.label <- c()}
    trees_ls <- mclapply(trees_ls, name_nodes, m_tax = tax_tbl, mc.cores = n_cores) # name nodes 
    
    targ_tr_ls <- trees_ls # set trees to impute
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) # get missing tips per tree
    
    for(k in 1:length(targ_tr_ls)){
      gens2add <- strsplit(tips2add_ls[[k]], '_', fixed = F) %>% # get the genera from the tips to be added in the target tree
        sapply("[", 1) %>% unique() %>% sort()
      gens_tree <- strsplit(targ_tr_ls[[k]]$tip.label, '_', fixed = F) %>% # get the genera in the target tree
        sapply("[", 1) %>% unique() %>% sort
      
      gen_diff <- setdiff(gens2add, gens_tree) %>% sort # get the missing genera in the target tree that need to be imputed first
      
      tps_gen_impt <- c()
      if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be imputed first
        for(j in 1:length(gen_diff)){
          nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(family) %>% # get node to impute genus at family level
            unique() %>% grepl(pattern = paste0("^", . , "$"), x = targ_tr_ls[[k]]$node.label) %>% which
          
          ifelse(test = !is_empty(nd_where), yes = print(paste0("Imputed at family level")), 
                 no = { nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(order) %>% # get node to impute genus at order
                   unique() %>% grepl(pattern = paste0("^", . , "$"), x = targ_tr_ls[[k]]$node.label) %>% 
                   which ; ifelse(!is_empty(nd_where), yes = print(paste0("Imputed at order level")),
                                  no = {nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(class) %>%  # get node to impute genus at class
                                    unique() %>% grepl(pattern = paste0("^", . , "$"), x = targ_tr_ls[[k]]$node.label) %>% 
                                    which ; ifelse(test = is_empty(nd_where), 
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
        #targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = impt_tips, where = "root")
        for(l in 1:length(impt_tips)){
          targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = impt_tips[l], where = where2impt)
        }
        
      } else { # otherwise impute species to genus already in the tree
        #targ_tr_ls[[k]] <- add.species.to.genus(targ_tr_ls[[k]], tips2add_ls[[k]], where = "random")
        for(l in 1:length(tips2add_ls[[k]])){
          targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = tips2add_ls[[k]][l], where = where2impt)
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
    for(l in 1:length(trees_ls)){trees_ls[[l]]$node.label <- c()}
    trees_ls <- mclapply(trees_ls, name_nodes, m_tax = tax_tbl, mc.cores = n_cores) # name nodes 
    
    targ_tr_ls <- trees_ls # set trees to impute
    tips2add_ls <- lapply(targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label)) # get missing tips per tree
    
    for(k in 1:length(targ_tr_ls)){
      gens2add <- strsplit(tips2add_ls[[k]], '_', fixed = F) %>% # get the genera from the tips to be added in the target tree
        sapply("[", 1) %>% unique() %>% sort()
      gens_tree <- strsplit(targ_tr_ls[[k]]$tip.label, '_', fixed = F) %>% # get the genera in the target tree
        sapply("[", 1) %>% unique() %>% sort
      
      gen_diff <- setdiff(gens2add, gens_tree) %>% sort # get the missing genera in the target tree that need to be imputed first
      
      tps_gen_impt <- c()
      if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be imputed first
        for(j in 1:length(gen_diff)){
          nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(family) %>% # get node to impute genus at family level
            unique() %>% grepl(pattern = paste0("^", . , "$"), x = targ_tr_ls[[k]]$node.label) %>% which
          
          ifelse(test = !is_empty(nd_where), yes = print(paste0("Imputed at family level")), 
                 no = { nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(order) %>% # get node to impute genus at order
                   unique() %>% grepl(pattern = paste0("^", . , "$"), x = targ_tr_ls[[k]]$node.label) %>% 
                   which ; ifelse(!is_empty(nd_where), yes = print(paste0("Imputed at order level")),
                                  no = {nd_where <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(class) %>%  # get node to impute genus at class
                                    unique() %>% grepl(pattern = paste0("^", . , "$"), x = targ_tr_ls[[k]]$node.label) %>% 
                                    which ; ifelse(test = is_empty(nd_where), 
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
        #targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = impt_tips, where = "root")
        for(l in 1:length(impt_tips)){
          targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = impt_tips[l], where = where2impt)
        }
        
      } else { # otherwise impute species to genus already in the tree
        #targ_tr_ls[[k]] <- add.species.to.genus(targ_tr_ls[[k]], tips2add_ls[[k]], where = "random")
        for(l in 1:length(tips2add_ls[[k]])){
          targ_tr_ls[[k]] <- add.species.to.genus(tree = targ_tr_ls[[k]], species = tips2add_ls[[k]][l], where = where2impt)
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

# vis comparative trees ####
library(ggtree) ; library(ggpubr)
t1 <- fern_tre %>% 
  ggtree(layout = "roundrect", ladderize = F) + 
  geom_nodepoint(size = .58, alpha = .55) +
  labs(subtitle = "Full tree")
t2 <- trees_ls[[10]] %>% 
  ggtree(layout = "roundrect", ladderize = F) +
  geom_nodepoint(size = .58, alpha = .55)  +
  labs(subtitle = "Pruned tree")
t3 <- targ_tr_ls[[10]] %>% 
  ggtree(layout = "roundrect", ladderize = F) +
  geom_nodepoint(size = .58, alpha = .55)  +
  labs(subtitle = "Imputed tree")

ggarrange(t1, t2, t3, ncol = 3 #, labels = c("A", "B", "C")
)
ggsave(filename = paste0("figures/comp_trees.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)

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
