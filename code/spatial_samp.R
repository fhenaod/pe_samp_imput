library(tidyverse)
library(geiger)
library(parallel)

source("code/pe_samp_funs.R")

# load data ####
# load community matrix
fern_com <- read.csv("data/nitta-et-al/japan_ferns_comm_full.csv", header = T) %>% 
  column_to_rownames("grids")

# load tree
fern_tre <- read.tree("output/nitta_node-named.tre") 
is.binary(fern_tre) ; is.ultrametric(fern_tre)

# load taxonomy table
tax_tbl <- read.csv("output/fern_taxonomy_ncbi.csv") %>% 
  mutate(tip_names = sub(pattern = "(\\_.*?)\\_", "\\1-", tip_names))

# spatial data
spp_dist_df <- fern_com %>% colSums() %>% data.frame(n_grid = .) %>% 
  rownames_to_column("tip_names") %>% 
  mutate(prop_grid = n_grid/dim(fern_com)[1],
         ln_prop_grd = log(prop_grid))

spp_dist_df %>% 
  gather(key, value, -tip_names) %>% 
  ggplot(aes(x = (value))) + geom_histogram() +
  facet_wrap(~ key, scales = "free")

# plot geogr range size hist
spp_dist_df %>% 
ggplot(aes(x = (n_grid))) + 
  geom_histogram(color = "white", fill = "gray10") +
  theme_classic(base_size = 15) +
  labs(title = "Ferns of Japan", 
       x = "Geographic range size", y = "Number of Species")
ggsave("figures/hist_fern_rang.png", width = 15, height = 15, units = "cm")

# general analysis controls ####
mc_cores <- 4

prune_trs <- F
imput_trs <- F
rn_canape <- F

# load simulated trees
yng_trees <- readRDS("output/sim_trees/sim_trees_100_age.rds")
old_trees <- readRDS("output/sim_trees/sim_trees_300_age.rds")
cnt_trees <- readRDS("output/sim_trees/sim_trees_217_age.rds")
exp_trees <- readRDS("output/sim_trees/sim_trees_exp_brlen.rds")

trs_ls <- list(fern_tre, yng_trees, old_trees, cnt_trees, exp_trees)
names(trs_ls) <- c("fern", "young", "old", "cnt", "exp")

names2load <- c("fern", "young", "old", "cnt", "exp")

# pruning ####
if(prune_trs){
dir.create("output/pruned_trees/")

for(m in 1:length(trs_ls)){
  targ_tre <- trs_ls[[m]]
  out_fold <- names(trs_ls)[[m]]
  
  dir.create(paste0("output/pruned_trees/", out_fold))

## lateral pruning ####
# classifies continuous values into five percentiles, and samples/prunes the first four in increments
lat_samp_prop <- c(.05, .15, .3, .5) 
lat_cat2filt <- c("below.05", ".05to.15", ".15to.3", ".3to.5", "above.5")

lat_tips2drop <- list()
for(j in 1:length(lat_samp_prop)){
  lat_tips2drop <- spp_dist_df %>% 
    select(tip_names, ln_prop_grd) %>% 
    mutate(across(where(is.numeric), lat_funct)) %>% 
    filter(ln_prop_grd %in% lat_cat2filt[1:j]) %>% pull(tip_names)
  
  tr_class <- class(targ_tre) ; print(tr_class)
  lat_trs <- switch(as.character(tr_class), 
                    list = mclapply(X = targ_tre, function(X) drop.tip(phy = X, tip = lat_tips2drop), mc.cores = mc_cores), 
                    phylo = drop.tip(phy = targ_tre, tip = lat_tips2drop)
                    )
  
  # print n tips
  switch(as.character(tr_class),
         list = print(sapply(lat_trs, Ntip)),
         phylo = print(Ntip(lat_trs))
         )
  
  lat_trs %>% 
    saveRDS(paste0("output/pruned_trees/", out_fold, "/lateral_spat_pruned_trs_", lat_samp_prop[j], ".rds"))
}

## marginal pruning ####
marg_samp_prop <- seq(from = .025, to = 0.975, by = .05)
cat2filt <- paste0(rep("perc_", length(marg_samp_prop)+1), c(marg_samp_prop, 1))

marg_cat2filt_ls <- list(c(head(cat2filt, 1), tail(cat2filt, 1)),
                         c(head(cat2filt, 2), tail(cat2filt, 2)),
                         c(head(cat2filt, 3), tail(cat2filt, 3)),
                         c(head(cat2filt, 5), tail(cat2filt, 5))
                         )

marg_trs <- list()
for(i in 1:length(marg_cat2filt_ls)){
  marg_tips2drop <- spp_dist_df %>%
    select(tip_names, ln_prop_grd) %>% 
    mutate(ln_prop_grd = margin_funct(ln_prop_grd)) %>% 
    filter(ln_prop_grd %in% marg_cat2filt_ls[[i]]) %>% pull(tip_names)
  
  tr_class <- class(targ_tre) ; print(tr_class)
  marg_trs <- switch(as.character(tr_class),
                    list = mclapply(X = targ_tre, function(X) drop.tip(phy = X, tip = marg_tips2drop), mc.cores = mc_cores),
                    phylo = drop.tip(phy = targ_tre, tip = marg_tips2drop)
                    )
  
  # print n tips
  switch(as.character(tr_class),
         list = print(sapply(marg_trs, Ntip)),
         phylo = print(Ntip(marg_trs))
  )
  marg_trs %>% 
    saveRDS(paste0("output/pruned_trees/", out_fold, "/marginal_spat_pruned_trs_", c(length(marg_cat2filt_ls[[i]])*.05), ".rds"))
  }
 }
} else {print("Spatial pruning skiped")}

# imputation ####

run_late_rand <- F ; run_marg_rand <- F # spatially pruned and randomly imputed trees
run_late_tax <- T  ; run_marg_tax <- T  # spatially pruned and taxonomically imputed trees

if(imput_trs){
# loop over trees folders, sampling fractions and impute tips 
for(m in 1:length(names2load)){
  to_load <- names2load[m]
  out_fold <- to_load
  
  late_fls <- list.files(paste0("output/pruned_trees/", to_load), "lateral_spat",  all.files = T, full.names = T)
  marg_fls <- list.files(paste0("output/pruned_trees/", to_load), "marginal_spat", all.files = T, full.names = T)
  
  late_frc <- readr::parse_number(late_fls)
  marg_frc <- readr::parse_number(marg_fls)
  
  dir.create(paste0("output/imputed_trees/", out_fold))

## add tips at random ####
# lateral pruned
if(run_late_rand) {
  for(j in 1:length(late_fls)){
    targ_tr_ls <- readRDS(late_fls[[j]])
    
    tr_class <- class(targ_tr_ls) ; print(tr_class)
    # get tips to impute
    tips2add_ls <- switch(as.character(tr_class), 
                          list = mclapply(X = targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label), mc.cores = mc_cores), 
                          phylo = setdiff(fern_tre$tip.label, targ_tr_ls$tip.label)
                          )
      
    # impute tips
    impt_tr_rmd <- switch(as.character(tr_class), 
                      list = mcmapply(add.random, tree = targ_tr_ls, tips = tips2add_ls, mc.cores = mc_cores, SIMPLIFY = F),
                      phylo = add.random(tree = targ_tr_ls, tips = tips2add_ls)
                      )
    
    # print n tips
    switch(as.character(tr_class),
           list = print(sapply(impt_tr_rmd, Ntip)),
           phylo = print(Ntip(impt_tr_rmd))
           )
    
    impt_tr_rmd %>% 
      saveRDS(paste0("output/imputed_trees/", out_fold, "/late_spat_pr_rand_imp_", late_frc[j], ".rds"))
  }
} else {print("Random imputation skiped on lateral pruning")}

# marginal pruned
if(run_marg_rand){
  for(j in 1:length(marg_fls)){
    targ_tr_ls <- readRDS(marg_fls[[j]])
    
    tr_class <- class(targ_tr_ls) ; print(tr_class)
    # get tips to impute
    tips2add_ls <- switch(as.character(tr_class), 
                          list = mclapply(X = targ_tr_ls, function(X) setdiff(fern_tre$tip.label, X$tip.label), mc.cores = mc_cores), 
                          phylo = setdiff(fern_tre$tip.label, targ_tr_ls$tip.label)
                          )
    
    # impute tips
    impt_tr_rmd <- switch(as.character(tr_class), 
                          list = mcmapply(add.random, tree = targ_tr_ls, tips = tips2add_ls, mc.cores = mc_cores, SIMPLIFY = F),
                          phylo = add.random(tree = targ_tr_ls, tips = tips2add_ls)
                          )
    
    # print n tips
    switch(as.character(tr_class),
           list = print(sapply(impt_tr_rmd, Ntip)),
           phylo = print(Ntip(impt_tr_rmd))
           )
    
    impt_tr_rmd %>% 
      saveRDS(paste0("output/imputed_trees/", out_fold, "/marg_spat_pr_rand_imp_", marg_frc[j], ".rds"))
  }
} else {print("Random imputation skiped on marginal pruning")}

## add tips to based on taxonomy ####

meth2impt <- c("root", "random") 

for(f in 1:length(meth2impt)){
  where2impt <- meth2impt[f]
  
  # lateral pruned
  if(run_late_tax){ 
    for(i in 1:length(late_fls)){ 
      trees_ls <- readRDS(late_fls[[i]]) # load set of trees
      
      tr_class <- class(trees_ls) ; print(tr_class)
      # impute tips
      targ_tr_ls <- switch(as.character(tr_class), 
                            list = impute_list_trees(tree_list = trees_ls, n_cores = mc_cores, where2impt = meth2impt[f]),
                            phylo = impute_tree(tree = trees_ls)
                           )
      
      # print n tips
      switch(as.character(tr_class),
             list = print(sapply(targ_tr_ls, Ntip)),
             phylo = print(Ntip(targ_tr_ls))
             )
      
      print(paste0("## file: ", i))
      targ_tr_ls %>% 
        saveRDS(paste0("output/imputed_trees/", out_fold, "/late_spat_pr_tax_imp_", where2impt, "_", late_frc[i], ".rds"))
    }
  } else {
    for(i in 1:length(late_fls)){
      targ_tr_ls <- 
        readRDS(paste0("output/imputed_trees/", out_fold, "/late_spat_pr_tax_imp_", where2impt, "_",  late_frc[i], ".rds"))
    }
  }
  
  # marginal pruned
  if(run_marg_tax){ 
    for(i in 1:length(marg_fls)){ 
      trees_ls <- readRDS(marg_fls[[i]]) # load set of trees
      
      tr_class <- class(trees_ls) ; print(tr_class)
      # impute tips
      targ_tr_ls <- switch(as.character(tr_class), 
                           list = impute_list_trees(tree_list = trees_ls, n_cores = mc_cores, where2impt = meth2impt[f]),
                           phylo = impute_tree(tree = trees_ls)
                           )
      
      # print n tips
      switch(as.character(tr_class),
             list = print(sapply(targ_tr_ls, Ntip)),
             phylo = print(Ntip(targ_tr_ls))
             )
       
      print(paste0("## file: ", i))
      targ_tr_ls %>% 
        saveRDS(paste0("output/imputed_trees/", out_fold, "/marg_spat_pr_tax_imp_", where2impt, "_", marg_frc[i], ".rds"))
    }
  } else {
    for(i in 1:length(marg_fls)){
      targ_tr_ls <- 
        readRDS(paste0("output/imputed_trees/", out_fold, "/marg_spat_pr_tax_imp_", where2impt, "_",  marg_frc[i], ".rds"))
    }
  }
  
  print(paste0("END ", where2impt, " imputation"))
  }
}
} else {print("Tips imputation skiped")}

# run canaper ####
library(canaper)
library(future)

fern_com <- read.csv("data/nitta-et-al/japan_ferns_comm_full.csv", header = T) %>% 
  column_to_rownames("grids")

colnames(fern_com) <- gsub(pattern = "\\.", replacement = "-", colnames(fern_com))
colnames(fern_com) <- sub(pattern = "(\\_.*?)\\_", "\\1-", colnames(fern_com))

# missing taxa in tree
fern_com %>% colnames() %in% fern_tre$tip.label %>% table()

## canaper control pars ####
run_iter <- F
run_late <- T ; run_marg <- T # pruned trees
run_late_rand <- T ; run_marg_rand <- T # pruned and randomly imputed trees
run_late_tax <- T  ; run_marg_tax <- T  # pruned and taxonomically imputed trees

canaper_n_reps <- 10 #1000 

n_cores <- floor(parallel::detectCores()*.75) # 3
plan(multisession, workers = n_cores)

# re-name data
targ_tr <- fern_tre
occ_mat <- fern_com

# run iterations to find saturation point for randomization
if(run_iter){
  iter_sim_res <- cpr_iter_sim(comm = occ_mat, 
                               null_model = "curveball", 
                               n_iterations = 20000, thin = 200)
  head(iter_sim_res) ; tail(iter_sim_res)
  iter_sim_res %>% ggplot(aes(x = iteration, y = similarity)) + geom_line()
  
  dydx <- diff(iter_sim_res$similarity) / diff(iter_sim_res$iteration)
  saturation_point <- iter_sim_res$iteration[which.min(abs(dydx))]
} else {
  saturation_point <- 10 #10600
}

if(rn_canape){
for(m in 1:length(names2load)){
  to_load <- names2load[m]
  out_fold <- to_load

## only pruned ####
dir.create(paste0("output/canaper/", out_fold))

# lateral pruning 
late_fls <- list.files(paste0("output/pruned_trees/", to_load), "lateral_spat",  all.files = T, full.names = T)
late_frc <- readr::parse_number(late_fls)

if(run_late){
  for(i in 1:length(late_fls)){
    late_canaper_ls <- list()
    targ_tr_ls <- readRDS(late_fls[i])
    
    tr_class <- class(targ_tr_ls) ; print(tr_class)
    late_canaper_ls <- switch(as.character(tr_class), 
                          list = run_canape_ls_trees(targ_tr_ls), 
                          phylo = run_cananpe_tree(targ_tr_ls)
                          )

    late_canaper_ls %>% 
      saveRDS(paste0("output/canaper/", out_fold, "/canaper_late_spat_pr_", late_frc[i],".rds"))
    gc()
  }
} else {} 

# marginal pruning 
marg_fls <- list.files(paste0("output/pruned_trees/", to_load), "marginal_spat", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_fls)

if(run_marg){
  for(i in 1:length(marg_fls)){
    marg_canaper_ls <- list()
    targ_tr_ls <- readRDS(marg_fls[i])
    
    tr_class <- class(targ_tr_ls) ; print(tr_class)
    marg_canaper_ls <- switch(as.character(tr_class), 
                              list = run_canape_ls_trees(targ_tr_ls), 
                              phylo = run_cananpe_tree(targ_tr_ls)
                              )
      
    marg_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
    print(paste0("iteration: ", i))
    
    marg_canaper_ls %>% 
      saveRDS(paste0("output/canaper/", out_fold, "/canaper_marg_spat_pr_", marg_frc[i],".rds"))
    gc()
  }
} else {}

## pruned and randomly imputed ####
# lateral pruning
pttr <- "late_spat_pr_rand_imp_"
late_fls <- list.files(paste0("output/imputed_trees/", to_load), pttr, all.files = T, full.names = T)
late_frc <- readr::parse_number(late_fls)

if(run_late_rand){
  for(i in 1:length(late_fls)){
    late_canaper_ls <- list()
    targ_tr_ls <- readRDS(late_fls[i])
    
    tr_class <- class(targ_tr_ls) ; print(tr_class)
    late_canaper_ls <- switch(as.character(tr_class), 
                              list = run_canape_ls_trees(targ_tr_ls), 
                              phylo = run_cananpe_tree(targ_tr_ls)
                              )
      
    late_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
    print(paste0("iteration: ", i))

    late_canaper_ls %>% 
      saveRDS(paste0("output/canaper/", out_fold, "/canaper_", pttr, late_frc[i],".rds"))
    gc()
  }
} else {}

# marginal pruning
pttr <- "marg_spat_pr_rand_imp_"
marg_fls <- list.files(paste0("output/imputed_trees/", to_load), pttr, all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_fls)

if(run_marg_rand){
  for(i in 1:length(marg_fls)){
    marg_canaper_ls <- list()
    targ_tr_ls <- readRDS(marg_fls[i])

    tr_class <- class(targ_tr_ls) ; print(tr_class)
    marg_canaper_ls <- switch(as.character(tr_class), 
                              list = run_canape_ls_trees(targ_tr_ls), 
                              phylo = run_cananpe_tree(targ_tr_ls)
                              )
      
    marg_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
    print(paste0("iteration: ", i))
    
    marg_canaper_ls %>% 
      saveRDS(paste0("output/canaper/", out_fold, "/canaper_", pttr, marg_frc[i],".rds"))
    gc()
  }
} else {}

## pruned and taxonomic imputed ####
# lateral pruning
pttr <- c("late_spat_pr_tax_imp_root", "late_spat_pr_tax_imp_random")

for(f in 1:length(pttr)){
  late_fls <- list.files(paste0("output/imputed_trees/", to_load), pttr[f], all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_fls)
  if(run_late_tax){
    for(i in 1:length(late_fls)){
      late_canaper_ls <- list()
      targ_tr_ls <- readRDS(late_fls[i])
      
      tr_class <- class(targ_tr_ls) ; print(tr_class)
      late_canaper_ls <- switch(as.character(tr_class), 
                                list = run_canape_ls_trees(targ_tr_ls), 
                                phylo = run_cananpe_tree(targ_tr_ls)
                                )
        
      late_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i))
      
      late_canaper_ls %>% 
        saveRDS(paste0("output/canaper/", out_fold, "/canaper_", pttr[f], "_", late_frc[i],".rds"))
      gc()
    }
  } else {}
  cat(paste0("END: ", pttr[f]))
}

# marginal pruning
pttr <- c("marg_spat_pr_tax_imp_root", "marg_spat_pr_tax_imp_random")
for(f in 1:length(pttr)){
  marg_fls <- list.files(paste0("output/imputed_trees/", to_load), pttr[f], all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_fls)
  if(run_marg_tax){
    for(i in 1:length(marg_fls)){
      marg_canaper_ls <- list()
      targ_tr_ls <- readRDS(marg_fls[i])
      
      tr_class <- class(targ_tr_ls) ; print(tr_class)
      marg_canaper_ls <- switch(as.character(tr_class), 
                                list = run_canape_ls_trees(targ_tr_ls), 
                                phylo = run_cananpe_tree(targ_tr_ls)
                                )
        
      marg_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i))
      
      marg_canaper_ls %>% 
        saveRDS(paste0("output/canaper/", out_fold, "/canaper_", pttr[f], "_",  marg_frc[i],".rds"))
      gc()
    }
  } else {}
  cat(paste0("END: ", pttr[f]))
    }
  }
} else {print("CANAPER skiped")}

