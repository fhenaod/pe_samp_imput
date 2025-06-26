library(tidyverse)
library(geiger)

source("code/pe_samp_funs.R")

# load data ####
# load community matrix
fern_com <- read.csv("data/nitta-et-al/japan_ferns_comm_full.csv", header = T) %>% 
  column_to_rownames("grids")

# load tree
fern_tre <- read.tree("output/nitta_node-named.tre") 
is.binary(fern_tre) ; is.ultrametric(fern_tre)
targ_tre <- fern_tre

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

dir.create("output/pruned_trees/")
# lateral pruning ####
# classifies continuous values into five percentiles, and samples/prunes the first four in increments
lat_samp_prop <- c(.05, .15, .3, .5) 
lat_cat2filt <- c("below.05", ".05to.15", ".15to.3", ".3to.5", "above.5")

lat_tips2drop <- list()
for(j in 1:length(lat_samp_prop)){
  lat_tips2drop <- spp_dist_df %>% 
    select(tip_names, ln_prop_grd) %>% 
    mutate(across(where(is.numeric), lat_funct)) %>% 
    filter(ln_prop_grd %in% lat_cat2filt[1:j]) %>% pull(tip_names)
  
  lat_trs <- drop.tip(phy = targ_tre, tip = lat_tips2drop) ; print(Ntip(lat_trs))
  lat_trs %>% saveRDS(paste0("output/pruned_trees/lateral_spat_pruned_trs_", lat_samp_prop[j], ".rds"))
}

# marginal pruning ####
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
  
  marg_trs <- drop.tip(phy = targ_tre, tip = marg_tips2drop) ; print(Ntip(marg_trs))
  marg_trs %>% 
    saveRDS(paste0("output/pruned_trees/marginal_spat_pruned_trs_", c(length(marg_cat2filt_ls[[i]])*.05), ".rds"))
}

# imputation ####
late_fls <- list.files("output/pruned_trees", "lateral_spat",  all.files = T, full.names = T)
marg_fls <- list.files("output/pruned_trees", "marginal_spat", all.files = T, full.names = T)

late_frc <- readr::parse_number(late_fls)
marg_frc <- readr::parse_number(marg_fls)

dir.create("output/imputed_trees/")

run_late_rand <- T ; run_marg_rand <- T # spatially pruned and randomly imputed trees
run_late_tax <- T  ; run_marg_tax <- T  # spatially pruned and taxonomically imputed trees

meth2impt <- c("root", "random") 
n_cores <- 3 #20

## add tips location at random ####
# lateral pruned
if(run_late_rand) {
  for(j in 1:length(late_fls)){
    targ_tr_ls <- readRDS(late_fls[[j]])
    tips2add_ls <- setdiff(fern_tre$tip.label, targ_tr_ls$tip.label) 
    length(targ_tr_ls) == length(tips2add_ls)
    
    impt_tr_rmd <- add.random(tree = targ_tr_ls, tips = tips2add_ls) 
    impt_tr_rmd %>% saveRDS(paste0("output/imputed_trees/late_spat_pr_rand_imp_", late_frc[j], ".rds"))
  }
}

# marginal pruned
if(run_marg_rand){
  for(j in 1:length(marg_fls)){
    targ_tr_ls <- readRDS(marg_fls[[j]])
    tips2add_ls <- setdiff(fern_tre$tip.label, targ_tr_ls$tip.label) 
    length(targ_tr_ls) == length(tips2add_ls)
    
    impt_tr_rmd <- add.random(tree = targ_tr_ls, tips = tips2add_ls) 
    impt_tr_rmd %>% saveRDS(paste0("output/imputed_trees/marg_spat_pr_rand_imp_", marg_frc[j], ".rds"))
  }
}

## add tips to based on taxonomy ####
for(f in 1:length(meth2impt)){
  where2impt <- meth2impt[f]
  
  # lateral pruned
  if(run_late_tax){ 
    for(i in 1:length(late_fls)){ 
      trees_ls <- readRDS(late_fls[[i]]) # load set of trees
      for(h in 1:length(trees_ls)){trees_ls$node.label <- c()}
      trees_ls <- name_nodes(tree = trees_ls, m_tax = tax_tbl) # name nodes 

      targ_tr_ls <- trees_ls # set trees to impute
      tips2add_ls <- setdiff(fern_tre$tip.label, targ_tr_ls$tip.label) # get missing tips per tree
      
        gens2add <- strsplit(tips2add_ls, '_', fixed = F) %>% # get the genera from the tips to be added in the target tree
          sapply("[", 1) %>% unique() %>% sort()
        gens_tree <- strsplit(targ_tr_ls$tip.label, '_', fixed = F) %>% # get the genera in the target tree
          sapply("[", 1) %>% unique() %>% sort
        
        gen_diff <- setdiff(gens2add, gens_tree) %>% sort # get the missing genera in the target tree that need to be imputed first
        
        tps_gen_impt <- c() ; impt_tps_vec <- c()
        if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be imputed first
          for(j in 1:length(gen_diff)){
            nd_where <- node2where(gen = gen_diff[j], tree = targ_tr_ls) # get the nodes to impute genus
            tps_2add <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(tip_names) # get tips to impute per genus
            
            # impute the first tip per genus to anchor the clade
            if(length(tps_2add) > 1) {  # in case of more than one tip
              targ_tr_ls <- TreeTools::AddTip(tree = targ_tr_ls, 
                                                   where = c(Ntip(targ_tr_ls) + nd_where), label = tps_2add[1]) %>% 
                castor::extend_tree_to_height() %>% .$tree
              for(l in 2:length(tps_2add)){ # then impute the other tips in the same genus
                tem_tr <- add.species.to.genus(tree = targ_tr_ls, 
                                               species = tps_2add[l], where = where2impt)
                if(is.ultrametric(tem_tr)) {targ_tr_ls <- tem_tr} else {targ_tr_ls <- force.ultrametric(tem_tr, method = "nnls", message = F)}
              }
            } else { # in case of one tip per genus
              targ_tr_ls <- TreeTools::AddTip(tree = targ_tr_ls, 
                                                   where = c(Ntip(targ_tr_ls) + nd_where), label = tps_2add) %>% 
                castor::extend_tree_to_height() %>% .$tree
            }
            impt_tps_vec <- c(impt_tps_vec, tps_2add)
          }
          
          tips2add_ls <- setdiff(tips2add_ls, impt_tps_vec) %>% sort
          for(l in 1:length(tips2add_ls)){
            tem_tr <- add.species.to.genus(tree = targ_tr_ls, 
                                           species = tips2add_ls[l], where = where2impt)
            if(is.ultrametric(tem_tr)) {targ_tr_ls <- tem_tr} else {targ_tr_ls <- force.ultrametric(tem_tr, method = "nnls", message = F)}
          }
          
        } else { # otherwise impute species to genus already in the tree
          for(l in 1:length(tips2add_ls)){
            tem_tr <- add.species.to.genus(tree = targ_tr_ls, 
                                           species = tips2add_ls[l], where = where2impt)
            if(is.ultrametric(tem_tr)) {targ_tr_ls <- tem_tr} else {targ_tr_ls <- force.ultrametric(tem_tr, method = "nnls", message = F)}
          }
        }
        ifelse(targ_tr_ls$tip.label %in% fern_tre$tip.label %>% sum == length(fern_tre$tip.label), 
               print(paste0("All tips imputed")), print(paste0("WARNING ! unequal number of tips")))
       
      print(paste0("## file: ", i))
      targ_tr_ls %>% 
        saveRDS(paste0("output/imputed_trees/late_spat_pr_tax_imp_", where2impt, "_", late_frc[i], ".rds"))
    }
  } else {
    for(i in 1:length(late_fls)){
      targ_tr_ls <- readRDS(paste0("output/imputed_trees/late_spat_pr_tax_imp_", where2impt, "_",  late_frc[i], ".rds"))
    }
  }
  
  # marginal pruned
  if(run_marg_tax){ 
    for(i in 1:length(marg_fls)){ 
      trees_ls <- readRDS(marg_fls[[i]]) # load set of trees
      for(h in 1:length(trees_ls)){trees_ls$node.label <- c()}
      trees_ls <- name_nodes(tree = trees_ls, m_tax = tax_tbl) # name nodes 
      
      targ_tr_ls <- trees_ls # set trees to impute
      tips2add_ls <- setdiff(fern_tre$tip.label, targ_tr_ls$tip.label) # get missing tips per tree
      
        gens2add <- strsplit(tips2add_ls, '_', fixed = F) %>% # get the genera from the tips to be added in the target tree
          sapply("[", 1) %>% unique() %>% sort()
        gens_tree <- strsplit(targ_tr_ls$tip.label, '_', fixed = F) %>% # get the genera in the target tree
          sapply("[", 1) %>% unique() %>% sort
        
        gen_diff <- setdiff(gens2add, gens_tree) %>% sort # get the missing genera in the target tree that need to be imputed first
        
        tps_gen_impt <- c() ; impt_tps_vec <- c()
        if(length(gen_diff) > 0) { # if the genus is not in the tree, it needs to be imputed first
          for(j in 1:length(gen_diff)){
            nd_where <- node2where(gen = gen_diff[j], tree = targ_tr_ls) # get the nodes to impute genus
            tps_2add <- tax_tbl %>% filter(genus == gen_diff[j]) %>% pull(tip_names) # get tips to impute per genus
            
            # impute the first tip per genus to anchor the clade
            if(length(tps_2add) > 1) {  # in case of more than one tip
              targ_tr_ls <- TreeTools::AddTip(tree = targ_tr_ls, 
                                                   where = c(Ntip(targ_tr_ls) + nd_where), label = tps_2add[1]) %>% 
                castor::extend_tree_to_height() %>% .$tree
              for(l in 2:length(tps_2add)){ # then impute the other tips in the same genus
                tem_tr <- add.species.to.genus(tree = targ_tr_ls, 
                                               species = tps_2add[l], where = where2impt)
                if(is.ultrametric(tem_tr)) {targ_tr_ls <- tem_tr} else {targ_tr_ls <- force.ultrametric(tem_tr, method = "nnls", message = F)}
              }
            } else { # in case of one tip per genus
              targ_tr_ls <- TreeTools::AddTip(tree = targ_tr_ls, 
                                                   where = c(Ntip(targ_tr_ls) + nd_where), label = tps_2add) %>% 
                castor::extend_tree_to_height() %>% .$tree
            }
            impt_tps_vec <- c(impt_tps_vec, tps_2add)
          }
          
          tips2add_ls <- setdiff(tips2add_ls, impt_tps_vec) %>% sort
          for(l in 1:length(tips2add_ls)){
            tem_tr <- add.species.to.genus(tree = targ_tr_ls, 
                                           species = tips2add_ls[l], where = where2impt)
            if(is.ultrametric(tem_tr)) {targ_tr_ls <- tem_tr} else {targ_tr_ls <- force.ultrametric(tem_tr, method = "nnls", message = F)}
          }
          
        } else { # otherwise impute species to genus already in the tree
          for(l in 1:length(tips2add_ls)){
            tem_tr <- add.species.to.genus(tree = targ_tr_ls, 
                                           species = tips2add_ls[l], where = where2impt)
            if(is.ultrametric(tem_tr)) {targ_tr_ls <- tem_tr} else {targ_tr_ls <- force.ultrametric(tem_tr, method = "nnls", message = F)}
            
          }
        }
        ifelse(targ_tr_ls$tip.label %in% fern_tre$tip.label %>% sum == length(fern_tre$tip.label), 
               print(paste0("All tips imputed")), print(paste0("WARNING ! unequal number of tips")))
       
      print(paste0("## file: ", i))
      targ_tr_ls %>% 
        saveRDS(paste0("output/imputed_trees/marg_spat_pr_tax_imp_", where2impt, "_", marg_frc[i], ".rds"))
    }
  } else {
    for(i in 1:length(marg_fls)){
      targ_tr_ls <- readRDS(paste0("output/imputed_trees/marg_spat_pr_tax_imp_", where2impt, "_",  marg_frc[i], ".rds"))
    }
  }
  
  print(paste0("END ", where2impt, " imputation"))
}





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

## only pruned ####
# lateral pruning 
late_fls <- list.files("output/pruned_trees/", "lateral_spat_pruned", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_fls)

if(run_late){
  for(i in 1:length(late_fls)){
    late_canaper_ls <- list()
    targ_tr_ls <- readRDS(late_fls[i])
    
    ifelse(Ntip(targ_tr_ls) %>% max() <= Ntip(fern_tre), "Ok n tips", stop(print("Error")))
    
      ifelse(sum(names(occ_mat) %in% targ_tr_ls$tip.label) < 30, next, "")
      
      occ_mat2 <- occ_mat[rowSums(occ_mat[, which(targ_tr_ls$tip.label %in% colnames(occ_mat))]) > 0, ]
      rand_test <- cpr_rand_test(occ_mat2, targ_tr_ls, 
                                 null_model = "curveball", 
                                 n_iterations = saturation_point,
                                 n_reps = canaper_n_reps,
                                 tbl_out = T)
      
      canape_res <- cpr_classify_endem(rand_test)
      sig_list <- list(canape_res, 
                       cpr_classify_signif(rand_test, metric = "pd")  %>% select(site, "pd_signif"), 
                       cpr_classify_signif(rand_test, metric = "pe")  %>% select(site, "pe_signif"), 
                       cpr_classify_signif(rand_test, metric = "rpd") %>% select(site, "rpd_signif"), 
                       cpr_classify_signif(rand_test, metric = "rpe") %>% select(site, "rpe_signif")
      )
      
      late_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i))

    late_canaper_ls %>% 
      saveRDS(paste0("output/canaper/canaper_late_spat_pr_", late_frc[i],".rds"))
    gc()
  }
} else {} 

# marginal pruning 
marg_fls <- list.files("output/pruned_trees/", "marginal_spat_pruned", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_fls)

if(run_marg){
  for(i in 1:length(marg_fls)){
    marg_canaper_ls <- list()
    targ_tr_ls <- readRDS(marg_fls[i])
    
      ifelse(sum(names(occ_mat) %in% targ_tr_ls$tip.label) < 30, next, "")
      
      occ_mat2 <- occ_mat[rowSums(occ_mat[, which(targ_tr_ls$tip.label %in% colnames(occ_mat))]) > 0, ]
      rand_test <- cpr_rand_test(occ_mat2, targ_tr_ls, 
                                 null_model = "curveball", 
                                 n_iterations = saturation_point,
                                 n_reps = canaper_n_reps,
                                 tbl_out = T)
      
      canape_res <- cpr_classify_endem(rand_test)
      sig_list <- list(canape_res, 
                       cpr_classify_signif(rand_test, metric = "pd")  %>% select(site, "pd_signif"), 
                       cpr_classify_signif(rand_test, metric = "pe")  %>% select(site, "pe_signif"), 
                       cpr_classify_signif(rand_test, metric = "rpd") %>% select(site, "rpd_signif"), 
                       cpr_classify_signif(rand_test, metric = "rpe") %>% select(site, "rpe_signif")
      )
      
      marg_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i))
    
    marg_canaper_ls %>% 
      saveRDS(paste0("output/canaper/canaper_marg_spat_pr_", marg_frc[i],".rds"))
    gc()
  }
} else {}

## pruned and randomly imputed ####
# lateral pruning
pttr <- "late_spat_pr_rand_imp_"
late_fls <- list.files("output/imputed_trees/", pttr, all.files = T, full.names = T)
late_frc <- readr::parse_number(late_fls)

if(run_late_rand){
  for(i in 1:length(late_fls)){
    late_canaper_ls <- list()
    targ_tr_ls <- readRDS(late_fls[i])
    
    ifelse(sapply(targ_tr_ls, Ntip) %>% max() <= Ntip(fern_tre), "Ok n tips", stop(print("Error")))
    
      ifelse(sum(names(occ_mat) %in% targ_tr_ls$tip.label) < 30, next, "")
      rand_test <- cpr_rand_test(occ_mat, targ_tr_ls, 
                                 null_model = "curveball", 
                                 n_iterations = saturation_point,
                                 n_reps = canaper_n_reps,
                                 tbl_out = T)
      
      canape_res <- cpr_classify_endem(rand_test)
      sig_list <- list(canape_res, 
                       cpr_classify_signif(rand_test, metric = "pd")  %>% select(site, "pd_signif"), 
                       cpr_classify_signif(rand_test, metric = "pe")  %>% select(site, "pe_signif")
                       #, cpr_classify_signif(rand_test, metric = "rpd") %>% select(site, "rpd_signif"), 
                       #cpr_classify_signif(rand_test, metric = "rpe") %>% select(site, "rpe_signif")
      )
      
      late_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i))

    late_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr, late_frc[i],".rds"))
    gc()
  }
} else {}

# marginal pruning
pttr <- "marg_spat_pr_rand_imp_"
marg_fls <- list.files("output/imputed_trees/", pttr, all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_fls)

if(run_marg_rand){
  for(i in 1:length(marg_fls)){
    marg_canaper_ls <- list()
    targ_tr_ls <- readRDS(marg_fls[i])

      ifelse(sum(names(occ_mat) %in% targ_tr_ls$tip.label) < 30, next, "")
      rand_test <- cpr_rand_test(occ_mat, targ_tr_ls, 
                                 null_model = "curveball", 
                                 n_iterations = saturation_point,
                                 n_reps = canaper_n_reps,
                                 tbl_out = T)
      
      canape_res <- cpr_classify_endem(rand_test)
      sig_list <- list(canape_res, 
                       cpr_classify_signif(rand_test, metric = "pd")  %>% select(site, "pd_signif"), 
                       cpr_classify_signif(rand_test, metric = "pe")  %>% select(site, "pe_signif")
                       #, cpr_classify_signif(rand_test, metric = "rpd") %>% select(site, "rpd_signif"), 
                       #cpr_classify_signif(rand_test, metric = "rpe") %>% select(site, "rpe_signif")
      )
      
      marg_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i))
    
    marg_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr, marg_frc[i],".rds"))
    gc()
  }
} else {}

## pruned and taxonomic imputed ####
# lateral pruning
pttr <- c("late_spat_pr_tax_imp_root", "late_spat_pr_tax_imp_random")

for(f in 1:length(pttr)){
  late_fls <- list.files("output/imputed_trees/", pttr[f], all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_fls)
  if(run_late_tax){
    for(i in 1:length(late_fls)){
      late_canaper_ls <- list()
      targ_tr_ls <- readRDS(late_fls[i])
      
      ifelse(sapply(targ_tr_ls, Ntip) %>% max() <= Ntip(fern_tre), "Ok n tips", stop(print("Error")))
      
        ifelse(sum(names(occ_mat) %in% targ_tr_ls$tip.label) < 30, next, "")
        rand_test <- cpr_rand_test(occ_mat, targ_tr_ls, 
                                   null_model = "curveball", 
                                   n_iterations = saturation_point,
                                   n_reps = canaper_n_reps,
                                   tbl_out = T)
        
        canape_res <- cpr_classify_endem(rand_test)
        sig_list <- list(canape_res, 
                         cpr_classify_signif(rand_test, metric = "pd")  %>% select(site, "pd_signif"), 
                         cpr_classify_signif(rand_test, metric = "pe")  %>% select(site, "pe_signif"), 
                         cpr_classify_signif(rand_test, metric = "rpd") %>% select(site, "rpd_signif"), 
                         cpr_classify_signif(rand_test, metric = "rpe") %>% select(site, "rpe_signif")
        )
        
        late_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
        print(paste0("iteration: ", i))
      
      late_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr[f], "_", late_frc[i],".rds"))
      gc()
    }
  } else {}
  cat(paste0("END: ", pttr[f]))
}

# marginal pruning
pttr <- c("marg_spat_pr_tax_imp_root", "marg_spat_pr_tax_imp_random")
for(f in 1:length(pttr)){
  marg_fls <- list.files("output/imputed_trees/", pttr[f], all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_fls)
  if(run_marg_tax){
    for(i in 1:length(marg_fls)){
      marg_canaper_ls <- list()
      targ_tr_ls <- readRDS(marg_fls[i])
      
        ifelse(sum(names(occ_mat) %in% targ_tr_ls$tip.label) < 30, next, "")
        rand_test <- cpr_rand_test(occ_mat, targ_tr_ls, 
                                   null_model = "curveball", 
                                   n_iterations = saturation_point,
                                   n_reps = canaper_n_reps,
                                   tbl_out = T)
        
        canape_res <- cpr_classify_endem(rand_test)
        sig_list <- list(canape_res, 
                         cpr_classify_signif(rand_test, metric = "pd")  %>% select(site, "pd_signif"), 
                         cpr_classify_signif(rand_test, metric = "pe")  %>% select(site, "pe_signif"), 
                         cpr_classify_signif(rand_test, metric = "rpd") %>% select(site, "rpd_signif"), 
                         cpr_classify_signif(rand_test, metric = "rpe") %>% select(site, "rpe_signif")
        )
        
        marg_canaper_ls <- plyr::join_all(sig_list, by = 'site', type = 'left')
        print(paste0("iteration: ", i))
      
      marg_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr[f], "_",  marg_frc[i],".rds"))
      gc()
    }
  } else {}
  cat(paste0("END: ", pttr[f]))
}
