library(tidyverse) ; library(ape)
library(canaper) ; library(future)

# load data ####
# empirical dataset
fern_tre <- read.tree("output/nitta_node-named.tre")
is.binary(fern_tre) ; is.ultrametric(fern_tre)

fern_com <- read.csv("data/nitta-et-al/japan_ferns_comm_full.csv", header = T) %>% 
  column_to_rownames("grids")

colnames(fern_com) <- gsub(pattern = "\\.", replacement = "-", colnames(fern_com))
colnames(fern_com) <- sub(pattern = "(\\_.*?)\\_", "\\1-", colnames(fern_com))

# missing taxa in tree
fern_com %>% colnames() %in% fern_tre$tip.label %>% table()

# canaper control pars ####
dir.create("output/canaper/")

run_iter <- F
run_emp <- F ; run_rand <- F ; run_late <- F ; run_marg <- F # pruned trees
run_rand_rand <- F ; run_late_rand <- F ; run_marg_rand <- F # pruned and randomly imputed trees
run_rand_tax <- T  ; run_late_tax <- T  ; run_marg_tax <- T  # pruned and taxonomically imputed trees

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

# empirical dataset ####
if(run_emp){
  rand_test <- cpr_rand_test(occ_mat, targ_tr, 
                                 null_model = "curveball", 
                                 n_iterations = saturation_point,
                                 n_reps = canaper_n_reps,
                                 tbl_out = T)

  rand_test[, 1:9] %>% head()

emp_canape_res <- cpr_classify_endem(rand_test)
sig_list <- list(emp_canape_res, 
                 cpr_classify_signif(rand_test, metric = "pd")  %>% select(site, "pd_signif"), 
                 cpr_classify_signif(rand_test, metric = "pe")  %>% select(site, "pe_signif"), 
                 cpr_classify_signif(rand_test, metric = "rpd") %>% select(site, "rpd_signif"), 
                 cpr_classify_signif(rand_test, metric = "rpe") %>% select(site, "rpe_signif")
)

plyr::join_all(sig_list, by = 'site', type = 'left') %>%
  saveRDS(paste0("output/canaper/canaper_empirical.rds"))
} else {}

# pruned dataset ####
# random pruning
rand_fls <- list.files("output/pruned_trees/", "rand_pruned", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rand_fls)

if(run_rand){
  for(i in 1:length(rand_fls)){
    rand_canaper_ls <- list()
    targ_tr_ls <- readRDS(rand_fls[i])
    ifelse(sapply(targ_tr_ls, Ntip) %>% max() <= Ntip(fern_tre), "Ok n tips", stop(print("Error")))
    
    for(j in 1:length(targ_tr_ls)){
      occ_mat2 <- occ_mat[rowSums(occ_mat[, which(targ_tr_ls[[j]]$tip.label %in% colnames(occ_mat))]) > 0, ]
      rand_test <- cpr_rand_test(occ_mat2, targ_tr_ls[[j]], 
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
      
      rand_canaper_ls[[j]] <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i, "-",  j))
    }
    rand_canaper_ls %>% 
      saveRDS(paste0("output/canaper/canaper_rand_pr_", rand_frc[i],".rds"))
    gc()
  }
} else {}

# lateral pruning
late_fls <- list.files("output/pruned_trees/", "lateral_pruned", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_fls)

if(run_late){
  for(i in 1:length(late_fls)){
    late_canaper_ls <- list()
    targ_tr_ls <- readRDS(late_fls[i])
    
    ifelse(sapply(targ_tr_ls, Ntip) %>% max() <= Ntip(fern_tre), "Ok n tips", stop(print("Error")))
    for(j in 1:length(targ_tr_ls)){
      ifelse(sum(names(occ_mat) %in% targ_tr_ls[[j]]$tip.label) < 30, next, "")
      
      occ_mat2 <- occ_mat[rowSums(occ_mat[, which(targ_tr_ls[[j]]$tip.label %in% colnames(occ_mat))]) > 0, ]
      rand_test <- cpr_rand_test(occ_mat2, targ_tr_ls[[j]], 
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
      
      late_canaper_ls[[j]] <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i, "-",  j))
    }
    late_canaper_ls %>% 
      saveRDS(paste0("output/canaper/canaper_late_pr_", late_frc[i],".rds"))
    gc()
  }
} else {} 

# marginal pruning
marg_fls <- list.files("output/pruned_trees/", "marginal_pruned", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_fls)

if(run_marg){
  for(i in 1:length(marg_fls)){
    marg_canaper_ls <- list()
    targ_tr_ls <- readRDS(marg_fls[i])
    for(j in 1:length(targ_tr_ls)){
      ifelse(sum(names(occ_mat) %in% targ_tr_ls[[j]]$tip.label) < 30, next, "")
      
      occ_mat2 <- occ_mat[rowSums(occ_mat[, which(targ_tr_ls[[j]]$tip.label %in% colnames(occ_mat))]) > 0, ]
      rand_test <- cpr_rand_test(occ_mat2, targ_tr_ls[[j]], 
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
      
      marg_canaper_ls[[j]] <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i, "-",  j))
    }
    marg_canaper_ls %>% 
      saveRDS(paste0("output/canaper/canaper_marg_pr_", marg_frc[i],".rds"))
    gc()
  }
} else {}

# imputed dataset ####
## randomly imputed tree tips #####
# random pruning
pttr <- "rand_pr_rand_imp_"
rand_fls <- list.files("output/imputed_trees/", pttr, all.files = T, full.names = T)
rand_frc <- readr::parse_number(rand_fls)

if(run_rand_rand){
  for(i in 1:length(rand_fls)){
    rand_canaper_ls <- list()
    targ_tr_ls <- readRDS(rand_fls[i])
    ifelse(sapply(targ_tr_ls, Ntip) %>% max() <= Ntip(fern_tre), "Ok n tips", stop(print("Error")))
    for(j in 1:length(targ_tr_ls)){
      rand_test <- cpr_rand_test(occ_mat, targ_tr_ls[[j]], 
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
      rand_canaper_ls[[j]] <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i, "-",  j))
    }
    rand_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr, rand_frc[i],".rds"))
    gc()
  }
} else {}

# lateral pruning
pttr <- "late_pr_rand_imp_"
late_fls <- list.files("output/imputed_trees/", pttr, all.files = T, full.names = T)
late_frc <- readr::parse_number(late_fls)

if(run_late_rand){
  for(i in 1:length(late_fls)){
    late_canaper_ls <- list()
    targ_tr_ls <- readRDS(late_fls[i])
    
    ifelse(sapply(targ_tr_ls, Ntip) %>% max() <= Ntip(fern_tre), "Ok n tips", stop(print("Error")))
    for(j in 1:length(targ_tr_ls)){
      ifelse(sum(names(occ_mat) %in% targ_tr_ls[[j]]$tip.label) < 30, next, "")
      rand_test <- cpr_rand_test(occ_mat, targ_tr_ls[[j]], 
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
      
      late_canaper_ls[[j]] <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i, "-",  j))
    }
    late_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr, late_frc[i],".rds"))
    gc()
  }
} else {}

# marginal pruning
pttr <- "marg_pr_rand_imp_"
marg_fls <- list.files("output/imputed_trees/", pttr, all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_fls)

if(run_marg_rand){
  for(i in 1:length(marg_fls)){
    marg_canaper_ls <- list()
    targ_tr_ls <- readRDS(marg_fls[i])
    for(j in 1:length(targ_tr_ls)){
      ifelse(sum(names(occ_mat) %in% targ_tr_ls[[j]]$tip.label) < 30, next, "")
      rand_test <- cpr_rand_test(occ_mat, targ_tr_ls[[j]], 
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
      
      marg_canaper_ls[[j]] <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i, "-",  j))
    }
    marg_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr, marg_frc[i],".rds"))
    gc()
  }
} else {}

## taxonomic imputed tree tips #####
# random pruning
pttr <- c("rand_pr_tax_imp_root", "rand_pr_tax_imp_random")
for(f in 1:length(pttr)){
  rand_fls <- list.files("output/imputed_trees", pttr[f], all.files = T, full.names = T)
  rand_frc <- readr::parse_number(rand_fls)
  
  if(run_rand_tax){
    for(i in 1:length(rand_fls)){
    rand_canaper_ls <- list()
    targ_tr_ls <- readRDS(rand_fls[i])
    ifelse(sapply(targ_tr_ls, Ntip) %>% max() <= Ntip(fern_tre), "Ok n tips", stop(print("Error")))
    for(j in 1:length(targ_tr_ls)){
      rand_test <- cpr_rand_test(occ_mat, targ_tr_ls[[j]], 
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
      
      rand_canaper_ls[[j]] <- plyr::join_all(sig_list, by = 'site', type = 'left')
      print(paste0("iteration: ", i, "-",  j))
    }
    rand_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr[f], "_", rand_frc[i],".rds"))
    gc()
    }
  } else {}
  cat(paste0("END: ", pttr[f]))
}

# lateral pruning
pttr <- c("late_pr_tax_imp_root", "late_pr_tax_imp_random")

for(f in 1:length(pttr)){
  late_fls <- list.files("output/imputed_trees/", pttr[f], all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_fls)
  if(run_late_tax){
    for(i in 1:length(late_fls)){
      late_canaper_ls <- list()
      targ_tr_ls <- readRDS(late_fls[i])
      
      ifelse(sapply(targ_tr_ls, Ntip) %>% max() <= Ntip(fern_tre), "Ok n tips", stop(print("Error")))
      for(j in 1:length(targ_tr_ls)){
        ifelse(sum(names(occ_mat) %in% targ_tr_ls[[j]]$tip.label) < 30, next, "")
        rand_test <- cpr_rand_test(occ_mat, targ_tr_ls[[j]], 
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
        
        late_canaper_ls[[j]] <- plyr::join_all(sig_list, by = 'site', type = 'left')
        print(paste0("iteration: ", i, "-",  j))
      }
      late_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr[f], "_", late_frc[i],".rds"))
      gc()
    }
  } else {}
  cat(paste0("END: ", pttr[f]))
}

# marginal pruning
pttr <- c("marg_pr_tax_imp_root", "marg_pr_tax_imp_random")
for(f in 1:length(pttr)){
  marg_fls <- list.files("output/imputed_trees/", pttr[f], all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_fls)
  if(run_marg_tax){
    for(i in 1:length(marg_fls)){
      marg_canaper_ls <- list()
      targ_tr_ls <- readRDS(marg_fls[i])
      for(j in 1:length(targ_tr_ls)){
        ifelse(sum(names(occ_mat) %in% targ_tr_ls[[j]]$tip.label) < 30, next, "")
        rand_test <- cpr_rand_test(occ_mat, targ_tr_ls[[j]], 
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
        
        marg_canaper_ls[[j]] <- plyr::join_all(sig_list, by = 'site', type = 'left')
        print(paste0("iteration: ", i, "-",  j))
      }
      marg_canaper_ls %>% saveRDS(paste0("output/canaper/canaper_", pttr[f], "_",  marg_frc[i],".rds"))
      gc()
    }
  } else {}
  cat(paste0("END: ", pttr[f]))
}
