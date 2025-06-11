library(tidyverse)

cols2sel <- c("site", 
              "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              "pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean")

# load pruned ####
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_[0-9]", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_[0-9]", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_[0-9]", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

# random
rnd_canape_df_ls <- list()
for(i in 1:length(rnd_canaper_fls)){
  rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
  rnd_canape_df_ls[[i]] <- rnd_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls))), 
        samp_frac = rand_frc[i], treat = "rand-pr", imputation = "none") %>% 
    bind_rows()
}
rnd_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_canape_all_df)

# lateral 
late_canape_df_ls <- list()
for(i in 1:length(late_canaper_fls)){
  late_canape_ls <- readRDS(late_canaper_fls[i])
  late_canape_df_ls[[i]] <- late_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls))), 
        samp_frac = late_frc[i], treat = "late-pr", imputation = "none") %>% 
    bind_rows()
}
late_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_canape_all_df)

# marginal
marg_canape_df_ls <- list()
for(i in 1:length(marg_canaper_fls)){
  marg_canape_ls <- readRDS(marg_canaper_fls[i])
  marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canape_ls))), 
        samp_frac = marg_frc[i], treat = "marg-pr", imputation = "none") %>% 
    bind_rows()
}
marg_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_canape_all_df)

prun_res_all <- rbind(rnd_canape_all_df, late_canape_all_df, marg_canape_all_df) ; dim(prun_res_all)

# load pruned  - random imputed ####
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_rand", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_rand", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_rand", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

# random
rnd_canape_df_ls <- list()
for(i in 1:length(rnd_canaper_fls)){
  rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
  rnd_canape_df_ls[[i]] <- rnd_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls))), 
        samp_frac = rand_frc[i], treat = "rand-pr", imputation = "random") %>% 
    bind_rows()
}
rnd_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_canape_all_df)

# lateral 
late_canape_df_fls <- list()
for(i in 1:length(late_canaper_fls)){
  late_canape_ls <- readRDS(late_canaper_fls[i])
  late_canape_df_ls[[i]] <- late_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls))), 
        samp_frac = late_frc[i], treat = "late-pr", imputation = "random") %>% 
    bind_rows()
}
late_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_canape_all_df)

# marginal
marg_canape_df_fls <- list()
for(i in 1:length(marg_canaper_fls)){
  marg_canape_ls <- readRDS(marg_canaper_fls[i])
  marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canape_ls))), 
        samp_frac = marg_frc[i], treat = "marg-pr", imputation = "random") %>% 
    bind_rows()
}
marg_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_canape_all_df)

rand_imp_res_all <- rbind(rnd_canape_all_df, late_canape_all_df, marg_canape_all_df) ; dim(rand_imp_res_all)

# load pruned - taxonomic imputed (node - clade) ####
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_tax_imp_root", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_tax_imp_root", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_tax_imp_root", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

# random
rnd_canape_df_ls <- list()
for(i in 1:length(rnd_canaper_fls)){
  rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
  rnd_canape_df_ls[[i]] <- rnd_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls))), 
        samp_frac = rand_frc[i], treat = "rand-pr", imputation = "taxonomy - node") %>% 
    bind_rows()
}
rnd_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_canape_all_df)

# lateral 
late_canape_df_fls <- list()
for(i in 1:length(late_canaper_fls)){
  late_canape_ls <- readRDS(late_canaper_fls[i])
  late_canape_df_ls[[i]] <- late_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls))), 
        samp_frac = late_frc[i], treat = "late-pr", imputation = "taxonomy - node") %>% 
    bind_rows()
}
late_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_canape_all_df)

# marginal
marg_canape_df_fls <- list()
for(i in 1:length(marg_canaper_fls)){
  marg_canape_ls <- readRDS(marg_canaper_fls[i])
  marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canape_ls))), 
        samp_frac = marg_frc[i], treat = "marg-pr", imputation = "taxonomy - node") %>% 
    bind_rows()
}
marg_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_canape_all_df)

tax_imp_node_res_all <- rbind(rnd_canape_all_df, late_canape_all_df, marg_canape_all_df) ; dim(tax_imp_node_res_all)

# load pruned - taxonomic imputed (within - random) ####
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_tax_imp_random", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_tax_imp_random", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_tax_imp_random", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

# random
rnd_canape_df_ls <- list()
for(i in 1:length(rnd_canaper_fls)){
  rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
  rnd_canape_df_ls[[i]] <- rnd_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls))), 
        samp_frac = rand_frc[i], treat = "rand-pr", imputation = "taxonomy - random") %>% 
    bind_rows()
}
rnd_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_canape_all_df)

# lateral 
late_canape_df_fls <- list()
for(i in 1:length(late_canaper_fls)){
  late_canape_ls <- readRDS(late_canaper_fls[i])
  late_canape_df_ls[[i]] <- late_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls))), 
        samp_frac = late_frc[i], treat = "late-pr", imputation = "taxonomy - random") %>% 
    bind_rows()
}
late_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_canape_all_df)

# marginal
marg_canape_df_fls <- list()
for(i in 1:length(marg_canaper_fls)){
  marg_canape_ls <- readRDS(marg_canaper_fls[i])
  marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canape_ls))), 
        samp_frac = marg_frc[i], treat = "marg-pr", imputation = "taxonomy - random") %>% 
    bind_rows()
}
marg_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_canape_all_df)

tax_imp_rand_res_all <- rbind(rnd_canape_all_df, late_canape_all_df, marg_canape_all_df) ; dim(tax_imp_rand_res_all)

# all data  results raw #####
all_canape_res <- rbind(prun_res_all, rand_imp_res_all, tax_imp_node_res_all, tax_imp_rand_res_all)

saveRDS(all_canape_res, "output/canaper_all_data_res.rds")
all_canape_res <- readRDS("output/canaper_all_data_res.rds")
