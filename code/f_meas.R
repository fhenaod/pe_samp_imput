library(tidyverse)

cols2sel <- c("site", "endem_type")

# load results using Nitta et al. 2022 dataset ####
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))
emp_df <- emp_canape_res %>% select(cols2sel) ; emp_df %>% head()

spat_df <- read.csv("data/nitta-et-al/results/japan_ferns_biodiv.csv", header = T) %>% 
  select(grids, lat, long) %>% mutate(site = as.character(grids))

# load results data - analyses and make tables ####
## pruned ####
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_[0-9]", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

get_tabs2res <- F
if(get_tabs2res){
  # random
  rnd_canape_df_ls <- list()
  for(i in 1:length(rnd_canaper_fls)){
    rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
    rnd_canape_df_ls[[i]] <- rnd_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls))), 
          samp_frac = rand_frc[i], treat = "rand-pr") %>% 
      bind_rows()
  }
  rnd_pr_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_pr_canape_all_df)
  
  # quatitative trait
  late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_[0-9]", all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_canaper_fls)
  
  marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_[0-9]", all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_canaper_fls)
  
  # lateral - quantitative trait
  late_canape_df_ls <- list()
  for(i in 1:length(late_canaper_fls)){
    late_canape_ls <- readRDS(late_canaper_fls[i])
    late_canape_df_ls[[i]] <- late_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canaper_fls))), 
          samp_frac = late_frc[i], treat = "late-pr") %>% 
      bind_rows()
  }
  late_pr_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_pr_canape_all_df)
  
  # marginal - quantitative trait
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_canaper_fls)){
    marg_canape_ls <- readRDS(marg_canaper_fls[i])
    marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canaper_fls))), 
          samp_frac = marg_frc[i], treat = "marg-pr") %>% 
      bind_rows()
  }
  marg_pr_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_pr_canape_all_df)
  
  # spatial
  late_spat_canaper_fls <- list.files("output/canaper", "canaper_late_spat_pr_[0-9]", all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_spat_canaper_fls)
  
  marg_spat_canaper_fls <- list.files("output/canaper", "canaper_marg_spat_pr_[0-9]", all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_spat_canaper_fls)
  
  # lateral - spatial
  late_canape_df_ls <- list()
  for(i in 1:length(late_spat_canaper_fls)){
    late_canape_ls <- readRDS(late_spat_canaper_fls[i])
    late_canape_df_ls[[i]] <- late_canape_ls %>% 
      select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = late_frc[i], treat = "late-spat_pr")
  }
  late_spat_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_spat_canape_all_df) %>% print()
  
  # marginal - spatial
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_spat_canaper_fls)){
    marg_canape_ls <- readRDS(marg_spat_canaper_fls[i])
    marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
      select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = marg_frc[i], treat = "marg-spat_pr")
  }
  marg_spat_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_spat_canape_all_df) %>% print()
  
  prun_res_all <- rbind(rnd_pr_canape_all_df, 
                        late_pr_canape_all_df, marg_pr_canape_all_df,
                        late_spat_canape_all_df, marg_spat_canape_all_df) ; head(prun_res_all)
  prun_res_all %>% saveRDS("output/prun_res_canaper.rds")
} else {
  prun_res_all <- readRDS("output/prun_res_canaper.rds") ; head(prun_res_all)
}

## random imputed ####
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_rand", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

get_tabs2res <- F
if(get_tabs2res){
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
  rnd_pr_rnd_imp_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_pr_rnd_imp_canape_all_df) %>% print()
  
  # quantitative trait
  late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_rand", all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_canaper_fls)
  
  marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_rand", all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_canaper_fls)
  
  # lateral - quatitative traitr
  late_canape_df_ls <- list()
  for(i in 1:length(late_canaper_fls)){
    late_canape_ls <- readRDS(late_canaper_fls[i])
    late_canape_df_ls[[i]] <- late_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canaper_fls))), 
          samp_frac = late_frc[i], treat = "late-pr", imputation = "random") %>% 
      bind_rows()
  }
  late_pr_rnd_imp_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_pr_rnd_imp_canape_all_df) %>% print()
  
  # marginal - quatitative traitr
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_canaper_fls)){
    marg_canape_ls <- readRDS(marg_canaper_fls[i])
    marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canaper_fls))), 
          samp_frac = marg_frc[i], treat = "marg-pr", imputation = "random") %>% 
      bind_rows()
  }
  marg_pr_rand_imp_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_pr_rand_imp_canape_all_df) %>% print()
  
  # spatial pruning ###
  late_spat_canaper_fls <- list.files("output/canaper", "canaper_late_spat_pr_rand", all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_canaper_fls)
  
  marg_spat_canaper_fls <- list.files("output/canaper", "canaper_marg_spat_pr_rand", all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_canaper_fls)
  
  # lateral - spatial
  late_canape_df_ls <- list()
  for(i in 1:length(late_spat_canaper_fls)){
    late_canape_ls <- readRDS(late_spat_canaper_fls[i])
    late_canape_df_ls[[i]] <- late_canape_ls %>% 
      select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = late_frc[i], treat = "late-spat_pr", imputation = "random")
  }
  late_spat_pr_rand_imp_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_spat_pr_rand_imp_canape_all_df) %>% print()
  
  # marginal - spatial
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_spat_canaper_fls)){
    marg_canape_ls <- readRDS(marg_spat_canaper_fls[i])
    marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
      select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = marg_frc[i], treat = "marg-spat_pr", imputation = "random")
  }
  marg_spat_pr_rand_imp_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_spat_pr_rand_imp_canape_all_df) %>% print()
  
  rand_imp_all <- rbind(rnd_pr_rnd_imp_canape_all_df, 
                        late_pr_rnd_imp_canape_all_df, marg_pr_rand_imp_canape_all_df,
                        late_spat_pr_rand_imp_canape_all_df, marg_spat_pr_rand_imp_canape_all_df)
  rand_imp_all %>% saveRDS("output/prun_rand_imp_res_canaper.rds")
} else {
  rand_imp_all <- readRDS("output/prun_rand_imp_res_canaper.rds"); head(rand_imp_all)
}

## taxonomic imputed #####
# node
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_tax_imp_root", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

get_tabs2res <- F
if(get_tabs2res){
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
  rnd_pr_tax_imp_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_pr_tax_imp_canape_all_df) %>% print()
  
  late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_tax_imp_root", all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_canaper_fls)
  
  marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_tax_imp_root", all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_canaper_fls)
  
  # quantitative trait pruning
  # lateral - quantitative trait
  late_canape_df_ls <- list()
  for(i in 1:length(late_canaper_fls)){
    late_canape_ls <- readRDS(late_canaper_fls[i])
    late_canape_df_ls[[i]] <- late_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canaper_fls))), 
          samp_frac = late_frc[i], treat = "late-pr", imputation = "taxonomy - node") %>% 
      bind_rows()
  }
  late_pr_tax_imp_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_pr_tax_imp_canape_all_df) %>% print()
  
  # marginal - quantitative trait
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_canaper_fls)){
    marg_canape_ls <- readRDS(marg_canaper_fls[i])
    marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canaper_fls))), 
          samp_frac = marg_frc[i], treat = "marg-pr", imputation = "taxonomy - node") %>% 
      bind_rows()
  }
  marg_pr_tax_imp_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_pr_tax_imp_canape_all_df) %>% print()
  
  # spatial pruning
  late_spat_canaper_fls <- list.files("output/canaper", "canaper_late_spat_pr_tax_imp_root", all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_canaper_fls)
  
  marg_spat_canaper_fls <- list.files("output/canaper", "canaper_marg_spat_pr_tax_imp_root", all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_canaper_fls)
  
  # lateral - spatial
  late_canape_df_ls <- list()
  for(i in 1:length(late_spat_canaper_fls)){
    late_canape_ls <- readRDS(late_spat_canaper_fls[i])
    late_canape_df_ls[[i]] <- late_canape_ls %>% 
      select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = late_frc[i], treat = "late-spat_pr", imputation = "taxonomy - node")
  }
  late_spat_pr_tax_imp_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_spat_pr_tax_imp_canape_all_df) %>% print()
  
  # marginal - spatial
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_spat_canaper_fls)){
    marg_canape_ls <- readRDS(marg_spat_canaper_fls[i])
    marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
      select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = marg_frc[i], treat = "marg-spat_pr", imputation = "taxonomy - node")
  }
  marg_spat_pr_tax_imp_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_spat_pr_tax_imp_canape_all_df) %>% print()
  
  
  prun_tax_imp_node_all <- rbind(rnd_pr_tax_imp_canape_all_df, 
                                 late_pr_tax_imp_canape_all_df, marg_pr_tax_imp_canape_all_df,
                                 late_spat_pr_tax_imp_canape_all_df, marg_spat_pr_tax_imp_canape_all_df)
  prun_tax_imp_node_all %>% saveRDS("output/prun_tax_imp-node_res_canaper.rds")
} else {
  prun_tax_imp_node_all <- readRDS("output/prun_tax_imp-node_res_canaper.rds") ; head(prun_tax_imp_node_all)
}

## random within
# node
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_tax_imp_random", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

get_tabs2res <- F
if(get_tabs2res){
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
  rnd_pr_tax_imp_rand_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_pr_tax_imp_rand_canape_all_df) %>% print()
  
  late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_tax_imp_random", all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_canaper_fls)
  
  marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_tax_imp_random", all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_canaper_fls)
  
  # quantitative trait pruning 
  # lateral - quantitative trait
  late_canape_df_ls <- list()
  for(i in 1:length(late_canaper_fls)){
    late_canape_ls <- readRDS(late_canaper_fls[i])
    late_canape_df_ls[[i]] <- late_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canaper_fls))), 
          samp_frac = late_frc[i], treat = "late-pr", imputation = "taxonomy - random") %>% 
      bind_rows()
  }
  late_pr_tax_imp_rand_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_pr_tax_imp_rand_canape_all_df) %>% print()
  
  # marginal - quantitative trait
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_canaper_fls)){
    marg_canape_ls <- readRDS(marg_canaper_fls[i])
    marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canaper_fls))), 
          samp_frac = marg_frc[i], treat = "marg-pr", imputation = "taxonomy - random") %>% 
      bind_rows()
  }
  marg_pr_tax_imp_rand_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_pr_tax_imp_rand_canape_all_df) %>% print()
  
  # spatial pruning ##
  late_spat_canaper_fls <- list.files("output/canaper", "canaper_late_spat_pr_tax_imp_random", all.files = T, full.names = T)
  late_frc <- readr::parse_number(late_canaper_fls)
  
  marg_spat_canaper_fls <- list.files("output/canaper", "canaper_marg_spat_pr_tax_imp_random", all.files = T, full.names = T)
  marg_frc <- readr::parse_number(marg_canaper_fls)
  
  # lateral - spatial
  late_canape_df_ls <- list()
  for(i in 1:length(late_spat_canaper_fls)){
    late_canape_ls <- readRDS(late_spat_canaper_fls[i])
    late_canape_df_ls[[i]] <- late_canape_ls %>% 
      select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = late_frc[i], treat = "late-spat_pr", imputation = "taxonomy - random")
  }
  late_spat_pr_tax_imp_rand_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_spat_pr_tax_imp_rand_canape_all_df) %>% print()
  
  # marginal - spatial
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_spat_canaper_fls)){
    marg_canape_ls <- readRDS(marg_spat_canaper_fls[i])
    marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
      select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = marg_frc[i], treat = "marg-spat_pr", imputation = "taxonomy - random")
  }
  marg_spat_pr_tax_imp_rand_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_spat_pr_tax_imp_rand_canape_all_df) %>% print()
  
  
  prun_tax_imp_rand_all <- rbind(rnd_pr_tax_imp_rand_canape_all_df, 
                                 late_pr_tax_imp_rand_canape_all_df, marg_pr_tax_imp_rand_canape_all_df,
                                 late_spat_pr_tax_imp_rand_canape_all_df, marg_spat_pr_tax_imp_rand_canape_all_df)
  prun_tax_imp_rand_all %>% saveRDS("output/prun_tax_imp-rand_res_canaper.rds")
} else {
  prun_tax_imp_rand_all <- readRDS("output/prun_tax_imp-rand_res_canaper.rds") #; head(prun_tax_imp_rand_all)
}

# make one big df
canaper_res_df <- rbind(prun_res_all %>% mutate(imputation = "none"), 
                        rand_imp_all, prun_tax_imp_node_all, prun_tax_imp_rand_all)

# join with empirical df
canaper_res_all_df <- left_join(canaper_res_df, 
                                emp_df %>% rename(emp_endem_type = endem_type), 
                                by = join_by(site))

# Estimate: Precision, recall and f-meas #####
library(yardstick)

canaper_stat <- canaper_res_all_df %>% 
  #filter(endem_type !=is.na(endem_type)) %>% 
  mutate(emp_endem_type = factor(emp_endem_type), endem_type = factor(endem_type)) %>%
  group_by(sim_n, samp_frac, treat, imputation) %>% 
  summarize(precision = precision_vec(truth = emp_endem_type, estimate = endem_type, na_rm = T),
            recall = recall_vec(truth = emp_endem_type, estimate = endem_type, na_rm = T),
            f_mes = f_meas_vec(truth = emp_endem_type, estimate = endem_type, na_rm = T))

canaper_mean_stats <- canaper_stat %>% 
  ungroup() %>% group_by(samp_frac, treat, imputation) %>% 
  summarise(mean_prec = mean(precision, na.rm = T),
            sd_prec = sd(precision, na.rm = T),
            mean_reca = mean(recall, na.rm = T),
            sd_reca = sd(recall, na.rm = T),
            mean_fmes = mean(f_mes, na.rm = T),
            sd_fmes = sd(f_mes, na.rm = T))

# figs #####
forms2print <- c(".png", ".pdf", ".svg") 

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  # F measure
  canaper_mean_stats %>% 
    mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
           cat2map = paste0(treat, "-" , imputation)) %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, color = imputation, shape = treat)) +
    geom_point() + 
    scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
    geom_line(aes(linetype = imputation)) + 
    labs(title = "Endemism type", 
         x = "Proportion of tips removed", y = "Mean F-score", 
         color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
    theme_classic(base_size = 12) + 
    ylim(c(0.35, .9)) +
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    theme(strip.background = element_blank(), legend.position = "right") + 
    facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2, scales = "fixed")
  ggsave(filename = paste0("figures/y_meas_canaper", fig_out), 
         width = 25, height = 15, units = "cm", dpi = 200)

  # with sd
  canaper_mean_stats %>% 
    mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
           cat2map = paste0(treat, "-" , imputation)) %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
               color = imputation, shape = treat)) +
    geom_pointrange(position = position_dodge2(width = .025, padding = 0.1)) + 
    scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
    geom_line(aes(linetype = imputation)) + 
    #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
    labs(#title = "Endemism categorization", 
         x = "Proportion of tips removed", y = "Mean F-score", 
         color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
    theme_classic(base_size = 12) + 
    ylim(c(0.35, .9)) +
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
           shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
    theme(strip.background = element_blank(), legend.position = "right") #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2, scales = "free_y")
  ggsave(filename = paste0("figures/y_meas_sd_canaper", fig_out), 
         width = 30, height = 15, units = "cm", dpi = 300)
  
  # precision
  canaper_mean_stats %>% 
    mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
           cat2map = paste0(treat, "-" , imputation)) %>%
    ggplot(aes(x = samp_frac, y = mean_prec, ymin = mean_prec - sd_prec, ymax = mean_prec + sd_prec,
               color = imputation, shape = treat)) +
    geom_pointrange(position = position_dodge2(width = .025, padding = 0.1)) + 
    scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
    geom_line(aes(linetype = imputation)) + 
    #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
    labs(title = "Endemism type", 
         x = "Proportion of tips removed", y = "Mean Precision", 
         color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
    theme_classic(base_size = 12) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    theme(strip.background = element_blank(), legend.position = "right") + 
    facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
  ggsave(filename = paste0("figures/prec_sd_canaper", fig_out), 
         width = 25, height = 15, units = "cm", dpi = 200)
  
  # recall
  canaper_mean_stats %>% 
    mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
           cat2map = paste0(treat, "-" , imputation)) %>%
    ggplot(aes(x = samp_frac, y = mean_reca, ymin = mean_reca - sd_reca, ymax = mean_reca + sd_reca,
               color = imputation, shape = treat)) +
    geom_pointrange(position = position_dodge2(width = .025, padding = 0.1)) + 
    scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
    geom_line(aes(linetype = imputation)) + 
    #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
    labs(title = "Endemism type", 
         x = "Proportion of tips removed", y = "Mean Recall", 
         color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
    theme_classic(base_size = 12) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    theme(strip.background = element_blank(), legend.position = "right") + 
    facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
  ggsave(filename = paste0("figures/reca_sd_canaper", fig_out), 
         width = 25, height = 15, units = "cm", dpi = 200)
}


canaper_mean_stats %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  filter(imputation != "random") %>% 
  ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .01, padding = 0.1)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(#title = "Endemism categorization", 
    x = "Proportion of tips removed", y = "Mean F-score", 
    color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  ylim(c(0.35, .9)) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3), shape = c(NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  theme(strip.background = element_blank(), legend.position = "right") + 
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2, scales = "free_y")
ggsave(filename = paste0("figures/y_meas_sd_canaper_2pnl", ".pdf"), 
       width = 30, height = 15, units = "cm", dpi = 300)

## animated #####
library(gganimate)

# with sd
f_meas_pt <- canaper_mean_stats %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  filter(imputation != "random") %>% 
  ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .01, padding = 0.1)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(#title = "Endemism categorization", 
    x = "Proportion of tips removed", y = "Mean F-score", 
    color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  ylim(c(0.35, .9)) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  guides(shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8))), 
         colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3), shape = c(NA, NA, NA)))) +
  theme(strip.background = element_blank(), legend.position = "right") + 
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2, scales = "free_y") +
  transition_reveal(samp_frac, keep_last = T)
animate(f_meas_pt, width = 30, height = 15, units = "cm", 
        res = 200, nframes = 60,  fps = 10, duration = 10, 
        start_pause = 2, end_pause = 60)
anim_save(filename = "figures/f_meas_anim.gif")


f_meas_pt2 <- canaper_mean_stats %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% head()
  filter(imputation != "random") %>% 
  ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .01, padding = 0.1)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(#title = "Endemism categorization", 
    x = "Proportion of tips removed", y = "Mean F-score", 
    color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  ylim(c(0.35, .9)) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  guides(shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8))), 
         colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3), shape = c(NA, NA, NA)))) +
  theme(strip.background = element_blank(), legend.position = "right") + 
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2, scales = "free_y") +
  transition_states(cat2map, transition_length = 2, state_length = 1, wrap = F) +
  ease_aes('sine-in-out')
animate(f_meas_pt2, width = 30, height = 15, units = "cm", 
        res = 200, nframes = 60,  fps = 10, duration = 10, 
        start_pause = 2, end_pause = 60)
anim_save(filename = "figures/f_meas_anim2.gif")

