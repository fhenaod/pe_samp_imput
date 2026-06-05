library(tidyverse)
library(yardstick)

# load and get canaper tables ####
# empirical full data
cols2sel <- c("site", 
              "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              "pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_obs_z", "pe_obs_z", "rpd_obs_z", "rpe_obs_z" #, "endem_type" 
              )

canaper_emp_raw <- readRDS("output/canaper/canaper_empirical.rds")
canaper_emp <- canaper_emp_raw %>% 
  dplyr::select(any_of(cols2sel)) %>% 
  mutate(treat = "empirical", samp_frac = 1, imputation = "none")

# results dataframe - analyses
out_fold <- "fern"
canaper_res_all <- readRDS(paste0("output/canaper/canaper_all_data_res_", out_fold,".rds"))

# mean metric across sims #### 
# get mean across 100 sims per samp fraction, pruning and imputation methods
canaper_res_df <- canaper_res_all %>% 
  group_by(site, samp_frac, treat, imputation) %>% 
  summarise(across(pd_obs:rpe_obs_z, ~mean(.x, na.rm = T), .names = "m_{.col}")) #%>% dplyr::select(names(canaper_emp))

# build df 
avg_df <- left_join(canaper_res_df %>% ungroup(), 
                    canaper_emp %>% ungroup() %>% 
                      dplyr::select(-c(treat:imputation)) %>% 
                      rename_if(is.numeric, ~paste0("emp_", .x)), by = join_by("site")) 

m_stats_df <- avg_df %>% 
  #filter(pd_obs !=is.na(pd_obs)) %>% 
  group_by(samp_frac, treat, imputation) %>% 
  summarise(mpe_pd  = mpe_vec(truth = emp_pd_obs, estimate = m_pd_obs, na_rm = T), 
            mpe_pe  = mpe_vec(truth = emp_pe_obs, estimate = m_pe_obs, na_rm = T), 
            mpe_rpd  = mpe_vec(truth = emp_rpd_obs, estimate = m_rpd_obs, na_rm = T), 
            mpe_rpe  = mpe_vec(truth = emp_rpe_obs, estimate = m_rpe_obs, na_rm = T), 
            mpe_pd_rand_mean  = mpe_vec(truth = emp_pd_rand_mean, estimate = m_pd_rand_mean, na_rm = T), 
            mpe_pe_rand_mean  = mpe_vec(truth = emp_pe_rand_mean, estimate = m_pe_rand_mean, na_rm = T), 
            mpe_rpd_rand_mean  = mpe_vec(truth = emp_rpd_rand_mean, estimate = m_rpd_rand_mean, na_rm = T), 
            mpe_rpe_rand_mean  = mpe_vec(truth = emp_rpe_rand_mean, estimate = m_rpe_rand_mean, na_rm = T),
            mpe_pd_obs_z = mpe_vec(truth = emp_pd_obs_z, estimate = m_pd_obs_z, na_rm = T),
            mpe_pe_obs_z = mpe_vec(truth = emp_pe_obs_z, estimate = m_pe_obs_z, na_rm = T),
            mpe_rpd_obs_z = mpe_vec(truth = emp_rpd_obs_z, estimate = m_rpd_obs_z, na_rm = T),
            mpe_rpe_obs_z = mpe_vec(truth = emp_rpe_obs_z, estimate = m_rpe_obs_z, na_rm = T),
            
            rmse_pd  = rmse_vec(truth = emp_pd_obs, estimate = m_pd_obs, na_rm = T), 
            rmse_pe  = rmse_vec(truth = emp_pe_obs, estimate = m_pe_obs, na_rm = T), 
            rmse_rpd  = rmse_vec(truth = emp_rpd_obs, estimate = m_rpd_obs, na_rm = T), 
            rmse_rpe  = rmse_vec(truth = emp_rpe_obs, estimate = m_rpe_obs, na_rm = T), 
            rmse_pd_rand_mean  = rmse_vec(truth = emp_pd_rand_mean, estimate = m_pd_rand_mean, na_rm = T), 
            rmse_pe_rand_mean  = rmse_vec(truth = emp_pe_rand_mean, estimate = m_pe_rand_mean, na_rm = T), 
            rmse_rpd_rand_mean  = rmse_vec(truth = emp_rpd_rand_mean, estimate = m_rpd_rand_mean, na_rm = T), 
            rmse_rpe_rand_mean  = rmse_vec(truth = emp_rpe_rand_mean, estimate = m_rpe_rand_mean, na_rm = T),
            rmse_pd_obs_z = rmse_vec(truth = emp_pd_obs_z, estimate = m_pd_obs_z, na_rm = T),
            rmse_pe_obs_z = rmse_vec(truth = emp_pe_obs_z, estimate = m_pe_obs_z, na_rm = T),
            rmse_rpd_obs_z = rmse_vec(truth = emp_rpd_obs_z, estimate = m_rpd_obs_z, na_rm = T),
            rmse_rpe_obs_z = rmse_vec(truth = emp_rpe_obs_z, estimate = m_rpe_obs_z, na_rm = T),
            
            mape_pd  = mape_vec(truth = emp_pd_obs, estimate = m_pd_obs, na_rm = T), 
            mape_pe  = mape_vec(truth = emp_pe_obs, estimate = m_pe_obs, na_rm = T), 
            mape_rpd  = mape_vec(truth = emp_rpd_obs, estimate = m_rpd_obs, na_rm = T), 
            mape_rpe  = mape_vec(truth = emp_rpe_obs, estimate = m_rpe_obs, na_rm = T), 
            mape_pd_rand_mean  = mape_vec(truth = emp_pd_rand_mean, estimate = m_pd_rand_mean, na_rm = T), 
            mape_pe_rand_mean  = mape_vec(truth = emp_pe_rand_mean, estimate = m_pe_rand_mean, na_rm = T), 
            mape_rpd_rand_mean  = mape_vec(truth = emp_rpd_rand_mean, estimate = m_rpd_rand_mean, na_rm = T), 
            mape_rpe_rand_mean  = mape_vec(truth = emp_rpe_rand_mean, estimate = m_rpe_rand_mean, na_rm = T),
            mape_pd_obs_z = mape_vec(truth = emp_pd_obs_z, estimate = m_pd_obs_z, na_rm = T),
            mape_pe_obs_z = mape_vec(truth = emp_pe_obs_z, estimate = m_pe_obs_z, na_rm = T),
            mape_rpd_obs_z = mape_vec(truth = emp_rpd_obs_z, estimate = m_rpd_obs_z, na_rm = T),
            mape_rpe_obs_z = mape_vec(truth = emp_rpe_obs_z, estimate = m_rpe_obs_z, na_rm = T),
            
            rsq_pd  = rsq_vec(truth = emp_pd_obs, estimate = m_pd_obs, na_rm = T),
            rsq_pe  = rsq_vec(truth = emp_pe_obs, estimate = m_pe_obs, na_rm = T),
            rsq_rpd  = rsq_vec(truth = emp_rpd_obs, estimate = m_rpd_obs, na_rm = T),
            rsq_rpe  = rsq_vec(truth = emp_rpe_obs, estimate = m_rpe_obs, na_rm = T),
            rsq_pd_rand_mean  = rsq_vec(truth = emp_pd_rand_mean, estimate = m_pd_rand_mean, na_rm = T),
            rsq_pe_rand_mean  = rsq_vec(truth = emp_pe_rand_mean, estimate = m_pe_rand_mean, na_rm = T),
            rsq_rpd_rand_mean = rsq_vec(truth = emp_rpd_rand_mean, estimate = m_rpd_rand_mean, na_rm = T),
            rsq_rpe_rand_mean = rsq_vec(truth = emp_rpe_rand_mean, estimate = m_rpe_rand_mean, na_rm = T),
            rsq_pd_obs_z  = rsq_vec(truth = emp_pd_obs_z, estimate = m_pd_obs_z, na_rm = T),
            rsq_pe_obs_z  = rsq_vec(truth = emp_pe_obs_z, estimate = m_pe_obs_z, na_rm = T),
            rsq_rpd_obs_z = rsq_vec(truth = emp_rpd_obs_z, estimate = m_rpd_obs_z, na_rm = T),
            rsq_rpe_obs_z = rsq_vec(truth = emp_rpe_obs_z, estimate = m_rpe_obs_z, na_rm = T)
            )

# metric per sim then mean ####
canaper_all_df <- left_join(canaper_res_all %>% 
                              dplyr::select(any_of(c(names(canaper_emp), "sim_n"))) %>% ungroup(), 
                            canaper_emp %>% ungroup() %>% 
                              dplyr::select(-c(treat:imputation)) %>% 
                              rename_if(is.numeric, ~paste0("emp_", .x)), by = join_by("site"))

stats_sim_df <- canaper_all_df %>% 
  group_by(sim_n, samp_frac, treat, imputation) %>%
  filter(pd_obs !=is.na(pd_obs)) %>% 
  summarise(mpe_pd  = mpe_vec(truth = emp_pd_obs, estimate = pd_obs, na_rm = T),
            mpe_pe  = mpe_vec(truth = emp_pe_obs, estimate = pe_obs, na_rm = T),
            mpe_rpd = mpe_vec(truth = emp_rpd_obs, estimate = rpd_obs, na_rm = T),
            mpe_rpe = mpe_vec(truth = emp_rpe_obs, estimate = rpe_obs, na_rm = T),
            mpe_pd_rand_mean  = mpe_vec(truth = emp_pd_rand_mean, estimate = pd_rand_mean, na_rm = T),
            mpe_pe_rand_mean  = mpe_vec(truth = emp_pe_rand_mean, estimate = pe_rand_mean, na_rm = T),
            mpe_rpd_rand_mean = mpe_vec(truth = emp_rpd_rand_mean, estimate = rpd_rand_mean, na_rm = T),
            mpe_rpe_rand_mean = mpe_vec(truth = emp_rpe_rand_mean, estimate = rpe_rand_mean, na_rm = T),
            mpe_pd_obs_z = mpe_vec(truth = emp_pd_obs_z, estimate = pd_obs_z, na_rm = T),
            mpe_pe_obs_z = mpe_vec(truth = emp_pe_obs_z, estimate = pe_obs_z, na_rm = T),
            mpe_rpd_obs_z = mpe_vec(truth = emp_rpd_obs_z, estimate = rpd_obs_z, na_rm = T),
            mpe_rpe_obs_z = mpe_vec(truth = emp_rpe_obs_z, estimate = rpe_obs_z, na_rm = T),
            
            rmse_pd  = rmse_vec(truth = emp_pd_obs, estimate = pd_obs, na_rm = T),
            rmse_pe  = rmse_vec(truth = emp_pe_obs, estimate = pe_obs, na_rm = T),
            rmse_rpd = rmse_vec(truth = emp_rpd_obs, estimate = rpd_obs, na_rm = T),
            rmse_rpe = rmse_vec(truth = emp_rpe_obs, estimate = rpe_obs, na_rm = T),
            rmse_pd_rand_mean  = rmse_vec(truth = emp_pd_rand_mean, estimate = pd_rand_mean, na_rm = T),
            rmse_pe_rand_mean  = rmse_vec(truth = emp_pe_rand_mean, estimate = pe_rand_mean, na_rm = T),
            rmse_rpd_rand_mean = rmse_vec(truth = emp_rpd_rand_mean, estimate = rpd_rand_mean, na_rm = T),
            rmse_rpe_rand_mean = rmse_vec(truth = emp_rpe_rand_mean, estimate = rpe_rand_mean, na_rm = T),
            rmse_pd_obs_z = rmse_vec(truth = emp_pd_obs_z, estimate = pd_obs_z, na_rm = T),
            rmse_pe_obs_z = rmse_vec(truth = emp_pe_obs_z, estimate = pe_obs_z, na_rm = T),
            rmse_rpd_obs_z = rmse_vec(truth = emp_rpd_obs_z, estimate = rpd_obs_z, na_rm = T),
            rmse_rpe_obs_z = rmse_vec(truth = emp_rpe_obs_z, estimate = rpe_obs_z, na_rm = T),
            
            mape_pd  = mape_vec(truth = emp_pd_obs, estimate = pd_obs, na_rm = T),
            mape_pe  = mape_vec(truth = emp_pe_obs, estimate = pe_obs, na_rm = T),
            mape_rpd = mape_vec(truth = emp_rpd_obs, estimate = rpd_obs, na_rm = T),
            mape_rpe = mape_vec(truth = emp_rpe_obs, estimate = rpe_obs, na_rm = T),
            mape_pd_rand_mean  = mape_vec(truth = emp_pd_rand_mean, estimate = pd_rand_mean, na_rm = T),
            mape_pe_rand_mean  = mape_vec(truth = emp_pe_rand_mean, estimate = pe_rand_mean, na_rm = T),
            mape_rpd_rand_mean = mape_vec(truth = emp_rpd_rand_mean, estimate = rpd_rand_mean, na_rm = T),
            mape_rpe_rand_mean = mape_vec(truth = emp_rpe_rand_mean, estimate = rpe_rand_mean, na_rm = T),
            mape_pd_obs_z = mape_vec(truth = emp_pd_obs_z, estimate = pd_obs_z, na_rm = T),
            mape_pe_obs_z = mape_vec(truth = emp_pe_obs_z, estimate = pe_obs_z, na_rm = T),
            mape_rpd_obs_z = mape_vec(truth = emp_rpd_obs_z, estimate = rpd_obs_z, na_rm = T),
            mape_rpe_obs_z = mape_vec(truth = emp_rpe_obs_z, estimate = rpe_obs_z, na_rm = T),
            
            rsq_pd  = rsq_vec(truth = emp_pd_obs, estimate = pd_obs, na_rm = T),
            rsq_pe  = rsq_vec(truth = emp_pe_obs, estimate = pe_obs, na_rm = T),
            rsq_rpd  = rsq_vec(truth = emp_rpd_obs, estimate = rpd_obs, na_rm = T),
            rsq_rpe  = rsq_vec(truth = emp_rpe_obs, estimate = rpe_obs, na_rm = T),
            rsq_pd_rand_mean  = rsq_vec(truth = emp_pd_rand_mean, estimate = pd_rand_mean, na_rm = T),
            rsq_pe_rand_mean  = rsq_vec(truth = emp_pe_rand_mean, estimate = pe_rand_mean, na_rm = T),
            rsq_rpd_rand_mean  = rsq_vec(truth = emp_rpd_rand_mean, estimate = rpd_rand_mean, na_rm = T),
            rsq_rpe_rand_mean  = rsq_vec(truth = emp_rpe_rand_mean, estimate = rpe_rand_mean, na_rm = T),
            rsq_pd_obs_z  = rsq_vec(truth = emp_pd_obs_z, estimate = pd_obs_z, na_rm = T),
            rsq_pe_obs_z  = rsq_vec(truth = emp_pe_obs_z, estimate = pe_obs_z, na_rm = T),
            rsq_rpd_obs_z = rsq_vec(truth = emp_rpd_obs_z, estimate = rpd_obs_z, na_rm = T),
            rsq_rpe_obs_z = rsq_vec(truth = emp_rpe_obs_z, estimate = rpe_obs_z, na_rm = T)
            )

# mean and sd tables ####
## mpe - relative bias ####
mpe_means <- stats_sim_df %>% 
  select(sim_n, samp_frac, treat, imputation, contains(match = c("mpe"))) %>% 
  ungroup() %>% group_by(samp_frac, treat, imputation) %>% 
  summarise(across(mpe_pd:mpe_rpe_obs_z, ~mean(.x, na.rm = T), .names = "mean_{.col}"),
            across(mpe_pd:mpe_rpe_obs_z, ~sd(.x, na.rm = T), .names = "sd_{.col}"))

mpe_means <- mpe_means %>% 
  mutate(prun_type = case_when(str_detect(treat, "late-pr") ~ "Lateral - Phylogenetic", 
                               str_detect(treat, "late-spa") ~ "Lateral - Spatial",
                               str_detect(treat, "marg-pr") ~ "Bilateral - Phylogenetic", 
                               str_detect(treat, "marg-spa") ~ "Bilateral - Spatial",
                               str_detect(treat, "rand-pr") ~ "Random", 
                               TRUE ~ "other" # Default value
                               ),
        prun3type = case_when(str_detect(treat, "late-pr") ~ "Lateral", 
                        str_detect(treat, "late-spa") ~ "Lateral",
                        str_detect(treat, "marg-pr") ~ "Bilateral", 
                        str_detect(treat, "marg-spa") ~ "Bilateral",
                        str_detect(treat, "rand-pr") ~ "Random", 
                        TRUE ~ "other" # Default value
                        ),
         prun_mhtd = case_when(str_detect(treat, "late-pr") ~ "Phylogenetic", 
                        str_detect(treat, "late-spa") ~ "Spatial",
                        str_detect(treat, "marg-pr") ~ "Phylogenetic", 
                        str_detect(treat, "marg-spa") ~ "Spatial",
                        str_detect(treat, "rand-pr") ~ "Random", 
                        TRUE ~ "other" # Default value
                        )
                     ) %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"), 
         is_rand = ifelse(treat == "rand-pr", "Random pruning", "Biased pruning"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  mutate(prun_type = factor(prun_type, levels = c("Random", "Trait", "Spatial"), labels = c("Random", "Phylogenetic", "Spatial")),
         treat = factor(treat, levels = c("rand-pr", "late-pr", "marg-pr", "late-spat_pr", "marg-spat_pr"), 
                        labels = c("Random", "Lateral - Phylogenetic", "Bilateral - Phylogenetic", "Lateral - Spatial", "Bilateral - Spatial")),
         imputation = factor(imputation, levels = c("none", "random", 'taxonomy - node', "taxonomy - random"), 
                             labels = c("No imputation", "Random", "At crown node", "Bellow crown node")) )

## rmse - absolute accuracy ####
rmse_means <- stats_sim_df %>% 
  select(sim_n, samp_frac, treat, imputation, contains(match = c("rmse"))) %>% 
  ungroup() %>% group_by(samp_frac, treat, imputation) %>% 
  summarise(across(rmse_pd:rmse_rpe_obs_z, ~mean(.x, na.rm = T), .names = "mean_{.col}"),
            across(rmse_pd:rmse_rpe_obs_z, ~sd(.x, na.rm = T), .names = "sd_{.col}"))

rmse_means <- rmse_means %>% 
  mutate(prun_type = case_when(str_detect(treat, "late-pr") ~ "Lateral - Phylogenetic", 
                               str_detect(treat, "late-spa") ~ "Lateral - Spatial",
                               str_detect(treat, "marg-pr") ~ "Bilateral - Phylogenetic", 
                               str_detect(treat, "marg-spa") ~ "Bilateral - Spatial",
                               str_detect(treat, "rand-pr") ~ "Random", 
                               TRUE ~ "other" # Default value
  ),
  prun3type = case_when(str_detect(treat, "late-pr") ~ "Lateral", 
                        str_detect(treat, "late-spa") ~ "Lateral",
                        str_detect(treat, "marg-pr") ~ "Bilateral", 
                        str_detect(treat, "marg-spa") ~ "Bilateral",
                        str_detect(treat, "rand-pr") ~ "Random", 
                        TRUE ~ "other" # Default value
  ),
  prun_mhtd = case_when(str_detect(treat, "late-pr") ~ "Phylogenetic", 
                        str_detect(treat, "late-spa") ~ "Spatial",
                        str_detect(treat, "marg-pr") ~ "Phylogenetic", 
                        str_detect(treat, "marg-spa") ~ "Spatial",
                        str_detect(treat, "rand-pr") ~ "Random", 
                        TRUE ~ "other" # Default value
  )
  ) %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"), 
         is_rand = ifelse(treat == "rand-pr", "Random pruning", "Biased pruning"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  mutate(prun_type = factor(prun_type, levels = c("Random", "Trait", "Spatial"), labels = c("Random", "Phylogenetic", "Spatial")),
         treat = factor(treat, levels = c("rand-pr", "late-pr", "marg-pr", "late-spat_pr", "marg-spat_pr"), 
                        labels = c("Random", "Lateral - Phylogenetic", "Bilateral - Phylogenetic", "Lateral - Spatial", "Bilateral - Spatial")),
         imputation = factor(imputation, levels = c("none", "random", 'taxonomy - node', "taxonomy - random"), 
                             labels = c("No imputation", "Random", "At crown node", "Bellow crown node")) )

## mape - relative accuracy ####
mape_means <- stats_sim_df %>% 
  select(sim_n, samp_frac, treat, imputation, contains(match = c("mape"))) %>% 
  ungroup() %>% group_by(samp_frac, treat, imputation) %>% 
  summarise(across(mape_pd:mape_rpe_obs_z, ~mean(.x, na.rm = T), .names = "mean_{.col}"),
            across(mape_pd:mape_rpe_obs_z, ~sd(.x, na.rm = T), .names = "sd_{.col}"))

mape_means <- mape_means %>% 
  mutate(prun_type = case_when(str_detect(treat, "late-pr") ~ "Lateral - Phylogenetic", 
                               str_detect(treat, "late-spa") ~ "Lateral - Spatial",
                               str_detect(treat, "marg-pr") ~ "Bilateral - Phylogenetic", 
                               str_detect(treat, "marg-spa") ~ "Bilateral - Spatial",
                               str_detect(treat, "rand-pr") ~ "Random", 
                               TRUE ~ "other" # Default value
  ),
  prun3type = case_when(str_detect(treat, "late-pr") ~ "Lateral", 
                        str_detect(treat, "late-spa") ~ "Lateral",
                        str_detect(treat, "marg-pr") ~ "Bilateral", 
                        str_detect(treat, "marg-spa") ~ "Bilateral",
                        str_detect(treat, "rand-pr") ~ "Random", 
                        TRUE ~ "other" # Default value
  ),
  prun_mhtd = case_when(str_detect(treat, "late-pr") ~ "Phylogenetic", 
                        str_detect(treat, "late-spa") ~ "Spatial",
                        str_detect(treat, "marg-pr") ~ "Phylogenetic", 
                        str_detect(treat, "marg-spa") ~ "Spatial",
                        str_detect(treat, "rand-pr") ~ "Random", 
                        TRUE ~ "other" # Default value
  )
  ) %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"), 
         is_rand = ifelse(treat == "rand-pr", "Random pruning", "Biased pruning"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  mutate(prun_type = factor(prun_type, levels = c("Random", "Trait", "Spatial"), labels = c("Random", "Phylogenetic", "Spatial")),
         treat = factor(treat, levels = c("rand-pr", "late-pr", "marg-pr", "late-spat_pr", "marg-spat_pr"), 
                        labels = c("Random", "Lateral - Phylogenetic", "Bilateral - Phylogenetic", "Lateral - Spatial", "Bilateral - Spatial")),
         imputation = factor(imputation, levels = c("none", "random", 'taxonomy - node', "taxonomy - random"), 
                             labels = c("No imputation", "Random", "At crown node", "Bellow crown node")) )

## r-squared - consistency ####
rsq_means <- stats_sim_df %>% 
  select(sim_n, samp_frac, treat, imputation, contains("rsq")) %>% 
  ungroup() %>% group_by(samp_frac, treat, imputation) %>%
  summarise(across(rsq_pd:rsq_rpe_obs_z, ~mean(.x, na.rm = T), .names = "mean_{.col}"),
            across(rsq_pd:rsq_rpe_obs_z, ~sd(.x, na.rm = T), .names = "sd_{.col}"),
            across(rsq_pd:rsq_rpe_obs_z, ~sum(!is.na(.x)), .names = "n_{.col}"))

rsq_means <- rsq_means %>% 
  mutate(prun_type = case_when(str_detect(treat, "late-pr") ~ "Lateral - Phylogenetic", 
                               str_detect(treat, "late-spa") ~ "Lateral - Spatial",
                               str_detect(treat, "marg-pr") ~ "Bilateral - Phylogenetic", 
                               str_detect(treat, "marg-spa") ~ "Bilateral - Spatial",
                               str_detect(treat, "rand-pr") ~ "Random", 
                               TRUE ~ "other" # Default value
  ),
  prun3type = case_when(str_detect(treat, "late-pr") ~ "Lateral", 
                        str_detect(treat, "late-spa") ~ "Lateral",
                        str_detect(treat, "marg-pr") ~ "Bilateral", 
                        str_detect(treat, "marg-spa") ~ "Bilateral",
                        str_detect(treat, "rand-pr") ~ "Random", 
                        TRUE ~ "other" # Default value
  ),
  prun_mhtd = case_when(str_detect(treat, "late-pr") ~ "Phylogenetic", 
                        str_detect(treat, "late-spa") ~ "Spatial",
                        str_detect(treat, "marg-pr") ~ "Phylogenetic", 
                        str_detect(treat, "marg-spa") ~ "Spatial",
                        str_detect(treat, "rand-pr") ~ "Random", 
                        TRUE ~ "other" # Default value
  )
  ) %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"), 
         is_rand = ifelse(treat == "rand-pr", "Random pruning", "Biased pruning"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  mutate(prun_type = factor(prun_type, levels = c("Random", "Trait", "Spatial"), labels = c("Random", "Phylogenetic", "Spatial")),
         treat = factor(treat, levels = c("rand-pr", "late-pr", "marg-pr", "late-spat_pr", "marg-spat_pr"), 
                        labels = c("Random", "Lateral - Phylogenetic", "Bilateral - Phylogenetic", "Lateral - Spatial", "Bilateral - Spatial")),
         imputation = factor(imputation, levels = c("none", "random", 'taxonomy - node', "taxonomy - random"), 
                             labels = c("No imputation", "Random", "At crown node", "Bellow crown node")) )

# Figures #####
library(patchwork)
dir.create(paste0("figures/", out_fold), recursive = T)

forms2print <- c(".png", ".pdf", ".svg") 

# Mean and sd ####
## MPE ####
### with random imputation #### 
## Observed
# PD
mpe_pd_sd_pt <- mpe_means %>%
  ggplot(aes(x = samp_frac, y = mean_mpe_pd, 
             ymin = mean_mpe_pd - sd_mpe_pd, 
             ymax = mean_mpe_pd + sd_mpe_pd,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = "PD", 
       x = "", y = "Mean Percent Error", 
    color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")
  
# PE
mpe_pe_sd_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pe, 
             ymin = mean_mpe_pe - sd_mpe_pe, 
             ymax = mean_mpe_pe + sd_mpe_pe,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = "PE", 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# RPD
mpe_rpd_sd_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpd, 
             ymin = mean_mpe_rpd - sd_mpe_rpd, 
             ymax = mean_mpe_rpd + sd_mpe_rpd,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = "RPD", 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# RPE
mpe_rpe_sd_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpe, 
             ymin = mean_mpe_rpe - sd_mpe_rpe, 
             ymax = mean_mpe_rpe + sd_mpe_rpe,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = "RPE", 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# SES
# PD ses
mpe_pd_ses_pt <- mpe_means %>% 
   ggplot(aes(x = samp_frac, y = mean_mpe_pd_obs_z, 
             ymin = mean_mpe_pd_obs_z - sd_mpe_pd_obs_z, 
             ymax = mean_mpe_pd_obs_z + sd_mpe_pd_obs_z,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("PD"["SES"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# PE ses
mpe_pe_ses_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pe_obs_z, 
             ymin = mean_mpe_pe_obs_z - sd_mpe_pe_obs_z, 
             ymax = mean_mpe_pe_obs_z + sd_mpe_pe_obs_z,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("PE"["SES"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# RPD ses
mpe_rpd_ses_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpd_obs_z, 
             ymin = mean_mpe_rpd_obs_z - sd_mpe_rpd_obs_z, 
             ymax = mean_mpe_rpd_obs_z + sd_mpe_rpd_obs_z,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("RPD"["SES"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# RPE ses
mpe_rpe_ses_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpe_obs_z, 
             ymin = mean_mpe_rpe_obs_z - sd_mpe_rpe_obs_z, 
             ymax = mean_mpe_rpe_obs_z + sd_mpe_rpe_obs_z,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("RPE"["SES"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# Randomized 
# randPD
mpe_pd_rand_sd_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pd_rand_mean, 
             ymin = mean_mpe_pd_rand_mean - sd_mpe_pd_rand_mean, 
             ymax = mean_mpe_pd_rand_mean + sd_mpe_pd_rand_mean,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("PD"["randomized"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# randPE
mpe_pe_rand_sd_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pe_rand_mean, 
             ymin = mean_mpe_pe_rand_mean - sd_mpe_pe_rand_mean, 
             ymax = mean_mpe_pe_rand_mean + sd_mpe_pe_rand_mean,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("PE"["randomized"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# randRPD
mpe_rpd_rand_sd_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpd_rand_mean, 
             ymin = mean_mpe_rpd_rand_mean - sd_mpe_rpd_rand_mean, 
             ymax = mean_mpe_rpd_rand_mean + sd_mpe_rpd_rand_mean,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("RPD"["randomized"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# randRPE
mpe_rpe_rand_sd_pt <- mpe_means %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpe_rand_mean, 
             ymin = mean_mpe_rpe_rand_mean - sd_mpe_rpe_rand_mean, 
             ymax = mean_mpe_rpe_rand_mean + sd_mpe_rpe_rand_mean,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("RPE"["randomized"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# panel plots
pl_long <- "
  AAEE
  AAEE
  BBFF
  BBFF
  CCGG
  CCGG
  DDHH
  DDHH"

# observed and ses metrics
(mpe_pd_sd_pt) + 
(mpe_pe_sd_pt) + 
(mpe_rpd_sd_pt) + 
(mpe_rpe_sd_pt) +
(mpe_pd_ses_pt) + 
(mpe_pe_ses_pt) + 
(mpe_rpd_ses_pt) + 
(mpe_rpe_ses_pt) +
plot_layout(guides = "collect", 
  design = pl_long, axis_titles = "collect", tag_level = "new") #+ plot_annotation(tag_levels = "a", tag_prefix = '(', tag_suffix = ")")

# loop over fig formats
for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/", out_fold, "/mpe_obs_ses_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}

### without random imputation #####
## Observed
# PD
mpe_pd_sd_pt <- mpe_means %>%
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pd, 
             ymin = mean_mpe_pd - sd_mpe_pd, 
             ymax = mean_mpe_pd + sd_mpe_pd,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = "PD", 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")

# PE
mpe_pe_sd_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pe, 
             ymin = mean_mpe_pe - sd_mpe_pe, 
             ymax = mean_mpe_pe + sd_mpe_pe,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = "PE", 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")

# RPD
mpe_rpd_sd_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpd, 
             ymin = mean_mpe_rpd - sd_mpe_rpd, 
             ymax = mean_mpe_rpd + sd_mpe_rpd,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = "RPD", 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")

# RPE
mpe_rpe_sd_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpe, 
             ymin = mean_mpe_rpe - sd_mpe_rpe, 
             ymax = mean_mpe_rpe + sd_mpe_rpe,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = "RPE", 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")

# SES
# PD ses
mpe_pd_ses_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pd_obs_z, 
             ymin = mean_mpe_pd_obs_z - sd_mpe_pd_obs_z, 
             ymax = mean_mpe_pd_obs_z + sd_mpe_pd_obs_z,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("PD"["SES"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")

# PE ses
mpe_pe_ses_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pe_obs_z, 
             ymin = mean_mpe_pe_obs_z - sd_mpe_pe_obs_z, 
             ymax = mean_mpe_pe_obs_z + sd_mpe_pe_obs_z,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("PE"["SES"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")

# RPD ses
mpe_rpd_ses_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpd_obs_z, 
             ymin = mean_mpe_rpd_obs_z - sd_mpe_rpd_obs_z, 
             ymax = mean_mpe_rpd_obs_z + sd_mpe_rpd_obs_z,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("RPD"["SES"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")

# RPE ses
mpe_rpe_ses_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpe_obs_z, 
             ymin = mean_mpe_rpe_obs_z - sd_mpe_rpe_obs_z, 
             ymax = mean_mpe_rpe_obs_z + sd_mpe_rpe_obs_z,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("RPE"["SES"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")

# Randomized 
# randPD
mpe_pd_rand_sd_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pd_rand_mean, 
             ymin = mean_mpe_pd_rand_mean - sd_mpe_pd_rand_mean, 
             ymax = mean_mpe_pd_rand_mean + sd_mpe_pd_rand_mean,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("PD"["randomized"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right")

# randPE
mpe_pe_rand_sd_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_pe_rand_mean, 
             ymin = mean_mpe_pe_rand_mean - sd_mpe_pe_rand_mean, 
             ymax = mean_mpe_pe_rand_mean + sd_mpe_pe_rand_mean,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("PE"["randomized"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")

# randRPD
mpe_rpd_rand_sd_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpd_rand_mean, 
             ymin = mean_mpe_rpd_rand_mean - sd_mpe_rpd_rand_mean, 
             ymax = mean_mpe_rpd_rand_mean + sd_mpe_rpd_rand_mean,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("RPD"["randomized"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")

# randRPE
mpe_rpe_rand_sd_pt <- mpe_means %>% 
  filter(imputation != "Random") %>% 
  ggplot(aes(x = samp_frac, y = mean_mpe_rpe_rand_mean, 
             ymin = mean_mpe_rpe_rand_mean - sd_mpe_rpe_rand_mean, 
             ymax = mean_mpe_rpe_rand_mean + sd_mpe_rpe_rand_mean,
             color = imputation, shape = treat)) +
  geom_line(linewidth = 0.9, alpha = 0.65) + 
  geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  labs(subtitle = expression("RPE"["randomized"]), 
       x = "", y = "Mean Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")

# panel plot
pl_long <- "
  AAEE
  AAEE
  BBFF
  BBFF
  CCGG
  CCGG
  DDHH
  DDHH"

# observed and ses metrics
(mpe_pd_sd_pt + theme(legend.position = "none")) +
(mpe_pe_sd_pt + theme(legend.position = "none")) +
(mpe_rpd_sd_pt + theme(legend.position = "none")) +
(mpe_rpe_sd_pt + theme(legend.position = "none")) +
(mpe_pd_ses_pt + theme(legend.position = "none")) +
(mpe_pe_ses_pt + theme(legend.position = "none")) +
(mpe_rpd_ses_pt + theme(legend.position = "none")) +
(mpe_rpe_ses_pt + theme(legend.position = "none")) +
  plot_layout(#guides = "keep", 
              design = pl_long, axis_titles = "collect", tag_level = "new")
  #guides(color = "none", shape = "none") #+ plot_annotation(tag_levels = "a", tag_prefix = '(', tag_suffix = ")")

# loop over fig formats
for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/", out_fold, "/mpe_obs_ses-rand_imp_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}

# all metrics panel
pl_long <- "
  AAEEII
  AAEEII
  BBFFJJ
  BBFFJJ
  CCGGKK
  CCGGKK
  DDHHLL
  DDHHLL"

# observed and ses metrics
(mpe_pd_sd_pt) + 
  (mpe_pe_sd_pt) + 
  (mpe_rpd_sd_pt) + 
  (mpe_rpe_sd_pt) +
  (mpe_pd_rand_sd_pt) + 
  (mpe_pe_rand_sd_pt) + 
  (mpe_rpd_rand_sd_pt) + 
  (mpe_rpe_rand_sd_pt) + 
  (mpe_pd_ses_pt) + 
  (mpe_pe_ses_pt) + 
  (mpe_rpd_ses_pt) + 
  (mpe_rpe_ses_pt) +
  plot_layout(guides = "collect", 
              design = pl_long, axis_titles = "collect", tag_level = "new")

# loop over fig formats
for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/", out_fold, "/all_metr_-rand_imp_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}

# MPE standard deviation trends 
mpe_means %>% 
  filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  select(samp_frac, treat, imputation, contains("sd_")) %>% 
  select(!contains("rand_")) %>% 
  filter(imputation != "none") %>% view()
  gather(key, value, -samp_frac, -treat, -imputation) %>% 
  ggplot(aes(x = samp_frac, y = value, 
             color = imputation, shape = treat)) +
  geom_line() + geom_point() +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  labs(x = "Sampling fraction", y = "MPE's SD") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~ key, scales = "free", ncol = 2)

# MAPE ####
### with random imputation #### 
# PD
mape_pd_sd_pt <- mape_means %>%
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_pd, 
             ymin = mean_mape_pd - sd_mape_pd, 
             ymax = mean_mape_pd + sd_mape_pd,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "PD", 
       x = "", y = "Mean Absolute Percentage Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
mape_pe_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_pe, 
             ymin = mean_mape_pe - sd_mape_pe, 
             ymax = mean_mape_pe + sd_mape_pe,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "PE", 
       x = "", y = "Mean Absolute Percentage Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
mape_rpd_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpd, 
             ymin = mean_mape_rpd - sd_mape_rpd, 
             ymax = mean_mape_rpd + sd_mape_rpd,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "RPD", 
       x = "", y = "Mean Absolute Percentage Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
mape_rpe_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpe, 
             ymin = mean_mape_rpe - sd_mape_rpe, 
             ymax = mean_mape_rpe + sd_mape_rpe,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "RPE", 
       x = "Proportion of tips removed", y = "Mean Absolute Percentage Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PD ses
mape_pd_ses_pt <- mape_means %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_pd_obs_z, 
             ymin = mean_mape_pd_obs_z - sd_mape_pd_obs_z, 
             ymax = mean_mape_pd_obs_z + sd_mape_pd_obs_z,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = expression("PD"["SES"]), 
       x = "", y = "Mean Absolute Percentage Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE ses
mape_pe_ses_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  ggplot(aes(x = samp_frac, y = mean_mape_pe_obs_z, 
             ymin = mean_mape_pe_obs_z - sd_mape_pe_obs_z, 
             ymax = mean_mape_pe_obs_z + sd_mape_pe_obs_z,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = expression("PE"["SES"]), 
       x = "", y = "Mean Absolute Percentage Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD ses
mape_rpd_ses_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpd_obs_z, 
             ymin = mean_mape_rpd_obs_z - sd_mape_rpd_obs_z, 
             ymax = mean_mape_rpd_obs_z + sd_mape_rpd_obs_z,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = expression("RPD"["SES"]), 
       x = "", y = "Mean Absolute Percentage Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) +
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE ses
mape_rpe_ses_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpe_obs_z, 
             ymin = mean_mape_rpe_obs_z - sd_mape_rpe_obs_z, 
             ymax = mean_mape_rpe_obs_z + sd_mape_rpe_obs_z,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = expression("RPE"["SES"]), 
       x = "Proportion of tips removed", y = "Mean Absolute Percentage Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)


pl_long <- "
  AAEE
  AAEE
  BBFF
  BBFF
  CCGG
  CCGG
  DDHH
  DDHH"

# observed and ses metrics
(mape_pd_sd_pt) + 
  (mape_pe_sd_pt) + 
  (mape_rpd_sd_pt) + 
  (mape_rpe_sd_pt) +
  (mape_pd_ses_pt) + 
  (mape_pe_ses_pt) + 
  (mape_rpd_ses_pt) + 
  (mape_rpe_ses_pt) +
  plot_layout(guides = "collect", 
              design = pl_long, axis_titles = "collect", tag_level = "new") #+ plot_annotation(tag_levels = "a", tag_prefix = '(', tag_suffix = ")")

# loop over fig formats
for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/", out_fold, "/mape_obs_ses_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}

## RMSE ####
### with random imputation #### 
# PD
rmse_pd_sd_pt <- rmse_means %>%
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_pd, 
             ymin = mean_rmse_pd - sd_rmse_pd, 
             ymax = mean_rmse_pd + sd_rmse_pd,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = "PD", 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
rmse_pe_sd_pt <- rmse_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_pe, 
             ymin = mean_rmse_pe - sd_rmse_pe, 
             ymax = mean_rmse_pe + sd_rmse_pe,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = "PE", 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
rmse_rpd_sd_pt <- rmse_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_rpd, 
             ymin = mean_rmse_rpd - sd_rmse_rpd, 
             ymax = mean_rmse_rpd + sd_rmse_rpd,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = "RPD", 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
rmse_rpe_sd_pt <- rmse_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_rpe, 
             ymin = mean_rmse_rpe - sd_rmse_rpe, 
             ymax = mean_rmse_rpe + sd_rmse_rpe,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = "RPE", 
       x = "Proportion of tips removed", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PD ses
rmse_pd_ses_pt <- rmse_means %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_pd_obs_z, 
             ymin = mean_rmse_pd_obs_z - sd_rmse_pd_obs_z, 
             ymax = mean_rmse_pd_obs_z + sd_rmse_pd_obs_z,
             color = imputation, shape = treat)) + 
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = expression("PD"["SES"]), 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE ses
rmse_pe_ses_pt <- rmse_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  ggplot(aes(x = samp_frac, y = mean_rmse_pe_obs_z, 
             ymin = mean_rmse_pe_obs_z - sd_rmse_pe_obs_z, 
             ymax = mean_rmse_pe_obs_z + sd_rmse_pe_obs_z,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = expression("PE"["SES"]), 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD ses
rmse_rpd_ses_pt <- rmse_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_rpd_obs_z, 
             ymin = mean_rmse_rpd_obs_z - sd_rmse_rpd_obs_z, 
             ymax = mean_rmse_rpd_obs_z + sd_rmse_rpd_obs_z,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = expression("RPD"["SES"]), 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) +
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE ses
rmse_rpe_ses_pt <- rmse_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_rpe_obs_z, 
             ymin = mean_rmse_rpe_obs_z - sd_rmse_rpe_obs_z, 
             ymax = mean_rmse_rpe_obs_z + sd_rmse_rpe_obs_z,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = expression("RPE"["SES"]), 
       x = "Proportion of tips removed", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

pl_long <- "
  AAEE
  AAEE
  BBFF
  BBFF
  CCGG
  CCGG
  DDHH
  DDHH"

# observed and ses metrics
(rmse_pd_sd_pt) + 
  (rmse_pe_sd_pt) + 
  (rmse_rpd_sd_pt) + 
  (rmse_rpe_sd_pt) +
  (rmse_pd_ses_pt) + 
  (rmse_pe_ses_pt) + 
  (rmse_rpd_ses_pt) + 
  (rmse_rpe_ses_pt) +
  plot_layout(guides = "collect", 
              design = pl_long, axis_titles = "collect", tag_level = "new") #+ plot_annotation(tag_levels = "a", tag_prefix = '(', tag_suffix = ")")

# loop over fig formats
for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/", out_fold, "/rmse_obs_ses_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}

## RMSE ####
### without random imputation #### 
# PD
rmse_pd_sd_pt <- rmse_means %>%
  filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_pd, 
             ymin = mean_rmse_pd - sd_rmse_pd, 
             ymax = mean_rmse_pd + sd_rmse_pd,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = "PD", 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3), shape = c(NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
rmse_pe_sd_pt <- rmse_means %>% 
  filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_pe, 
             ymin = mean_rmse_pe - sd_rmse_pe, 
             ymax = mean_rmse_pe + sd_rmse_pe,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = "PE", 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3), shape = c(NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
rmse_rpd_sd_pt <- rmse_means %>% 
  filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_rpd, 
             ymin = mean_rmse_rpd - sd_rmse_rpd, 
             ymax = mean_rmse_rpd + sd_rmse_rpd,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = "RPD", 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3), shape = c(NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
rmse_rpe_sd_pt <- rmse_means %>% 
  filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_rpe, 
             ymin = mean_rmse_rpe - sd_rmse_rpe, 
             ymax = mean_rmse_rpe + sd_rmse_rpe,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = "RPE", 
       x = "Proportion of tips removed", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3), shape = c(NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PD ses
rmse_pd_ses_pt <- rmse_means %>% 
  filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_pd_obs_z, 
             ymin = mean_rmse_pd_obs_z - sd_rmse_pd_obs_z, 
             ymax = mean_rmse_pd_obs_z + sd_rmse_pd_obs_z,
             color = imputation, shape = treat)) + 
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = expression("PD"["SES"]), 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3), shape = c(NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE ses
rmse_pe_ses_pt <- rmse_means %>% 
  filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  ggplot(aes(x = samp_frac, y = mean_rmse_pe_obs_z, 
             ymin = mean_rmse_pe_obs_z - sd_rmse_pe_obs_z, 
             ymax = mean_rmse_pe_obs_z + sd_rmse_pe_obs_z,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = expression("PE"["SES"]), 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3), shape = c(NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD ses
rmse_rpd_ses_pt <- rmse_means %>% 
  filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_rpd_obs_z, 
             ymin = mean_rmse_rpd_obs_z - sd_rmse_rpd_obs_z, 
             ymax = mean_rmse_rpd_obs_z + sd_rmse_rpd_obs_z,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = expression("RPD"["SES"]), 
       x = "", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3), shape = c(NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) +
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE ses
rmse_rpe_ses_pt <- rmse_means %>% 
  filter(imputation != "random") %>%
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rmse_rpe_obs_z, 
             ymin = mean_rmse_rpe_obs_z - sd_rmse_rpe_obs_z, 
             ymax = mean_rmse_rpe_obs_z + sd_rmse_rpe_obs_z,
             color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  labs(subtitle = expression("RPE"["SES"]), 
       x = "Proportion of tips removed", y = "Root Mean Square Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2, 3), shape = c(NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

pl_long <- "
  AAEE
  AAEE
  BBFF
  BBFF
  CCGG
  CCGG
  DDHH
  DDHH"

# observed and ses metrics
(rmse_pd_sd_pt) + 
  (rmse_pe_sd_pt) + 
  (rmse_rpd_sd_pt) + 
  (rmse_rpe_sd_pt) +
  (rmse_pd_ses_pt) + 
  (rmse_pe_ses_pt) + 
  (rmse_rpd_ses_pt) + 
  (rmse_rpe_ses_pt) +
  plot_layout(guides = "collect", 
              design = pl_long, axis_titles = "collect", tag_level = "new") #+ plot_annotation(tag_levels = "a", tag_prefix = '(', tag_suffix = ")")

# loop over fig formats
for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/", out_fold, "/rmse-rand_imp_ses_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}

## R2 ####
### with random imputation #### 
# PD
rsq_pd_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_pd, 
             ymin = mean_rsq_pd - sd_rsq_pd, 
             ymax = mean_rsq_pd + sd_rsq_pd,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "PD", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) +
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
rsq_pe_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_pe, 
             ymin = mean_rsq_pe - sd_rsq_pe, 
             ymax = mean_rsq_pe + sd_rsq_pe,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "PE", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
rsq_rpd_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpd, 
             ymin = mean_rsq_rpd - sd_rsq_rpd, 
             ymax = mean_rsq_rpd + sd_rsq_rpd,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "RPD", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
rsq_rpe_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpe, 
             ymin = mean_rsq_rpe - sd_rsq_rpe, 
             ymax = mean_rsq_rpe + sd_rsq_rpe,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "RPE", 
       x = "Proportion of tips removed", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PD ses
rsq_pd_ses_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  ggplot(aes(x = samp_frac, y = mean_rsq_pd_obs_z, 
             ymin = mean_rsq_pd_obs_z - sd_rsq_pd_obs_z, 
             ymax = mean_rsq_pd_obs_z + sd_rsq_pd_obs_z,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "PD ses", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE ses
rsq_pe_ses_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_pe_obs_z, 
             ymin = mean_rsq_pe_obs_z - sd_rsq_pe_obs_z, 
             ymax = mean_rsq_pe_obs_z + sd_rsq_pe_obs_z,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "PE ses", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD ses
rsq_rpd_ses_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpd_obs_z, 
             ymin = mean_rsq_rpd_obs_z - sd_rsq_rpd_obs_z, 
             ymax = mean_rsq_rpd_obs_z + sd_rsq_rpd_obs_z,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "RPD ses", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE ses
rsq_rpe_ses_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpe_obs_z, 
             ymin = mean_rsq_rpe_obs_z - sd_rsq_rpe_obs_z, 
             ymax = mean_rsq_rpe_obs_z + sd_rsq_rpe_obs_z,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "RPE ses", 
       x = "Proportion of tips removed", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPD
rsq_pd_rand_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  ggplot(aes(x = samp_frac, y = mean_rsq_pd_rand_mean, 
             ymin = mean_rsq_pd_rand_mean - sd_rsq_pd_rand_mean, 
             ymax = mean_rsq_pd_rand_mean + sd_rsq_pd_rand_mean,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Mean randomized PD", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPE
rsq_pe_rand_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_pe_rand_mean, 
             ymin = mean_rsq_pe_rand_mean - sd_rsq_pe_rand_mean, 
             ymax = mean_rsq_pe_rand_mean + sd_rsq_pe_rand_mean,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Mean randomized PE", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPD
rsq_rpd_rand_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpd_rand_mean, 
             ymin = mean_rsq_rpd_rand_mean - sd_rsq_rpd_rand_mean, 
             ymax = mean_rsq_rpd_rand_mean + sd_rsq_rpd_rand_mean,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Mean randomized RPD", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPE
rsq_rpe_rand_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpe_rand_mean, 
             ymin = mean_rsq_rpe_rand_mean - sd_rsq_rpe_rand_mean, 
             ymax = mean_rsq_rpe_rand_mean + sd_rsq_rpe_rand_mean,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Mean randomized RPE", 
       x = "Proportion of tips removed", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

pl_long <- 
 "AAEE
  AAEE
  BBFF
  BBFF
  CCGG
  CCGG
  DDHH
  DDHH"

(rsq_pd_sd_pt)  + #theme(legend.position = "none") +
(rsq_pe_sd_pt)  + #theme(legend.position = "none") +
(rsq_rpd_sd_pt) + #theme(legend.position = "none") +
(rsq_rpe_sd_pt) + #theme(legend.position = "none") +
(rsq_pd_ses_pt)  + #theme(legend.position = "none") +
(rsq_pe_ses_pt)  + #theme(legend.position = "none") +
(rsq_rpd_ses_pt) + #theme(legend.position = "none") +
(rsq_rpe_ses_pt) + #theme(legend.position = "none") +
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
ggsave(filename = paste0("figures/", out_fold, "/rsq_mean_sd-all_pnl", fig_out), 
       width = 30, height = 40, units = "cm", dpi = 200)
}

pl_long <- 
 "AAEEII
  AAEEII
  BBFFJJ
  BBFFJJ
  CCGGKK
  CCGGKK
  DDHHLL
  DDHHLL"

(rsq_pd_sd_pt)  + #theme(legend.position = "none") +
(rsq_pe_sd_pt)  + #theme(legend.position = "none") +
(rsq_rpd_sd_pt) + #theme(legend.position = "none") +
(rsq_rpe_sd_pt) + #theme(legend.position = "none") +
(rsq_pd_rand_sd_pt)  + #theme(legend.position = "none") +
(rsq_pe_rand_sd_pt)  + #theme(legend.position = "none") +
(rsq_rpd_rand_sd_pt) + #theme(legend.position = "none") +
(rsq_rpe_rand_sd_pt) + #theme(legend.position = "none") +
(rsq_pd_ses_pt)  + #theme(legend.position = "none") +
(rsq_pe_ses_pt)  + #theme(legend.position = "none") +
(rsq_rpd_ses_pt) + #theme(legend.position = "none") +
(rsq_rpe_ses_pt) + #theme(legend.position = "none") +
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")

# Only means ####
## MPE #####
# PD
mpe_pd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mpe_pd, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Phylogenetic Diversity (PD)", 
       x = "", y = "Mean Absolute Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  #scale_shape_discrete(guide = "none") +
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
mpe_pe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mpe_pe, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Phylogenetic Endemism (PE)", 
       x = "", y = "Mean Absolute Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
mpe_rpd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mpe_rpd, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Relative Phylogenetic Diversity (RPD)", 
       x = "", y = "Mean Absolute Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
mpe_rpe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mpe_rpe, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Relative Phylogenetic Endemism (RPE)", 
       x = "Proportion of tips removed", y = "Mean Absolute Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPD
mpe_rand_pd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mpe_pd_rand, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Randomized Mean Phylogenetic Diversity (randPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPE
mpe_rand_pe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mpe_pe_rand, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Randomized Mean Phylogenetic Endemism (randPE)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPD
mpe_rand_rpd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mpe_rpd_rand, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Randomized Mean Relative Phylogenetic Diversity (randRPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPE
mpe_rand_rpe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mpe_rpe_rand, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Randomized Mean Relative Phylogenetic Endemism (randRPE)", 
       x = "Proportion of tips removed", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)
#facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

pl_long <- "
  AAEE
  AAEE
  BBFF
  BBFF
  CCGG
  CCGG
  DDHH
  DDHH
  "
free(mpe_pd_pt) + free(mpe_pe_pt) + free(mpe_rpd_pt) + free(mpe_rpe_pt) +
  free(mpe_rand_pd_pt) + free(mpe_rand_pe_pt) + free(mpe_rand_rpd_pt) + free(mpe_rand_rpe_pt) +
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures", out_fold, "/mpe_mean-all_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}
  
## R2 ####
# PD
rsq_pd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = rsq_pd, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Phylogenetic Diversity (PD)", 
       x = "", y = expression("R"^2),
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
rsq_pe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = rsq_pe, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Phylogenetic Endemism (PE)", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
rsq_rpd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = rsq_rpd, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Relative Phylogenetic Diversity (RPD)", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
rsq_rpe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = rsq_rpe, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Relative Phylogenetic Endemism (RPE)", 
       x = "Proportion of tips removed", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPD
rsq_rand_pd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = rsq_pd_rand, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Randomized Mean Phylogenetic Diversity (randPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPE
rsq_rand_pe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = rsq_pe_rand, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Randomized Mean Phylogenetic Endemism (randPE)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPD
rsq_rand_rpd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = rsq_rpd_rand, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Randomized Mean Relative Phylogenetic Diversity (randRPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPE
rsq_rand_rpe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = rsq_rpe_rand, 
             color = imputation, shape = treat)) +
  geom_point() + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Randomized Mean Relative Phylogenetic Endemism (randRPE)", 
       x = "Proportion of tips removed", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# one panel
pl_long <- "
  AAEE
  AAEE
  BBFF
  BBFF
  CCGG
  CCGG
  DDHH
  DDHH
  "
free(rsq_pd_pt) + free(rsq_pe_pt) + free(rsq_rpd_pt) + free(rsq_rpe_pt) +
  free(rsq_rand_pd_pt) + free(rsq_rand_pe_pt) + free(rsq_rand_rpd_pt) + free(rsq_rand_rpe_pt) +
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
ggsave(filename = paste0("figures", out_fold, "/rsq_mean-all_pnl", fig_out), 
       width = 30, height = 40, units = "cm", dpi = 200)
}

mpe_means %>% 
  dplyr::select(samp_frac, treat, imputation, starts_with("mean"))%>% # mean sd
  gather(key, value, -samp_frac, -treat, -imputation) %>% 
  #filter(treat != "rand-pr" #, imputation == "none" #, samp_frac == 0.3) %>% 
  ggplot(aes(x = as.factor(samp_frac), y = value, col = treat)) +
  geom_boxplot(outlier.colour = "red", varwidth = T, na.rm = T) +
  #geom_point(col = "coral", cex = .5) + 
  theme_bw() + 
  labs(x = "Proportion of tips removed", y = "mean MPE") +
  geom_hline(yintercept = 20, col = "darkgray") +
  #geom_hline(yintercept = 10, col = "steelblue", alpha = .5) +
  facet_wrap(~key, scales = "fixed", ncol = 2)

rsq_means %>% 
  dplyr::select(samp_frac, treat, imputation, starts_with("mean"))%>% # mean sd
  gather(key, value, -samp_frac, -treat, -imputation) %>% 
  #filter(treat != "rand-pr" #, imputation == "none" #, samp_frac == 0.3) %>%
  ggplot(aes(x = as.factor(samp_frac), y = value, col = treat)) +
  geom_boxplot(outlier.colour = "red", varwidth = T, na.rm = T) +
  #geom_point(col = "coral", cex = .5) + 
  theme_bw() + 
  labs(x = "Proportion of tips removed", y = "mean R2") +
  geom_hline(yintercept = 0.6, col = "darkgray") +
  facet_wrap(~key, scales = "fixed", ncol = 2)

# get descriptive numbers ####
mpe_means %>% 
  dplyr::select(samp_frac, treat, imputation, starts_with("mean")) %>% # mean sd
  gather(key, value, -samp_frac, -treat, -imputation) %>% 
  filter(#treat != "rand-pr" 
    imputation == "none" #, samp_frac == 0.3
         ) %>% #head()
  #ggplot(aes(x = samp_frac, y = value)) + geom_point()
  ungroup() %>% 
  group_by(key , samp_frac, 
    treat#, imputation
           ) %>% 
  summarise(min = min(value, na.rm = T) %>% round(1), 
            max = max(value, na.rm = T) %>% round(1),
            mean = mean(value, na.rm = T) %>% round(1)) %>% 
  filter(grepl("_obs_z", key) ) %>% view

rsq_means %>% 
  dplyr::select(samp_frac, treat, imputation, starts_with("mean"))%>% # mean
  gather(key, value, -samp_frac, -treat, -imputation) %>%
  filter(#treat != "rand-pr" , #imputation == "none" #, samp_frac == 0.3
         ) %>% ungroup() %>% 
  group_by(key, samp_frac, treat #, imputation
           ) %>%
  summarize(min = min(value, na.rm = T), 
            max = max(value, na.rm = T),
            mean = mean(value, na.rm = T)) %>% 
    filter(samp_frac == 0.5) %>% view
