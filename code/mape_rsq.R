library(tidyverse)
library(yardstick)

# load and get canaper tables ####
# empirical full data
cols2sel <- c("site", 
              "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              "pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean")

canaper_emp_raw <- readRDS("output/canaper/canaper_empirical.rds")
canaper_emp <- canaper_emp_raw %>% 
  dplyr::select(any_of(cols2sel)) %>% 
  mutate(treat = "empirical", samp_frac = 1, imputation = "none")

# results data - analyses
canaper_res_all <- readRDS("output/canaper_all_data_res.rds")

# one value per sim ####
canaper_all_df <- left_join(canaper_res_all %>% ungroup(), 
                            canaper_emp %>% ungroup() %>% 
                              dplyr::select(-c(treat:imputation)) %>% 
                              rename_if(is.numeric, ~paste0("emp_", .x)), by = join_by("site"))

stats_sim_df <- canaper_all_df %>% 
  group_by(sim_n, samp_frac, treat, imputation) %>% 
  filter(pd_obs !=is.na(pd_obs)) %>% 
  summarise(mape_pd  = mape_vec(truth = emp_pd_obs, estimate = pd_obs, na_rm = T),
            mape_pe  = mape_vec(truth = emp_pe_obs, estimate = pe_obs, na_rm = T),
            mape_rpd  = mape_vec(truth = emp_rpd_obs, estimate = rpd_obs, na_rm = T),
            mape_rpe  = mape_vec(truth = emp_rpe_obs, estimate = rpe_obs, na_rm = T),
            mape_pd_rand  = mape_vec(truth = emp_pd_rand_mean, estimate = pd_rand_mean, na_rm = T),
            mape_pe_rand  = mape_vec(truth = emp_pe_rand_mean, estimate = pe_rand_mean, na_rm = T),
            mape_rpd_rand  = mape_vec(truth = emp_rpd_rand_mean, estimate = rpd_rand_mean, na_rm = T),
            mape_rpe_rand  = mape_vec(truth = emp_rpe_rand_mean, estimate = rpe_rand_mean, na_rm = T),
            
            rsq_pd  = rsq_vec(truth = emp_pd_obs, estimate = pd_obs, na_rm = T),
            rsq_pe  = rsq_vec(truth = emp_pe_obs, estimate = pe_obs, na_rm = T),
            rsq_rpd  = rsq_vec(truth = emp_rpd_obs, estimate = rpd_obs, na_rm = T),
            rsq_rpe  = rsq_vec(truth = emp_rpe_obs, estimate = rpe_obs, na_rm = T),
            rsq_pd_rand  = rsq_vec(truth = emp_pd_rand_mean, estimate = pd_rand_mean, na_rm = T),
            rsq_pe_rand  = rsq_vec(truth = emp_pe_rand_mean, estimate = pe_rand_mean, na_rm = T),
            rsq_rpd_rand  = rsq_vec(truth = emp_rpd_rand_mean, estimate = rpd_rand_mean, na_rm = T),
            rsq_rpe_rand  = rsq_vec(truth = emp_rpe_rand_mean, estimate = rpe_rand_mean, na_rm = T)
  )

# one panel for all
# MAPE
stats_sim_df %>% 
  dplyr::select(sim_n, samp_frac, treat, imputation, starts_with("mape") #& !ends_with("rand")
                ) %>%
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  gather(key, value, -sim_n, -samp_frac, -treat, -imputation, -is_imp, -cat2map) %>% 
  filter(is_imp == "Only pruned") %>% 
  ggplot(aes(x = (samp_frac), y = value, color = treat, shape = imputation)) +
  stat_summary(fun.data = "mean_sdl", geom = "pointrange", position = position_jitter(0.025)) +
  scale_shape_manual(values = c(0, 1, 2, 5)) + geom_line(aes(linetype = imputation)) + 
  #geom_smooth(aes(linetype = imputation), linewidth = .25, se = F) +
  labs(#title = "Phylogenetic Diversity (PD)", 
       x = "Proportion of tips removed", y = "Mean Absolute Percent Error", color = "Pruning\nmethod", shape = "Imputation", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  scale_colour_viridis_d(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), option = "inferno", begin = .8, end = .05, alpha = .95) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  facet_wrap(~key*is_imp, labeller = , ncol = 2, scales = "free")

#  R2
stats_sim_df %>% 
  dplyr::select(sim_n, samp_frac, treat, imputation, starts_with("rsq") & !ends_with("rand")) %>%
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  gather(key, value, -sim_n, -samp_frac, -treat, -imputation, -is_imp, -cat2map) %>%
  ggplot(aes(x = (samp_frac), y = value, color = treat, shape = imputation)) +
  stat_summary(fun.data = "mean_sdl", geom = "pointrange", position = position_jitter(0.09)) +
  scale_shape_manual(values = c(0, 1, 2, 5)) + #geom_line(aes(linetype = imputation)) + 
  #geom_smooth() +
  labs(#title = "Phylogenetic Diversity (PD)", 
    x = "Proportion of tips removed", y = "R2", color = "Pruning\nmethod", shape = "Imputation", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  scale_colour_viridis_d(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), option = "inferno", begin = .8, end = .05, alpha = .95) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  facet_wrap(~key*is_imp, labeller = , ncol = 2, scales = "free")

# mean and sd table
mape_means <- stats_sim_df %>% 
  select(sim_n, samp_frac, treat,imputation, contains("mape")) %>% 
  ungroup() %>% group_by(samp_frac, treat, imputation) %>%
  summarise(across(mape_pd:mape_rpe_rand, ~mean(.x, na.rm = T), .names = "mean_{.col}"),
            across(mape_pd:mape_rpe_rand, ~sd(.x, na.rm = T), .names = "sd_{.col}"))

rsq_means <- stats_sim_df %>% 
  select(sim_n, samp_frac, treat,imputation, contains("rsq")) %>% 
  ungroup() %>% group_by(samp_frac, treat, imputation) %>%
  summarise(across(rsq_pd:rsq_rpe_rand, ~mean(.x, na.rm = T), .names = "mean_{.col}"),
            across(rsq_pd:rsq_rpe_rand, ~sd(.x, na.rm = T), .names = "sd_{.col}"))

# mean values across sims #### 
# get mean across 100 sims per samp fraction, pruning and imputation methods
canaper_res_df <- canaper_res_all %>% 
  group_by(site, samp_frac, treat, imputation) %>%
  summarise(across(pd_obs:rpe_rand_mean, ~mean(.x, na.rm = T), .names = "m_{.col}")) #%>% dplyr::select(names(canaper_emp))

# build df 
avg_df <- left_join(canaper_res_df %>% ungroup(), 
                    canaper_emp %>% ungroup() %>% 
                      dplyr::select(-c(treat:imputation)) %>% 
                      rename_if(is.numeric, ~paste0("emp_", .x)), by = join_by("site")) 

m_stats_df <- avg_df %>% 
  #filter(pd_obs !=is.na(pd_obs)) %>% 
  group_by(samp_frac, treat, imputation) %>% 
  summarise(mape_pd  = mape_vec(truth = emp_pd_obs, estimate = m_pd_obs, na_rm = T), 
            mape_pe  = mape_vec(truth = emp_pe_obs, estimate = m_pe_obs, na_rm = T), 
            mape_rpd  = mape_vec(truth = emp_rpd_obs, estimate = m_rpd_obs, na_rm = T), 
            mape_rpe  = mape_vec(truth = emp_rpe_obs, estimate = m_rpe_obs, na_rm = T), 
            mape_pd_rand  = mape_vec(truth = emp_pd_rand_mean, estimate = m_pd_rand_mean, na_rm = T), 
            mape_pe_rand  = mape_vec(truth = emp_pe_rand_mean, estimate = m_pe_rand_mean, na_rm = T), 
            mape_rpd_rand  = mape_vec(truth = emp_rpd_rand_mean, estimate = m_rpd_rand_mean, na_rm = T), 
            mape_rpe_rand  = mape_vec(truth = emp_rpe_rand_mean, estimate = m_rpe_rand_mean, na_rm = T),
            
            rsq_pd  = rsq_vec(truth = emp_pd_obs, estimate = m_pd_obs, na_rm = T),
            rsq_pe  = rsq_vec(truth = emp_pe_obs, estimate = m_pe_obs, na_rm = T),
            rsq_rpd  = rsq_vec(truth = emp_rpd_obs, estimate = m_rpd_obs, na_rm = T),
            rsq_rpe  = rsq_vec(truth = emp_rpe_obs, estimate = m_rpe_obs, na_rm = T),
            rsq_pd_rand  = rsq_vec(truth = emp_pd_rand_mean, estimate = m_pd_rand_mean, na_rm = T),
            rsq_pe_rand  = rsq_vec(truth = emp_pe_rand_mean, estimate = m_pe_rand_mean, na_rm = T),
            rsq_rpd_rand  = rsq_vec(truth = emp_rpd_rand_mean, estimate = m_rpd_rand_mean, na_rm = T),
            rsq_rpe_rand  = rsq_vec(truth = emp_rpe_rand_mean, estimate = m_rpe_rand_mean, na_rm = T)
  )

# Figures #####
library(patchwork)

forms2print <- c(".png", ".pdf", ".svg") 

# Mean and sd vals ####
## MAPE ####
# PD
mape_pd_sd_pt <- mape_means %>%
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_pd, ymin = mean_mape_pd - sd_mape_pd, ymax = mean_mape_pd + sd_mape_pd,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Phylogenetic Diversity (PD)", 
       x = "", y = "Mean Absolute Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
mape_pe_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_pe, ymin = mean_mape_pe - sd_mape_pe, ymax = mean_mape_pe + sd_mape_pe,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Phylogenetic Endemism (PE)", 
       x = "", y = "Mean Absolute Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
mape_rpd_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpd, ymin = mean_mape_rpd - sd_mape_rpd, ymax = mean_mape_rpd + sd_mape_rpd,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Relative Phylogenetic Diversity (RPD)", 
       x = "", y = "Mean Absolute Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
mape_rpe_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpe, ymin = mean_mape_rpe - sd_mape_rpe, ymax = mean_mape_rpe + sd_mape_rpe,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Relative Phylogenetic Endemism (RPE)", 
       x = "Proportion of tips removed", y = "Mean Absolute Percent Error", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPD
mape_pd_rand_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_pd_rand, ymin = mean_mape_pd_rand - sd_mape_pd_rand, ymax = mean_mape_pd_rand + sd_mape_pd_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(subtitle = "Randomized Mean Phylogenetic Diversity (randPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPE
mape_pe_rand_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_pe_rand, ymin = mean_mape_pe_rand - sd_mape_pe_rand, ymax = mean_mape_pe_rand + sd_mape_pe_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Phylogenetic Endemism (randPE)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPD
mape_rpd_rand_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpd_rand, ymin = mean_mape_rpd_rand - sd_mape_rpd_rand, ymax = mean_mape_rpd_rand + sd_mape_rpd_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Relative Phylogenetic Diversity (randRPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPE
mape_rpe_rand_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpe_rand, ymin = mean_mape_rpe_rand - sd_mape_rpe_rand, ymax = mean_mape_rpe_rand + sd_mape_rpe_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Relative Phylogenetic Endemism (randRPE)", 
       x = "Proportion of tips removed", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPD
mape_pd_rand_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_pd_rand, ymin = mean_mape_pd_rand - sd_mape_pd_rand, ymax = mean_mape_pd_rand + sd_mape_pd_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Phylogenetic Diversity (randPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPE
mape_pe_rand_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_pe_rand, ymin = mean_mape_pe_rand - sd_mape_pe_rand, ymax = mean_mape_pe_rand + sd_mape_pe_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Phylogenetic Endemism (randPE)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPD
mape_rpd_rand_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpd_rand, ymin = mean_mape_rpd_rand - sd_mape_rpd_rand, ymax = mean_mape_rpd_rand + sd_mape_rpd_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Relative Phylogenetic Diversity (randRPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPE
mape_rpe_rand_sd_pt <- mape_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_mape_rpe_rand, ymin = mean_mape_rpe_rand - sd_mape_rpe_rand, ymax = mean_mape_rpe_rand + sd_mape_rpe_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Relative Phylogenetic Endemism (randRPE)", 
       x = "Proportion of tips removed", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

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
free(mape_pd_sd_pt) + free(mape_pe_sd_pt) + free(mape_rpd_sd_pt) + free(mape_rpe_sd_pt) +
  free(mape_pd_rand_sd_pt) + free(mape_pe_rand_sd_pt) + free(mape_rpd_rand_sd_pt) + free(mape_rpe_rand_sd_pt) +
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/mape_mean_sd-all_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}
  
## R2 ####
# PD
rsq_pd_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_pd, ymin = mean_rsq_pd - sd_rsq_pd, ymax = mean_rsq_pd + sd_rsq_pd,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Phylogenetic Diversity (PD)", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
rsq_pe_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_pe, ymin = mean_rsq_pe - sd_rsq_pe, ymax = mean_rsq_pe + sd_rsq_pe,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Phylogenetic Endemism (PE)", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
rsq_rpd_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpd, ymin = mean_rsq_rpd - sd_rsq_rpd, ymax = mean_rsq_rpd + sd_rsq_rpd,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Relative Phylogenetic Diversity (RPD)", 
       x = "", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
rsq_rpe_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpe, ymin = mean_rsq_rpe - sd_rsq_rpe, ymax = mean_rsq_rpe + sd_rsq_rpe,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Relative Phylogenetic Endemism (RPE)", 
       x = "Proportion of tips removed", y = expression("R"^2), 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right", axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPD
rsq_pd_rand_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_pd_rand, ymin = mean_rsq_pd_rand - sd_rsq_pd_rand, ymax = mean_rsq_pd_rand + sd_rsq_pd_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Phylogenetic Diversity (randPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPE
rsq_pe_rand_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_pe_rand, ymin = mean_rsq_pe_rand - sd_rsq_pe_rand, ymax = mean_rsq_pe_rand + sd_rsq_pe_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Phylogenetic Endemism (randPE)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPD
rsq_rpd_rand_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpd_rand, ymin = mean_rsq_rpd_rand - sd_rsq_rpd_rand, ymax = mean_rsq_rpd_rand + sd_rsq_rpd_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Relative Phylogenetic Diversity (randRPD)", 
       x = "", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPE
rsq_rpe_rand_sd_pt <- rsq_means %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mean_rsq_rpe_rand, ymin = mean_rsq_rpe_rand - sd_rsq_rpe_rand, ymax = mean_rsq_rpe_rand + sd_rsq_rpe_rand,
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) +
  labs(subtitle = "Randomized Mean Relative Phylogenetic Endemism (randRPE)", 
       x = "Proportion of tips removed", y = "", 
       color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3, 4), shape = c(NA, NA, NA, NA))), 
         shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8)))) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

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

free(rsq_pd_sd_pt)  + #theme(legend.position = "none") +
free(rsq_pe_sd_pt)  + #theme(legend.position = "none") +
free(rsq_rpd_sd_pt) + #theme(legend.position = "none") +
free(rsq_rpe_sd_pt) + #theme(legend.position = "none") +
free(rsq_pd_rand_sd_pt)  + #theme(legend.position = "none") +
free(rsq_pe_rand_sd_pt)  + #theme(legend.position = "none") +
free(rsq_rpd_rand_sd_pt) + #theme(legend.position = "none") +
free(rsq_rpe_rand_sd_pt) + #theme(legend.position = "none") +
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
ggsave(filename = paste0("figures/rsq_mean_sd-all_pnl", fig_out), 
       width = 30, height = 40, units = "cm", dpi = 200)
}

# Only mean vals ####
## MAPE #####
# PD
mape_pd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mape_pd, 
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
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
mape_pe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mape_pe, 
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
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
mape_rpd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mape_rpd, 
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
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
mape_rpe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mape_rpe, 
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
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) #+ facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPD
mape_rand_pd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mape_pd_rand, 
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
mape_rand_pe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mape_pe_rand, 
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
mape_rand_rpd_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mape_rpd_rand, 
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
mape_rand_rpe_pt <- m_stats_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = mape_rpe_rand, 
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
free(mape_pd_pt) + free(mape_pe_pt) + free(mape_rpd_pt) + free(mape_rpe_pt) +
  free(mape_rand_pd_pt) + free(mape_rand_pe_pt) + free(mape_rand_rpd_pt) + free(mape_rand_rpe_pt) +
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/mape_mean-all_pnl", fig_out), 
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
ggsave(filename = paste0("figures/rsq_mean-all_pnl", fig_out), 
       width = 30, height = 40, units = "cm", dpi = 200)
}

# get numbers and descriptives ####
mape_means %>% 
  dplyr::select(samp_frac, treat, imputation, starts_with("mean"))%>% # mean sd
  gather(key, value, -samp_frac, -treat, -imputation) %>%
  filter(imputation != "none" #, samp_frac == 0.3
         ) %>% #view()
  #ggplot(aes(x = samp_frac, y = value)) + geom_point()
  group_by(key) %>% 
  summarise(min = min(value, na.rm = T), 
            max = max(value, na.rm = T)) %>% arrange(desc(max))

rsq_means %>% 
  dplyr::select(samp_frac, treat, imputation, starts_with("mean"))%>% # mean
  gather(key, value, -samp_frac, -treat, -imputation) %>%
  filter(imputation == "none") %>% 
  ggplot(aes(x = samp_frac, y = value)) + geom_point()
  group_by(key) %>% 
  summarise(min = min(value, na.rm = T), 
            max = max(value, na.rm = T)) %>% arrange(desc(max))
