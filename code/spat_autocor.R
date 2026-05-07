library(tidyverse)

# Nitta et al. spatial data
spat_df <- read.csv("data/nitta-et-al/results/japan_ferns_biodiv.csv", header = T) %>%
   mutate(site = as.character(grids)) %>% dplyr::select(site, lat, long)

# load and get canaper tables ####
# empirical full data
cols2sel <- c("site", 
              "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              "pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_obs_z", "pe_obs_z", "rpd_obs_z", "rpe_obs_z")

canaper_emp_raw <- readRDS("output/canaper/canaper_empirical.rds")
canaper_emp <- canaper_emp_raw %>% 
  dplyr::select(any_of(cols2sel)) %>% 
  mutate(treat = "empirical", samp_frac = 0, imputation = "none") %>% 
  left_join(., spat_df, by = join_by("site"))
head(canaper_emp)

# results data - analyses 
canaper_res_all <- readRDS("output/canaper_all_data_res.rds")

# Moran's I ####
# estimate moran's I across pruning methods, fractions and imputation methods
# empirical data
moran_emp <- canaper_emp %>% 
  filter(!is.na(long)) %>%
  group_by(samp_frac, treat, imputation) %>% 
  summarise(across(pd_obs:rpe_obs_z, 
                   ~ape::Moran.I(.x, weight = as.matrix(1/dist(cbind(long, lat))), na.rm = T)$observed, 
                   .names = "moran_{.col}"))

# get moran's I for sim per samp fraction, pruning and imputation methods
#canaper_mean_res_df <- canaper_res_all %>%
#  group_by(site, samp_frac, treat, imputation) %>%
#  summarise(across(pd_obs:rpe_rand_mean, ~mean(.x, na.rm = T), .names = "m_{.col}"))
#
#canaper_mean_spat_df <- left_join(canaper_mean_res_df %>% ungroup(), spat_df, by = join_by("site"))


# get estimate's means across 100 sims per samp fraction, pruning and imputation methods
canaper_mean_res_df <- canaper_res_all %>%
  group_by(site, samp_frac, treat, imputation) %>%
  summarise(across(pd_obs:rpe_obs_z, ~mean(.x, na.rm = T), .names = "m_{.col}"))

canaper_mean_spat_df <- left_join(canaper_mean_res_df %>% ungroup(), spat_df, by = join_by("site"))

# get mean moran's I
run_moran = F
if(run_moran){
autocor_mean_df <- canaper_mean_spat_df %>% 
  filter(!is.na(long)) %>% 
  group_by(samp_frac, treat, imputation) %>% 
  summarise(across(m_pd_obs:m_rpe_obs_z, 
                ~ape::Moran.I(.x, weight = as.matrix(1/dist(cbind(long, lat))), na.rm = T)$observed, 
                .names = "moran_{.col}"))
autocor_mean_df %>% saveRDS("output/autocor_mean_df.rds") ; autocor_mean_df %>% head()
} else {
  autocor_mean_df <- readRDS("output/autocor_mean_df.rds")
}

# add categories for plots
autocor_mean_df <- autocor_mean_df %>% 
  mutate(prun4type = case_when(str_detect(treat, "late-pr") ~ "Lateral", 
                               str_detect(treat, "late-spa") ~ "Lateral - Spatial",
                               str_detect(treat, "marg-pr") ~ "Bilateral", 
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
  prun_type = case_when(str_detect(treat, "late-pr") ~ "Trait", 
                        str_detect(treat, "late-spa") ~ "Spatial",
                        str_detect(treat, "marg-pr") ~ "Trait", 
                        str_detect(treat, "marg-spa") ~ "Spatial",
                        str_detect(treat, "rand-pr") ~ "Random", 
                        TRUE ~ "other" # Default value
                        )) %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"), 
         is_rand = ifelse(treat == "rand-pr", "Random pruning", "Biased pruning"),
         cat2map = paste0(treat, "-" , imputation))


# figures #####
library(patchwork)
forms2print <- c(".png", ".pdf", ".svg") 

# mean and sd vals ####

# mean value across sims ####
# plots by metric
p_size = 3.5

# PD
moran_pd_pt <- autocor_mean_df %>%
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_pd_obs, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pd_obs), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Phylogenetic Diversity (PD)", 
       x = "", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)


# PE
moran_pe_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_pe_obs, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pe_obs), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Phylogenetic Endemism (PE)", 
       x = "", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)


# RPD
moran_rpd_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_rpd_obs, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpd_obs), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Relative Phylogenetic Diversity (RPD)", 
       x = "", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
  facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)


# RPE
moran_rpe_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_rpe_obs, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpe_obs), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Relative Phylogenetic Endemism (RPE)", 
       x = "Proportion of tips removed", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)

# randPD
moran_rand_pd_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_pd_rand_mean, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pd_rand_mean), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Mean Randomized Phylogenetic Diversity (randPD)", 
       x = "", y = "", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPE
moran_rand_pe_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_pe_rand_mean, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pe_rand_mean), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Mean Randomized Phylogenetic Endemism (randPE)", 
       x = "", y = "", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPD
moran_rand_rpd_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_rpd_rand_mean, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpd_rand_mean), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Mean Randomized Relative Phylogenetic Diversity (randRPD)", 
       x = "", y = "", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randRPE
moran_rand_rpe_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_rpe_rand_mean, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpe_rand_mean), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Mean Randomized Relative Phylogenetic Endemism (randRPE)", 
       x = "Proportion of tips removed", y = "", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PD ses
moran_pd_ses_pt <- autocor_mean_df %>%
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  ggplot(aes(x = samp_frac, y = moran_m_pd_obs_z, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pd_obs_z), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "PD ses", 
       x = "", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)

# PE ses
moran_pe_ses_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_pe_obs_z, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pe_obs_z), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "PE ses", 
       x = "", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)

# RPD ses
moran_rpd_ses_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_rpd_obs_z, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpd_obs_z), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "RPD ses", 
       x = "", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)

# RPE ses
moran_rpe_ses_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_rpe_obs_z, color = imputation, shape = treat)) +
  geom_line(aes(linetype = imputation), alpha = 0.5) + 
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpe_obs_z), colour = "black", pch = 18, size = p_size) +
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "RPE ses", 
       x = "Proportion of tips removed", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  #facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)
facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), labeller = , ncol = 3)

# observed and randomized values
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
(moran_pd_pt) + 
  (moran_pe_pt) + 
  (moran_rpd_pt) + 
  (moran_rpe_pt) +
  (moran_rand_pd_pt) + 
  (moran_rand_pe_pt) + 
  (moran_rand_rpd_pt) + 
  (moran_rand_rpe_pt) +
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")

# observed and ses values
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
(moran_pd_pt) + 
  (moran_pe_pt) + 
  (moran_rpd_pt) + 
  (moran_rpe_pt) +
  (moran_pd_ses_pt) + 
  (moran_pe_ses_pt) + 
  (moran_rpd_ses_pt) + 
  (moran_rpe_ses_pt) +
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/", "fern", "/moran_mean-all_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}

