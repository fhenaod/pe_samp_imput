library(tidyverse)

# Nitta et al. spatial data
spat_df <- read.csv("data/nitta-et-al/results/japan_ferns_biodiv.csv", header = T) %>%
   mutate(site = as.character(grids)) %>% dplyr::select(site, lat, long)

# load and get canaper tables ####
# empirical full data
cols2sel <- c("site", 
              "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              "pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean")

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
  summarise(across(pd_obs:rpe_rand_mean, 
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
  summarise(across(pd_obs:rpe_rand_mean, ~mean(.x, na.rm = T), .names = "m_{.col}"))

canaper_mean_spat_df <- left_join(canaper_mean_res_df %>% ungroup(), spat_df, by = join_by("site"))

# get mean moran's I
run_moran = F
if(run_moran){
autocor_mean_df <- canaper_mean_spat_df %>% 
  filter(!is.na(long)) %>% 
  group_by(samp_frac, treat, imputation) %>%
  summarise(across(m_pd_obs:m_rpe_rand_mean, 
                ~ape::Moran.I(.x, weight = as.matrix(1/dist(cbind(long, lat))), na.rm = T)$observed, 
                .names = "moran_{.col}"))
autocor_mean_df %>% saveRDS("output/autocor_mean_df.rds") ; autocor_mean_df %>% head()
} else {
  autocor_mean_df <- readRDS("output/autocor_mean_df.rds")
}

# figures #####
library(patchwork)
forms2print <- c(".png", ".pdf", ".svg") 

# mean and sd vals ####

# mean value across sims ####
# one panel
autocor_mean_df %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
          cat2map = paste0(treat, "-" , imputation)) %>% 
  gather(key, value, -samp_frac, -treat, -imputation, -is_imp, -cat2map) %>%
  ggplot(aes(x = (samp_frac), y = value, color = treat, shape = imputation)) +
  geom_point() + 
  scale_shape_manual(values = c(16, 8, 17, 15)) +
  geom_line(aes(linetype = imputation)) + 
  labs(#title = "Phylogenetic Diversity (PD)", 
    x = "Proportion of tips removed", y = "Moran's I", color = "Pruning\nmethod", shape = "Imputation", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), option = "inferno", begin = .8, end = .05, alpha = .95) + 
  facet_wrap(~key*is_imp, labeller = , ncol = 2, scales = "free")

# plots by metric
p_size = 3.5

# PD
moran_pd_pt <- autocor_mean_df %>%
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_pd_obs, color = imputation, shape = treat)) +
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pd_obs), colour = "black", pch = 18, size = p_size) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Phylogenetic Diversity (PD)", 
       x = "", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# PE
moran_pe_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_pe_obs, color = imputation, shape = treat)) +
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pe_obs), colour = "black", pch = 18, size = p_size) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Phylogenetic Endemism (PE)", 
       x = "", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPD
moran_rpd_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_rpd_obs, color = imputation, shape = treat)) +
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpd_obs), colour = "black", pch = 18, size = p_size) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Relative Phylogenetic Diversity (RPD)", 
       x = "", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# RPE
moran_rpe_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_rpe_obs, color = imputation, shape = treat)) +
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpe_obs), colour = "black", pch = 18, size = p_size) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Relative Phylogenetic Endemism (RPE)", 
       x = "Proportion of tips removed", y = "Moran's I", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

# randPD
moran_rand_pd_pt <- autocor_mean_df %>% 
  #filter(imputation != "random") %>% 
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
         cat2map = paste0(treat, "-" , imputation)) %>%
  ggplot(aes(x = samp_frac, y = moran_m_pd_rand_mean, color = imputation, shape = treat)) +
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pd_rand_mean), colour = "black", pch = 18, size = p_size) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Randomized Mean Phylogenetic Diversity (randPD)", 
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
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_pe_rand_mean), colour = "black", pch = 18, size = p_size) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Randomized Mean Phylogenetic Endemism (randPE)", 
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
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpd_rand_mean), colour = "black", pch = 18, size = p_size) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Randomized Mean Relative Phylogenetic Diversity (randRPD)", 
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
  geom_point() + 
  geom_point(aes(x = 0, y = moran_emp$moran_rpe_rand_mean), colour = "black", pch = 18, size = p_size) +
  geom_line(aes(linetype = imputation)) + 
  #scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  labs(subtitle = "Randomized Mean Relative Phylogenetic Endemism (randRPE)", 
       x = "Proportion of tips removed", y = "", color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "right") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  facet_wrap(~factor(is_imp, c("Only pruned", "Imputed")), labeller = , ncol = 2)

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

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  ggsave(filename = paste0("figures/moran_mean-all_pnl", fig_out), 
         width = 30, height = 40, units = "cm", dpi = 200)
}
