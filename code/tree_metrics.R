library(tidyverse)
library(RPANDA) 
library(parallel)

# load data ####
# Nitta et al. tree
fern_tre <- read.tree("output/nitta_node-named.tre") 

# other trees
trs_meta <- read.csv("data/henaodiaz_etal_trees-master/henao_diaz_etal_trees_metadata.csv", header = T) #%>% filter(study.id  != 'Testo&Sunduel2016')
tr_paths <- list.files(paste0("data/henaodiaz_etal_trees-master"), ".txt", all.files = T, full.names = T)

trs_codes <- str_extract(gsub(pattern = "data/henaodiaz_etal_trees-master/", replacement = "", x = tr_paths), "(?<=_).*?(?=_)")

# read Henao-Diaz EA trees
hd_tr_ls <- lapply(tr_paths, read.tree)

# merge Nitta EA and Henao-Diaz EA trees
test_trs <- c(list(fern_tre), hd_tr_ls)

# compute tree metrics ####

# Jensen–Shannon distance between phylogenies
jsd_dist_norm1 <- JSDtree(test_trs, meth = "normal1")
jsd_dist_stand <- JSDtree(test_trs, meth = "standard")

#clusters phylogenies using hierarchical and k-medoids clustering
JSDtree_cluster(jsd_dist_norm1, alpha = 0.9, draw = T)
JSDtree_cluster(jsd_dist_stand, alpha = 0.9, draw = T)

# principal eigenvalue, skewness, and kurtosis of the spectral density profile, and the eigengap
spec_norm <- mclapply(X = test_trs, function(X) spectR(phylo = X, meth = "normal"), mc.cores = 5)
spec_stan <- mclapply(X = test_trs, function(X) spectR(phylo = X, meth = "standard"), mc.cores = 5)

# into dataframes
spect_norm_df <- purrr::map_dfr(spec_norm, ~ .x[c("principal_eigenvalue", "asymmetry", "peakedness", "eigengap")]) %>% 
  mutate(tree_name = c("N2022", trs_codes))

spect_stan_df <- purrr::map_dfr(spec_stan, ~ .x[c("principal_eigenvalue", "asymmetry", "peakedness", "eigengap")]) %>% 
  mutate(tree_name = c("N2022", trs_codes))

# density plots
spec_nor_p <- spect_norm_df %>% 
  dplyr::select(-eigengap) %>%
  gather(key, value, -tree_name) %>% 
  ggplot(aes(x = value)) + 
  #geom_density(fill = "lightgray", color = "darkgray", alpha = 0.7) + 
  geom_histogram(fill = "lightgray", color = "darkgray", alpha = 0.7, bins = 10) +
  geom_vline(
    data = spect_norm_df %>% dplyr::select(-eigengap) %>% 
      filter(tree_name == "N2022") %>%
               gather(key, value, -tree_name), 
    aes(xintercept = value), col = "coral", alpha = 0.7, lty = 2) +
  theme_minimal() +
  labs(title = "normalized MGL", x = "") +
  facet_wrap(~key, scales = "free")
spec_nor_p
ggsave(filename = paste0("figures/tre_met_nmgl.png"), 
       width = 30, height = 10, units = "cm", dpi = 300)

##
spec_sta_p <- spect_stan_df %>% 
  dplyr::select(-eigengap) %>%
  gather(key, value, -tree_name) %>% 
  ggplot(aes(x = value)) + 
  #geom_density(fill = "lightgray", color = "darkgray", alpha = 0.7) + 
  geom_histogram(fill = "lightgray", color = "darkgray", alpha = 0.7, bins = 10) +
  geom_vline(
    data = spect_stan_df %>% dplyr::select(-eigengap) %>% 
      filter(tree_name == "N2022") %>%
      gather(key, value, -tree_name), 
    aes(xintercept = value), col = "coral", alpha = 0.7, lty = 2) +
  theme_minimal() +
  labs(title = "MGL", x = "") +
  facet_wrap(~key, scales = "free")
spec_sta_p
ggsave(filename = paste0("figures/tre_met_mgl.png"), 
       width = 30, height = 10, units = "cm", dpi = 300)

# panel with both
ggpubr::ggarrange(spec_nor_p, spec_sta_p, nrow = 2)
ggsave(filename = paste0("figures/tre_met_mgl_panel.png"), 
       width = 30, height = 15, units = "cm", dpi = 300)
