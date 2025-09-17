library(tidyverse)
library(phangorn)
library(TreeDist)
library(progress)

source("code/pe_samp_funs.R")

# load tree
fern_tre <- read.tree("output/nitta_node-named.tre") 
ape::is.binary(fern_tre) ; ape::is.ultrametric(fern_tre)
targ_tre <- fern_tre

# pruned trees
pr_trs_p <- dir("output/pruned_trees/", full.names = T)

# imputed trees
im_trs_p <- dir("output/imputed_trees/", full.names = T)
#im_trs_p <- im_trs_p[-which(grepl("_rand_imp_", im_trs_p),im_trs_p)]


# distance based - phangorn ####
run <- F
if(run){
n_imp_ts <- list()
for(i in 1:length(im_trs_p)){
  temp_trs <- readRDS(im_trs_p[i])
  n_temp_trs <- temp_trs %>% lengths()
  
  if(length(temp_trs) > 5){
    n_imp_ts[[i]] <- sapply(X = temp_trs, 
                            function(X) treedist(targ_tre, X, check.labels = T), simplify = T) %>% 
      t() %>% data.frame(., im_trs_p[i])
                            } else {
    n_imp_ts[[i]] <- treedist(targ_tre, temp_trs, check.labels = T) %>% data.frame(., im_trs_p[i])
                            }
  }
  tree_m_df <- n_imp_ts %>% bind_rows()
  saveRDS(tree_m_df, "output/phangor_metrics_df.rds")
  system("say Just finished!") 
  } else {
  tree_m_df <- readRDS("output/phangor_metrics_df.rds")
  head(tree_m_df)
}

treestats_df <- tree_m_df %>% 
  mutate(im_trs_p.i. = gsub(pattern = "output/imputed_trees//", replacement = "", x = im_trs_p.i.),
         im_trs_p.i. = gsub(pattern = ".rds", replacement = "", x = im_trs_p.i.),
         treat = str_extract(im_trs_p.i., ".*_pr"),
         samp_fract = parse_number(im_trs_p.i.),
         imputation = str_extract(im_trs_p.i., "(?<=pr_)(.*?)(?=[0-9])"))
  
treestats_df %>% 
  gather(key = stat, value, -c(im_trs_p.i.:imputation)) %>%
  group_by(treat, samp_fract, imputation, stat) %>% 
  summarize(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) %>% 
  #filter(imputation == 'rand_imp_') %>% 
  ggplot(aes(x = (samp_fract), y = mean, ymin = mean - sd, ymax = mean + sd, color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  #scale_shape_manual(labels = c("Lateral", "Lateral - Spatial", "Bilateral", "Bilateral - Spatial", "Random"), values = c(16, 17, 15, 18, 8)) +
  geom_line(aes(linetype = imputation)) + 
  labs(title = "Tree comparison metrics", 
       x = "Proportion of tips removed", y = "", 
      # color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation"
      ) +
  theme_classic(base_size = 12) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  #guides(shape = guide_legend(override.aes = list(linetype = c(0, 0, 0, 0, 0), shape = c(16, 17, 15, 18, 8))), 
  #       colour = guide_legend(override.aes = list(linetype = c(1, 2 ,3), shape = c(NA, NA, NA)))) +
  theme(strip.background = element_blank(), legend.position = "right") + 
  facet_wrap(~stat, nrow = 2, scales = "free")
ggsave(filename = paste0("figures/tree_metrs_1pnl", ".png"), 
       width = 20, height = 20, units = "cm", dpi = 200)

# TreeDist ####
## tree distance ####
run <- F
if(run){
  n_iterations <- length(im_trs_p)
  pb <- progress_bar$new(total = n_iterations, format = "[:bar] :percent :eta")
  
  nn_imp_ts <- list()
  for(i in 1:length(im_trs_p)){
    temp_trs <- readRDS(im_trs_p[i])
    nn_temp_trs <- temp_trs %>% lengths()
    
    if(length(nn_temp_trs) > 5){
      nn_imp_ts[[i]] <- sapply(X = temp_trs, 
                              function(X) get_t2t_metrics(ref_tree = targ_tre, compar_tree = X, sim = T, sum = T, norm = T), simplify = T) %>% 
        t() %>% data.frame(., im_trs_p[i])
    } else {
      nn_imp_ts[[i]] <- get_t2t_metrics(ref_tree = targ_tre, compar_tree = temp_trs, sim = T, sum = T, norm = T) %>% data.frame(., im_trs_p[i])
    }
    pb$tick()
  }
  
  treedist_m_df <- nn_imp_ts %>% do.call(rbind, .) %>%
    mutate(across(.cols = rob.fou_norm:shared.phy_info, .fns = unlist))
  saveRDS(treedist_m_df, "output/treedist_metrics_df.rds")
  system("say Just finished!")
} else {
  treedist_m_df <- readRDS("output/treedist_metrics_df.rds")
  head(treedist_m_df)
}

get_t2t_metrics(ref_tree = targ_tre, compar_tree = (targ_tre), sim = T, sum = T, norm = T)

# clean data frame
treestats_df <- treedist_m_df %>% 
  mutate(im_trs_p.i. = gsub(pattern = "output/imputed_trees//", replacement = "", x = im_trs_p.i.),
         im_trs_p.i. = gsub(pattern = ".rds", replacement = "", x = im_trs_p.i.),
         treat = str_extract(im_trs_p.i., ".*_pr"),
         samp_fract = parse_number(im_trs_p.i.),
         imputation = str_extract(im_trs_p.i., "(?<=pr_)(.*?)(?=[0-9])"),
         across(rob.fou_norm:phylo_info_dist, ~as.numeric(.x)))

# get mean and sd vals per group
treestats_df %>% 
  gather(key = stat, value, -c(im_trs_p.i.:imputation)) %>% 
  group_by(treat, samp_fract, imputation, stat) %>% 
  summarize(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) %>% 
  ggplot(aes(x = (samp_fract), y = mean, ymin = mean - sd, ymax = mean + sd, 
             color = imputation, shape = treat
             #color = treat, shape = imputation
             )) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  geom_point(alpha = .5) + #geom_smooth(se = F) +
  geom_line(aes(linetype = imputation)) + 
  labs(title = "Tree comparison metrics", 
       x = "Proportion of tips removed", y = "Similarity", 
       #color = "Imputation", shape = "Pruning\nmethod", linetype = "Imputation"
       ) +
  theme_classic(base_size = 12) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  theme(strip.background = element_blank(), legend.position = "right") + 
  facet_wrap(~stat, nrow = 1, scales = "free")
ggsave(filename = paste0("figures/treedist_metrs_1pnl", ".png"), 
       width = 20, height = 20, units = "cm", dpi = 200)
  

## tree info ####
all_trees_path <- c(pr_trs_p, im_trs_p)

run <- F
if(run){
  n_iterations <- length(all_trees_path)
  pb <- progress_bar$new(total = n_iterations, format = "[:bar] :percent :eta")
  
  nn_imp_ts <- list()
  for(i in 1:n_iterations){
    temp_trs <- readRDS(all_trees_path[i])
    nn_temp_trs <- temp_trs %>% lengths()
    
    if(length(nn_temp_trs) > 5){
      nn_imp_ts[[i]] <- sapply(X = temp_trs, function(X) get_tree_metrics(tree = X,  sum = T), simplify = T) %>% 
        t() %>% data.frame(., all_trees_path[i])
    } else {
      nn_imp_ts[[i]] <- data.frame(get_tree_metrics(tree = temp_trs, sum = T), all_trees_path[i])
    }
    pb$tick()
  }
  
  tree_metrics_df <- nn_imp_ts %>% do.call(rbind, .) %>%
    mutate(across(.cols = splitwise_info:clustering_info, .fns = unlist))
  saveRDS(treedist_m_df, "output/tree_info_metrics_df.rds")
  system("say Just finished!")
} else {
  tree_metrics_df <- readRDS("output/tree_info_metrics_df.rds")
  head(treedist_m_df)
}

get_tree_metrics(tree = targ_tre, sum = T)

# clean data frame
tree_info_df <- tree_metrics_df %>% 
  mutate(all_trees_path.i. = gsub(pattern = "output/", replacement = "", x = all_trees_path.i.),
         all_trees_path.i. = gsub(pattern = "pruned_trees//", replacement = "", x = all_trees_path.i.),
         all_trees_path.i. = gsub(pattern = "imputed_trees//", replacement = "", x = all_trees_path.i.),
         all_trees_path.i. = gsub(pattern = ".rds", replacement = "", x = all_trees_path.i.),
         treat = str_extract(all_trees_path.i., ".*_pr"),
         samp_fract = parse_number(all_trees_path.i.),
         imputation = str_extract(all_trees_path.i., "(?<=pr_)(.*?)(?=[0-9])"),
         imputation = replace_na(imputation, "none"))

tree_info_df %>% 
  gather(key = stat, value, -c(all_trees_path.i.:imputation)) %>% 
  group_by(treat, samp_fract, imputation, stat) %>% 
  summarize(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) %>% 
  ggplot(aes(x = (samp_fract), y = (mean), ymin = mean - sd, ymax = mean + sd, 
             color = imputation, shape = treat)) +
  geom_pointrange(position = position_dodge2(width = .03, padding = 0.15)) + 
  #geom_point() + #geom_smooth(se = F) +
  geom_line(aes(linetype = imputation)) + 
  labs(title = "Tree comparison metrics", 
       x = "Proportion of tips removed", y = "", 
       #color = "Pruning\nmethod", shape = "Pruning\nmethod", linetype = "Pruning\nmethod"
       ) +
  theme_classic(base_size = 12) + 
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  theme(strip.background = element_blank(), legend.position = "right") + 
  facet_wrap(~stat, nrow = 1, scales = "free")
