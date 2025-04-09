library(tidyverse)
library(reshape) ; library(ggthemes)

source("code/pe_samp_funs.R")

# load results using Nitta et al. 2022 dataset ####
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))

cols2sel <- c("site", 
              #"pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              #"pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_signif", "pe_signif" ,  "rpd_signif", "rpe_signif"
              )

emp_df <- emp_canape_res %>% select(all_of(cols2sel)) ; emp_df %>% dim()

spat_df <- read.csv("data/nitta-et-al/results/japan_ferns_biodiv.csv", header = T) %>%
  select(grids, lat, long) %>% mutate(grids = as.character(grids))

# get japans layout
world <- map_data("world")
japan <- world[world$region == "Japan", ]

ggplot(data = japan, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightgray", color = "darkgray") +
  coord_sf(xlim = c(120, 150), ylim = c(20, 50)) +
  theme_map() +
  labs(title = "Japan's Coastline")

# full empirical ####
## PD prop significance ####
emp_pd_sig_df <-  # compute stats per grid
  left_join(
    emp_df %>% 
      select(site, pd_signif, pe_signif,  rpd_signif, rpe_signif) %>%
      mutate(pd_sign_cat = signif2cats(pd_signif),
             pe_sign_cat = signif2cats(pe_signif),
             rpd_sign_cat = signif2cats(rpd_signif),
             rpe_sign_cat = signif2cats(rpe_signif)), 
    spat_df, 
    by = join_by(site == grids)) %>% 
  select(site, lat, long, pd_sign_cat, pe_sign_cat, rpd_sign_cat, rpe_sign_cat) %>% 
  gather(metric, sign_cat, -site, -lat, -long) %>% 
  mutate(sign_cat = case_match(sign_cat, 
                               "non_sig" ~ "Non-sig", 
                               "sig_high" ~ "Sig. High", 
                               "sig_low" ~ "Sig. Low"),
         metric = case_match(metric,
                             "pd_sign_cat" ~ "PD",
                             "pe_sign_cat" ~ "PE",
                             "rpd_sign_cat" ~ "RPD",
                             "rpe_sign_cat" ~ "RPE"))

emp_pd_sig_df %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .75) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Empirical full", subtitle = "", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  facet_wrap(~metric, ncol = 2)
ggsave(filename = paste0("figures/map_sig_emp.pdf"), 
       width = 25, height = 25, units = "cm", dpi = 300)

# pruned #### 
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
        samp_frac = rand_frc[i], treat = "rand-pr") %>% 
    bind_rows()
}
rnd_pr_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_pr_canape_all_df)

# lateral
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

# marginal
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

prun_res_all <- rbind(rnd_pr_canape_all_df, late_pr_canape_all_df, marg_pr_canape_all_df)

## PD prop significance ####
pr_pd_sig_df <-  # compute stats per grid
  left_join( 
    prun_res_all %>%
      mutate(sign_cat = signif2cats(pd_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, pd_signif) %>%
      mutate(sign_cat_emp = ifelse(pd_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
pr_pd_sig_spat_df <- left_join(pr_pd_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_pd_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
        ) + 
  labs(title = "Pruned trees", subtitle = "Phylogenetic Diversity (PD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_pr_pd.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)

## PE prop significance ####
pr_pe_sig_df <- 
  left_join(
    prun_res_all %>%
      mutate(sign_cat = signif2cats(pe_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
  emp_df %>%
  select(site, pe_signif) %>%
  mutate(sign_cat_emp = ifelse(pe_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
  by = join_by(site)
  )

#  dataset to plot
pr_pe_sig_spat_df <- left_join(pr_pe_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_pe_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Pruned trees", subtitle = "Phylogenetic Endemism (PE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_pr_pe.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)
        
# randomly imputed #### 
cols2sel <- c("site", 
              #"pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              #"pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_signif", "pe_signif" ,  "rpd_signif", "rpe_signif"
)

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
rnd_pr_rnd_imp_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_pr_rnd_imp_canape_all_df)

# lateral
late_canape_df_ls <- list()
for(i in 1:length(late_canaper_fls)){
  late_canape_ls <- readRDS(late_canaper_fls[i])
  late_canape_df_ls[[i]] <- late_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canaper_fls))), 
        samp_frac = late_frc[i], treat = "late-pr", imputation = "random") %>% 
    bind_rows()
}
late_pr_rnd_imp__canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_pr_rnd_imp__canape_all_df)

# marginal
marg_canape_df_ls <- list()
for(i in 1:length(marg_canaper_fls)){
  marg_canape_ls <- readRDS(marg_canaper_fls[i])
  marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canaper_fls))), 
        samp_frac = marg_frc[i], treat = "marg-pr", imputation = "random") %>% 
    bind_rows()
}
marg_pr_rand_imp_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_pr_rand_imp_canape_all_df)

prun_rand_imp_all <- rbind(rnd_pr_rnd_imp_canape_all_df, late_pr_rnd_imp__canape_all_df, marg_pr_rand_imp_canape_all_df)

# compute stats per grid
## PD prop significance ####
pr_pd_sig_df <- 
  left_join(
    prun_rand_imp_all %>%
      mutate(sign_cat = signif2cats(pd_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, pd_signif) %>%
      mutate(sign_cat_emp = ifelse(pd_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
pr_pd_sig_spat_df <- left_join(pr_pd_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_pd_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Random imputation", subtitle = "Phylogenetic Diversity (PD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_rnd_imp_pd.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)

## PE prop significance ####
pr_pe_sig_df <- 
  left_join(
    prun_rand_imp_all %>%
      mutate(sign_cat = signif2cats(pe_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, pe_signif) %>%
      mutate(sign_cat_emp = ifelse(pe_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
pr_pe_sig_spat_df <- left_join(pr_pe_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_pe_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Random imputation", subtitle = "Phylogenetic Endemism (PE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_rnd_imp_pe.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)

## RPD prop significance ####
pr_rpd_sig_df <- 
  left_join(
    prun_rand_imp_all %>%
      mutate(sign_cat = signif2cats(rpd_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, rpd_signif) %>%
      mutate(sign_cat_emp = ifelse(rpd_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
pr_rpd_sig_spat_df <- left_join(pr_rpd_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_rpd_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Random imputation", subtitle = "Relative Phylogenetic Diversity (RPD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_rnd_imp_rpd.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)

## RPE prop significance ####
pr_rpe_sig_df <- 
  left_join(
    prun_rand_imp_all %>%
      mutate(sign_cat = signif2cats(rpe_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, rpe_signif) %>%
      mutate(sign_cat_emp = ifelse(rpe_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
pr_rpe_sig_spat_df <- left_join(pr_rpe_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_rpe_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Random imputation", subtitle = "Relative Phylogenetic Endemism (RPE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_rnd_imp_rpe.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)

# taxonomically imputed #### 
cols2sel <- c("site", 
              #"pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              #"pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_signif", "pe_signif" ,  "rpd_signif", "rpe_signif"
)

rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_tax", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_tax", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_tax", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

# random 
rnd_canape_df_ls <- list()
for(i in 1:length(rnd_canaper_fls)){
  rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
  rnd_canape_df_ls[[i]] <- rnd_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls))), 
        samp_frac = rand_frc[i], treat = "rand-pr", imputation = "taxonomic") %>% 
    bind_rows()
}
rnd_pr_tax_imp_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_pr_rnd_imp_canape_all_df)

# lateral
late_canape_df_ls <- list()
for(i in 1:length(late_canaper_fls)){
  late_canape_ls <- readRDS(late_canaper_fls[i])
  late_canape_df_ls[[i]] <- late_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canaper_fls))), 
        samp_frac = late_frc[i], treat = "late-pr", imputation = "taxonomic") %>% 
    bind_rows()
}
late_pr_tax_imp_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_pr_rnd_imp__canape_all_df)

# marginal
marg_canape_df_ls <- list()
for(i in 1:length(marg_canaper_fls)){
  marg_canape_ls <- readRDS(marg_canaper_fls[i])
  marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canaper_fls))), 
        samp_frac = marg_frc[i], treat = "marg-pr", imputation = "taxonomic") %>% 
    bind_rows()
}
marg_pr_tax_imp_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_pr_rand_imp_canape_all_df)

prun_tax_imp_all <- rbind(rnd_pr_tax_imp_canape_all_df, late_pr_tax_imp_canape_all_df, marg_pr_tax_imp_canape_all_df)

# compute stats per grid
## PD prop significance ####
pr_pd_sig_df <- 
  left_join(
    prun_tax_imp_all %>%
      mutate(sign_cat = signif2cats(pd_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, pd_signif) %>%
      mutate(sign_cat_emp = ifelse(pd_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
pr_pd_sig_spat_df <- left_join(pr_pd_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_pd_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Taxonomic imputation", subtitle = "Phylogenetic Diversity (PD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_pd.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)

## PE prop significance ####
pr_pe_sig_df <- 
  left_join(
    prun_tax_imp_all %>%
      mutate(sign_cat = signif2cats(pe_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, pe_signif) %>%
      mutate(sign_cat_emp = ifelse(pe_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
pr_pe_sig_spat_df <- left_join(pr_pe_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_pe_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Taxonomic imputation", subtitle = "Phylogenetic Endemism (PE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_pe.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)

## RPD prop significance ####
pr_rpd_sig_df <- 
  left_join(
    prun_tax_imp_all %>%
      mutate(sign_cat = signif2cats(rpd_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, rpd_signif) %>%
      mutate(sign_cat_emp = ifelse(rpd_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
pr_rpd_sig_spat_df <- left_join(pr_rpd_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_rpd_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Taxonomic imputation", subtitle = "Relative Phylogenetic Diversity (RPD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_rpd.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)

## RPE prop significance ####
pr_rpe_sig_df <- 
  left_join(
    prun_tax_imp_all %>%
      mutate(sign_cat = signif2cats(rpe_signif)) %>% 
      group_by(site, treat, samp_frac, sign_cat) %>% 
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>%
      spread(sign_cat, freq, fill = 0) %>% 
      rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, rpe_signif) %>%
      mutate(sign_cat_emp = ifelse(rpe_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
pr_rpe_sig_spat_df <- left_join(pr_rpe_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Marginal",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
pr_rpe_sig_spat_df %>% 
  filter(samp_frac != .2) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + 
  theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "Taxonomic imputation", subtitle = "Relative Phylogenetic Endemism (RPE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_rpe.pdf"), 
       width = 35, height = 25, units = "cm", dpi = 300)
