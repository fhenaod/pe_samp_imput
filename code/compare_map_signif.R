library(tidyverse)
library(reshape) ; library(ggthemes)

source("code/pe_samp_funs.R")

# load results using Nitta et al. 2022 dataset ####
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))

cols2sel <- c("site", 
              #"pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              #"pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_signif", "pe_signif", "rpd_signif", "rpe_signif")

emp_df <- emp_canape_res %>% select(all_of(cols2sel)) ; emp_df %>% dim()

spat_df <- read.csv("data/nitta-et-al/results/japan_ferns_biodiv.csv", header = T) %>%
  select(grids, lat, long) %>% mutate(grids = as.character(grids))

# get japan's layout
world <- map_data("world")
japan <- world[world$region == "Japan", ]

ggplot(data = japan, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightgray", color = "darkgray") +
  coord_sf(xlim = c(120, 150), ylim = c(20, 50)) +
  theme_map() +
  labs(title = "Japan's Coastline")

frt_out <- ".png" #".pdf"
pos_do <- position_dodge(0.1)

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
ggsave(filename = paste0("figures/map_sig_emp", frt_out), 
       width = 25, height = 25, units = "cm", dpi = 300)

# pruned #### 
cols2sel <- c("site", 
              #"pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              #"pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_signif", "pe_signif"#, "rpd_signif", "rpe_signif"
              )

names2load <- list.dirs(path = "output/canaper/", full.names = F, recursive = F)
print(names2load)

to_load <- "fern"

rnd_canaper_fls <- list.files(paste0("output/canaper/", to_load), "canaper_rand_pr_[0-9]", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files(paste0("output/canaper/", to_load), "canaper_late_pr_[0-9]", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files(paste0("output/canaper/", to_load), "canaper_marg_pr_[0-9]", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

get_tabs2res <- T
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
prun_res_all %>% saveRDS(file = paste0("output/canaper/canaper_prun_res_all_", to_load ,".rds"))
} else {
  prun_res_all <- readRDS(paste0("output/canaper/canaper_prun_res_all_", to_load, ".rds")) ; head(prun_res_all)
}

## PD prop significance ####
# compute stats per grid
pr_pd_sig_df <-  
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

# dataset to plot
pr_pd_sig_spat_df <- left_join(pr_pd_sig_df, spat_df, 
                               by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
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
ggsave(filename = paste0("figures/map_sig_pr_pd", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
pr_pd_res_sig_cats <- left_join(
  prun_res_all %>%
    mutate(sign_cat = signif2cats(pd_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, pd_signif) %>% 
    mutate(sign_cat_emp = ifelse(pd_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(pd_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(pd_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

pd_res_type_pr <- pr_pd_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") +
  ylim(0, 1) +
  labs(title = "Phylogenetic Diversity (PD) significance", 
       subtitle = "Pruned-only",
       x = "", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
pd_res_type_pr

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
                            "marg-pr" ~ "Bilateral",
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
ggsave(filename = paste0("figures/map_sig_pr_pe", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)
    
# grid prop significance by result type
pr_pe_res_sig_cats <- left_join(
  prun_res_all %>%
    mutate(sign_cat = signif2cats(pe_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, pe_signif) %>% 
    mutate(sign_cat_emp = ifelse(pe_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(pe_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(pe_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

pe_res_type_pr <- pr_pe_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(title = "Phylogenetic Endemism (PE) significance", 
       subtitle = "Pruned-only",
       x = "", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
pe_res_type_pr

# randomly imputed #### 
cols2sel <- c("site", 
              #"pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              #"pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_signif", "pe_signif", "rpd_signif", "rpe_signif"
)

rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_rand", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_rand", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_rand", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

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
  late_pr_rnd_imp_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_pr_rnd_imp_canape_all_df)

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

  rand_imp_all <- rbind(rnd_pr_rnd_imp_canape_all_df, late_pr_rnd_imp_canape_all_df, marg_pr_rand_imp_canape_all_df)
  rand_imp_all %>% saveRDS("output/prun_rand_imp_all.rds")
} else {
  rand_imp_all <- readRDS("output/prun_rand_imp_all.rds"); head(rand_imp_all)
}

# compute stats per grid
## PD prop significance ####
rand_imp_pd_sig_df <- 
  left_join(
    rand_imp_all %>%
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
rand_imp_pd_sig_spat_df <- left_join(rand_imp_pd_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
rand_imp_pd_sig_spat_df %>% 
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
ggsave(filename = paste0("figures/map_sig_rnd_imp_pd", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
rand_imp_pd_res_sig_cats <- left_join(
  rand_imp_all %>%
    mutate(sign_cat = signif2cats(pd_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, pd_signif) %>% 
    mutate(sign_cat_emp = ifelse(pd_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(pd_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(pd_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

pd_res_type_rand_imp <- rand_imp_pd_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Phylogenetic Diversity (PD) significance",
       subtitle = "Imputed: random",
       x = "", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
pd_res_type_rand_imp

## PE prop significance ####
rand_imp_pe_sig_df <- 
  left_join(
    rand_imp_all %>%
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
rand_imp_pe_sig_spat_df <- left_join(rand_imp_pe_sig_df, spat_df, 
                                     by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
rand_imp_pe_sig_spat_df %>% 
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
ggsave(filename = paste0("figures/map_sig_rnd_imp_pe", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
rand_imp_pe_res_sig_cats <- left_join(
  rand_imp_all %>%
    mutate(sign_cat = signif2cats(pe_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, pe_signif) %>% 
    mutate(sign_cat_emp = ifelse(pe_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(pe_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(pe_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
  )

pe_res_type_rand_imp <- rand_imp_pe_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Phylogenetic Endemism (PE) significance",
       subtitle = "Imputed: random",
       x = "Proportion of tips pruned", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
pe_res_type_rand_imp

## RPD prop significance ####
rand_imp_rpd_sig_df <- 
  left_join(
    rand_imp_all %>% 
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
rand_imp_rpd_sig_spat_df <- left_join(rand_imp_rpd_sig_df, spat_df, 
                                      by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
rand_imp_rpd_sig_spat_df %>% 
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
  labs(title = "Random imputation", 
       subtitle = "Relative Phylogenetic Diversity (RPD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_rnd_imp_rpd", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
rand_imp_rpd_res_sig_cats <- left_join(
  rand_imp_all %>%
    mutate(sign_cat = signif2cats(rpd_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, rpd_signif) %>% 
    mutate(sign_cat_emp = ifelse(rpd_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(rpd_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(rpd_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

rpd_res_type_rand_imp <- rand_imp_rpd_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(title = "Relative Phylogenetic Diversity (RPD) significance",
       subtitle = "Imputed: random",
       x = "", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
rpd_res_type_rand_imp

## RPE prop significance ####
rand_rpe_sig_df <- 
  left_join(
    rand_imp_all %>%
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
rand_imp_rpe_sig_spat_df <- left_join(rand_rpe_sig_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
rand_imp_rpe_sig_spat_df %>% 
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
  labs(title = "Random imputation", 
       subtitle = "Relative Phylogenetic Endemism (RPE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_rnd_imp_rpe", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
rand_imp_rpe_res_sig_cats <- left_join(
  rand_imp_all %>%
    mutate(sign_cat = signif2cats(rpe_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, rpe_signif) %>% 
    mutate(sign_cat_emp = ifelse(rpe_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(rpe_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(rpe_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

rpe_res_type_rand_imp <- rand_imp_rpe_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(title = "Relative Phylogenetic Endemism (RPE) significance",
       subtitle = "Imputed: random",
       x = "", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
rpe_res_type_rand_imp

# taxonomic imputed #### 
## rank's node ####
cols2sel <- c("site", 
              #"pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              #"pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_signif", "pe_signif" ,  "rpd_signif", "rpe_signif"
)

rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_tax_imp_root", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_tax_imp_root", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_tax_imp_root", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

get_tabs2res <- F
if(get_tabs2res){
# random 
  rnd_canape_df_ls <- list()
  for(i in 1:length(rnd_canaper_fls)){
  rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
  rnd_canape_df_ls[[i]] <- rnd_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls))), 
        samp_frac = rand_frc[i], treat = "rand-pr", imputation = "taxonomic - node") %>% 
    bind_rows()
}
  rnd_pr_tax_imp_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_pr_tax_imp_canape_all_df) %>% print()

# lateral
  late_canape_df_ls <- list()
  for(i in 1:length(late_canaper_fls)){
  late_canape_ls <- readRDS(late_canaper_fls[i])
  late_canape_df_ls[[i]] <- late_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canaper_fls))), 
        samp_frac = late_frc[i], treat = "late-pr", imputation = "taxonomic - node") %>% 
    bind_rows()
}
  late_pr_tax_imp_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_pr_tax_imp_canape_all_df) %>% print()

# marginal
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_canaper_fls)){
  marg_canape_ls <- readRDS(marg_canaper_fls[i])
  marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canaper_fls))), 
        samp_frac = marg_frc[i], treat = "marg-pr", imputation = "taxonomic - node") %>% 
    bind_rows()
}
  marg_pr_tax_imp_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_pr_tax_imp_canape_all_df) %>% print()

  prun_tax_imp_node_all <- rbind(rnd_pr_tax_imp_canape_all_df, late_pr_tax_imp_canape_all_df, marg_pr_tax_imp_canape_all_df)
  prun_tax_imp_node_all %>% saveRDS("output/prun_tax_imp-node_all.rds")
} else {
  prun_tax_imp_node_all <- readRDS("output/prun_tax_imp-node_all.rds") ; head(prun_tax_imp_node_all)
}

# compute stats per grid
### PD prop significance ####
root_imp_pd_sig_df <- 
  left_join(
    prun_tax_imp_node_all %>%
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
root_imp_pd_sig_spat_df <- left_join(root_imp_pd_sig_df, spat_df, 
                                     by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral", 
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"), 
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
root_imp_pd_sig_spat_df %>% 
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
  labs(title = "Taxonomic imputation - node", 
       subtitle = "Phylogenetic Diversity (PD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_node_pd", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
root_imp_pd_res_sig_cats <- left_join(
  prun_tax_imp_node_all %>%
    mutate(sign_cat = signif2cats(pd_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, pd_signif) %>% 
    mutate(sign_cat_emp = ifelse(pd_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(pd_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(pd_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

pd_res_type_root_imp <- root_imp_pd_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Phylogenetic Diversity (PD) significance",
       subtitle = "Imputed: taxonomy - node",
       x = "Proportion of tips pruned", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
pd_res_type_root_imp

### PE prop significance ####
root_imp_pe_sig_df <- 
  left_join(
    prun_tax_imp_node_all %>%
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
root_imp_pe_sig_spat_df <- left_join(root_imp_pe_sig_df, spat_df, 
                                     by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
root_imp_pe_sig_spat_df %>% 
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
  labs(title = "Taxonomic imputation - node", 
       subtitle = "Phylogenetic Endemism (PE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_node_pe", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
root_imp_pe_res_sig_cats <- left_join(prun_tax_imp_node_all %>% 
                                        mutate(sign_cat = signif2cats(pe_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, pe_signif) %>% 
    mutate(sign_cat_emp = ifelse(pe_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(pe_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(pe_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

pe_res_type_root_imp <- root_imp_pe_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) +
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Phylogenetic Endemism (PE) significance",
       subtitle = "Imputed: taxonomy - node",
       x = "Proportion of tips pruned", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
pe_res_type_root_imp

### RPD prop significance ####
root_imp_rpd_sig_df <- left_join(prun_tax_imp_node_all %>% 
                                   mutate(sign_cat = signif2cats(rpd_signif)) %>% 
                                   group_by(site, treat, samp_frac, sign_cat) %>% 
                                   summarise(n = n()) %>% mutate(freq = n / sum(n)) %>% 
                                   select(-n) %>% group_by(samp_frac) %>% 
                                   spread(sign_cat, freq, fill = 0) %>% rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, rpd_signif) %>%
      mutate(sign_cat_emp = ifelse(rpd_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

# dataset to plot
root_imp_rpd_sig_spat_df <- left_join(root_imp_rpd_sig_df, spat_df, 
                                      by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
root_imp_rpd_sig_spat_df %>% 
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
  labs(title = "Taxonomic imputation - node", 
       subtitle = "Relative Phylogenetic Diversity (RPD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_node_rpd", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
root_imp_rpd_res_sig_cats <- left_join(prun_tax_imp_node_all %>% 
                                         mutate(sign_cat = signif2cats(rpd_signif)) %>% 
                                         group_by(site, treat, samp_frac, sign_cat) %>% 
                                         summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, rpd_signif) %>% 
    mutate(sign_cat_emp = ifelse(rpd_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(rpd_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(rpd_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

rpd_res_type_root_imp <- root_imp_rpd_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Relative Phylogenetic Diversity (RPD) significance",
       subtitle = "Imputed: taxonomy - node", 
       x = "Proportion of tips pruned", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
rpd_res_type_root_imp

### RPE prop significance ####
root_imp_rpe_sig_df <- left_join(prun_tax_imp_node_all %>% 
                                   mutate(sign_cat = signif2cats(rpe_signif)) %>% 
                                   group_by(site, treat, samp_frac, sign_cat) %>% 
                                   summarise(n = n()) %>% mutate(freq = n / sum(n)) %>% 
                                   select(-n) %>% group_by(samp_frac) %>% 
                                   spread(sign_cat, freq, fill = 0) %>% rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, rpe_signif) %>%
      mutate(sign_cat_emp = ifelse(rpe_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
root_imp_rpe_sig_spat_df <- left_join(root_imp_rpe_sig_df, spat_df, 
                                      by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
root_imp_rpe_sig_spat_df %>% 
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
  labs(title = "Taxonomic imputation - node", 
       subtitle = "Relative Phylogenetic Endemism (RPE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_node_rpe", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
tax_imp_rpe_res_sig_cats <- left_join(prun_tax_imp_node_all %>%
                                        mutate(sign_cat = signif2cats(rpe_signif)) %>% 
                                        group_by(site, treat, samp_frac, sign_cat) %>% 
                                        summarise(n = n()) %>% mutate(freq = n / sum(n)) %>% 
                                        select(-n) %>% group_by(samp_frac) %>% 
                                        spread(sign_cat, freq, fill = 0) %>% rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, rpe_signif) %>% 
    mutate(sign_cat_emp = ifelse(rpe_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(rpe_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(rpe_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

rpe_res_type_root_imp <- tax_imp_rpe_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Relative Phylogenetic Endemism (RPE) significance",
       subtitle = "Imputed: taxonomy - node",
       x = "Proportion of tips pruned", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
rpe_res_type_root_imp

## random within subtree ####
cols2sel <- c("site", 
              #"pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              #"pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
              "pd_signif", "pe_signif" ,  "rpd_signif", "rpe_signif"
)

rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_tax_imp_random", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_tax_imp_random", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_tax_imp_random", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

get_tabs2res <- F
if(get_tabs2res){
  # random 
  rnd_canape_df_ls <- list()
  for(i in 1:length(rnd_canaper_fls)){
    rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
    rnd_canape_df_ls[[i]] <- rnd_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls))), 
          samp_frac = rand_frc[i], treat = "rand-pr", imputation = "taxonomic - random") %>% 
      bind_rows()
  }
  rnd_pr_tax_imp_canape_all_df <- data.table::rbindlist(rnd_canape_df_ls) ; dim(rnd_pr_tax_imp_canape_all_df) %>% print()
  
  # lateral
  late_canape_df_ls <- list()
  for(i in 1:length(late_canaper_fls)){
    late_canape_ls <- readRDS(late_canaper_fls[i])
    late_canape_df_ls[[i]] <- late_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canaper_fls))), 
          samp_frac = late_frc[i], treat = "late-pr", imputation = "taxonomic - random") %>% 
      bind_rows()
  }
  late_pr_tax_imp_canape_all_df <- data.table::rbindlist(late_canape_df_ls) ; dim(late_pr_tax_imp_canape_all_df) %>% print()
  
  # marginal
  marg_canape_df_ls <- list()
  for(i in 1:length(marg_canaper_fls)){
    marg_canape_ls <- readRDS(marg_canaper_fls[i])
    marg_canape_df_ls[[i]] <- marg_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canaper_fls))), 
          samp_frac = marg_frc[i], treat = "marg-pr", imputation = "taxonomic - random") %>% 
      bind_rows()
  }
  marg_pr_tax_imp_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_pr_tax_imp_canape_all_df) %>% print()
  
  prun_tax_imp_rand_all <- rbind(rnd_pr_tax_imp_canape_all_df, late_pr_tax_imp_canape_all_df, marg_pr_tax_imp_canape_all_df)
  prun_tax_imp_rand_all %>% saveRDS("output/prun_tax_imp-random_all.rds")
} else {
  prun_tax_imp_rand_all <- readRDS("output/prun_tax_imp-random_all.rds") ; head(prun_tax_imp_rand_all)
}

# compute stats per grid
### PD prop significance ####
rand_imp_pd_sig_df <- 
  left_join(prun_tax_imp_rand_all %>% 
              mutate(sign_cat = signif2cats(pd_signif)) %>% 
              group_by(site, treat, samp_frac, sign_cat) %>% 
              summarise(n = n()) %>% mutate(freq = n / sum(n)) %>% 
              select(-n) %>% group_by(samp_frac) %>% 
              spread(sign_cat, freq, fill = 0) %>% rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, pd_signif) %>%
      mutate(sign_cat_emp = ifelse(pd_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
tax_rand_imp_pd_sig_spat_df <- left_join(rand_imp_pd_sig_df, spat_df, 
                                     by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
tax_rand_imp_pd_sig_spat_df %>% 
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
  labs(title = "Taxonomic imputation - random", 
       subtitle = "Phylogenetic Diversity (PD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp-random_pd", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
tax_rand_imp_pd_res_sig_cats <- left_join(prun_tax_imp_rand_all %>% 
                                        mutate(sign_cat = signif2cats(pd_signif)) %>% 
                                        group_by(site, treat, samp_frac, sign_cat) %>% 
                                        summarise(n = n()) %>% 
                                        mutate(freq = n / sum(n)) %>% select(-n) %>% 
                                        group_by(samp_frac) %>% spread(sign_cat, freq, fill = 0) %>% 
                                        rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, pd_signif) %>% 
    mutate(sign_cat_emp = ifelse(pd_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(pd_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(pd_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

pd_res_type_tax_rand_imp <- tax_rand_imp_pd_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Phylogenetic Diversity (PD) significance",
    subtitle = "Imputed: taxonomy - random",
    x = "Proportion of tips pruned", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
pd_res_type_tax_rand_imp

### PE prop significance ####
tax_rand_imp_pe_sig_df <- left_join(prun_tax_imp_rand_all %>% 
                                  mutate(sign_cat = signif2cats(pe_signif)) %>% 
                                  group_by(site, treat, samp_frac, sign_cat) %>% 
                                  summarise(n = n()) %>% mutate(freq = n / sum(n)) %>% 
                                  select(-n) %>% group_by(samp_frac) %>% 
                                  spread(sign_cat, freq, fill = 0) %>% rowwise() %>%
      mutate(sign_cat = ifelse(non_sig >= .6, "Non-sig", 
                               ifelse(sig_high > .2, "Sig. High", 
                                      ifelse(sig_low > .2, "Sig. Low", NA)))), 
    emp_df %>%
      select(site, pe_signif) %>%
      mutate(sign_cat_emp = ifelse(pe_signif %in% c("< 0.01", "< 0.025", "> 0.975", "> 0.99"), "*", "")),
    by = join_by(site)
  )

#  dataset to plot
tax_rand_imp_pe_sig_spat_df <- left_join(tax_rand_imp_pe_sig_df, spat_df, 
                                         by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
tax_rand_imp_pe_sig_spat_df %>% 
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
  labs(title = "Taxonomic imputation - random", 
       subtitle = "Phylogenetic Endemism (PE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_rand_pe", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
tax_rand_imp_pe_res_sig_cats <- left_join(
  prun_tax_imp_rand_all %>%
    mutate(sign_cat = signif2cats(pe_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, pe_signif) %>% 
    mutate(sign_cat_emp = ifelse(pe_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(pe_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(pe_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

pe_res_type_tax_rand_imp <- tax_rand_imp_pe_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) +
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Phylogenetic Endemism (PE) significance",
    subtitle = "Imputed: taxonomy - random",
    x = "Proportion of tips pruned", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
pe_res_type_tax_rand_imp

### RPD prop significance ####
tax_rand_imp_rpd_sig_df <- 
  left_join(
    prun_tax_imp_rand_all %>%
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
tax_rand_imp_rpd_sig_spat_df <- left_join(tax_rand_imp_rpd_sig_df, spat_df, 
                                      by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
tax_rand_imp_rpd_sig_spat_df %>% 
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
  labs(title = "Taxonomic imputation - random", 
       subtitle = "Relative Phylogenetic Diversity (RPD)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_rand_rpd", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
tax_rand_imp_rpd_res_sig_cats <- left_join(
  prun_tax_imp_rand_all %>%
    mutate(sign_cat = signif2cats(rpd_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, rpd_signif) %>% 
    mutate(sign_cat_emp = ifelse(rpd_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(rpd_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(rpd_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

rpd_res_type_tax_rand_imp <- tax_rand_imp_rpd_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Relative Phylogenetic Diversity (RPD) significance",
    subtitle = "Imputed: taxonomy - random", 
    x = "Proportion of tips pruned", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
rpd_res_type_tax_rand_imp

### RPE prop significance ####
tax_rand_imp_rpe_sig_df <- 
  left_join(
    prun_tax_imp_rand_all %>%
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
tax_rand_imp_rpe_sig_spat_df <- left_join(tax_rand_imp_rpe_sig_df, spat_df, 
                                      by = join_by(site == grids)) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random"),
         treat = paste0(treat, " - ", samp_frac))

# plot all the maps
tax_rand_imp_rpe_sig_spat_df %>% 
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
  labs(title = "Taxonomic imputation - random", 
       subtitle = "Relative Phylogenetic Endemism (RPE)", 
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) +
  facet_wrap(~treat, ncol = 4)
ggsave(filename = paste0("figures/map_sig_tax_imp_rand_rpe", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# grid prop significance by result type
tax_rand_imp_rpe_res_sig_cats <- left_join(
  prun_tax_imp_rand_all %>%
    mutate(sign_cat = signif2cats(rpe_signif)) %>% 
    group_by(site, treat, samp_frac, sign_cat) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>% select(-n) %>% group_by(samp_frac) %>% 
    spread(sign_cat, freq, fill = 0) %>% 
    rowwise() %>% 
    mutate(sign_cat = ifelse(non_sig >= .9, "non_sig", 
                             ifelse(sig_high >= .9, "sig_high", 
                                    ifelse(sig_low >= .9, "sig_low", "n_a")))),
  emp_df %>%
    select(site, rpe_signif) %>% 
    mutate(sign_cat_emp = ifelse(rpe_signif %in% c("< 0.01", "< 0.025"), "sig_low", 
                                 ifelse(rpe_signif %in% c("> 0.975", "> 0.99"), "sig_high", 
                                        ifelse(rpe_signif == "not significant", "non_sig", "n_a")))), 
  by = join_by(site)
)

rpe_res_type_tax_rand_imp <- tax_rand_imp_rpe_res_sig_cats %>% 
  mutate(res_type = res_type(sign_cat, sign_cat_emp),
         treat = case_match(treat, 
                            "late-pr" ~ "Lateral",
                            "marg-pr" ~ "Bilateral",
                            "rand-pr" ~ "Random")) %>% #filter(is.na(res_type))
  group_by(treat, samp_frac, res_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = samp_frac, y = freq, group = res_type, col = res_type)) +
  geom_point(cex = .85, position = pos_do) + 
  geom_line(position = pos_do) +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        legend.position = "right") + ylim(0, 1) +
  labs(#title = "Relative Phylogenetic Endemism (RPE) significance",
    subtitle = "Imputed: taxonomy - random",
    x = "Proportion of tips pruned", y = "Frequency of grid cells", color = "") +
  facet_wrap(~ treat, scales = "fixed")
rpe_res_type_tax_rand_imp

# map figs ####
## PD ####
# full data
map_emp_pd <- emp_pd_sig_df %>% 
  filter(metric == "PD") %>%
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .95) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + theme_bw(base_size = 10) +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        #strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(#title = "Phylogengetic Diversity (PD)",
       title = "Full dataset",  
       color = "Significance", fill = "Significance", 
       x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato"))

# pruned
map_prun_pd <- pr_pd_sig_spat_df %>% 
  filter(samp_frac == .5) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + theme_bw(base_size = 10) +
  theme(legend.position = "none",
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        #legend.key = element_rect(fill = "gray80"),
        # legend.key.size = unit(.5, "cm"),
        #legend.text = element_text(size = 10),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        #strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(#title = "", 
    subtitle = "Pruned-only", 
    color = "Significance", fill = "Significance", 
    x = "", y = "") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  facet_wrap(~treat, ncol = 3)

# random imputation
map_rand_pd <- rand_imp_pd_sig_spat_df %>% 
  filter(samp_frac == .5) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + theme_bw(base_size = 10) +
  theme(legend.position = "none", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        #legend.key = element_rect(fill = "gray80"),
        #legend.key.size = unit(.5, "cm"),
        #legend.text = element_text(size = 10),
        strip.text = element_blank(),
        strip.background = element_blank(),
        #strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(#title = "Random imputation", 
    subtitle = "Random imputation", 
    color = "Significance", fill = "Significance", 
    x = "", y = "") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  facet_wrap(~treat, ncol = 3)

# imputed root
map_imp_tax_root_pd <- root_imp_pd_sig_spat_df %>% 
  filter(samp_frac == .5) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + theme_bw(base_size = 10) +
  theme(legend.position = "none", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        #legend.key = element_rect(fill = "gray80"),
        #legend.key.size = unit(.5, "cm"),
        #legend.text = element_text(size = 10),
        strip.text = element_blank(),
        strip.background = element_blank(),
        #strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(#title = "Taxonomic imputation", 
    subtitle = "Taxonomic imputation at node", 
    color = "Significance", fill = "Significance", 
    x = "", y = "") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  facet_wrap(~treat, ncol = 3)

# imputed random
map_imp_tax_rand_pd <- tax_rand_imp_pd_sig_spat_df %>% 
  filter(samp_frac == .5) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + theme_bw(base_size = 10) +
  theme(legend.position = "none",  
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        #legend.key = element_rect(fill = "gray80"),
        #legend.key.size = unit(.5, "cm"),
        #legend.text = element_text(size = 10),
        strip.text = element_blank(),
        strip.background = element_blank(),
        #strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(#title = "Taxonomic imputation", 
    subtitle = "Taxonomic imputation random within subtree", 
    color = "Significance", fill = "Significance", 
    x = "", y = "") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  facet_wrap(~treat, ncol = 3)

## RPE ####
# full data
map_emp_rpe <- emp_pd_sig_df %>% 
  filter(metric == "RPE") %>%
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .95) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + theme_bw(base_size = 10) +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        #strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(#title = "Phylogengetic Diversity (PD)",
    title = "Full dataset",  
    color = "Significance", fill = "Significance", 
    x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato"))

# random imputation
map_rand_rpe <- rand_imp_rpe_sig_spat_df %>% 
  filter(samp_frac == .5) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + theme_bw(base_size = 10) +
  theme(legend.position = "none", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        #legend.key = element_rect(fill = "gray80"),
        #legend.key.size = unit(.5, "cm"),
        #legend.text = element_text(size = 10),
        strip.text = element_blank(),
        strip.background = element_blank(),
        #strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(#title = "Random imputation", 
    subtitle = "Random imputation", 
    color = "Significance", fill = "Significance", 
    x = "", y = "") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  facet_wrap(~treat, ncol = 3)

# imputed root
map_imp_tax_root_rpe <- root_imp_rpe_sig_spat_df %>% 
  filter(samp_frac == .5) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + theme_bw(base_size = 10) +
  theme(legend.position = "none", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        #legend.key = element_rect(fill = "gray80"),
        #legend.key.size = unit(.5, "cm"),
        #legend.text = element_text(size = 10),
        strip.text = element_blank(),
        strip.background = element_blank(),
        #strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(#title = "Taxonomic imputation", 
    subtitle = "Taxonomic imputation at node", 
    color = "Significance", fill = "Significance", 
    x = "", y = "") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  facet_wrap(~treat, ncol = 3)

# imputed random
map_imp_tax_rand_rpe <- tax_rand_imp_rpe_sig_spat_df %>% 
  filter(samp_frac == .5) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = sign_cat, fill = sign_cat), shape = 20, size = .5) + 
  #geom_text(aes(label = sign_cat_emp), color = "black", size = 1) +
  coord_equal() + theme_bw(base_size = 10) +
  theme(legend.position = "none",  
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        #legend.key = element_rect(fill = "gray80"),
        #legend.key.size = unit(.5, "cm"),
        #legend.text = element_text(size = 10),
        strip.text = element_blank(),
        strip.background = element_blank(),
        #strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(#title = "Taxonomic imputation", 
    subtitle = "Taxonomic imputation random within subtree", 
    color = "Significance", fill = "Significance", 
    x = "", y = "") +
  scale_color_manual(values = c("wheat", "steelblue", "tomato")) +
  facet_wrap(~treat, ncol = 3)

# fig panels ####
library(patchwork)

# map panels
pl_lay <- "
###BBBB
AAACCCC
AAADDDD
###EEEE
"
free(map_emp_pd) + 
  (map_prun_pd) + 
  (map_rand_pd) + 
  (map_imp_tax_root_pd) + 
  (map_imp_tax_rand_pd) + 
  plot_layout(guides = "keep", design = pl_lay, axis_titles = "keep", 
              ) +
  plot_annotation(
    title = 'Phylogenetic Diversity (PD)'
    #caption = 'Disclaimer: None of these plots are insightful'
  )
ggsave(filename = paste0("figures/map_pd_comp-treats", frt_out), 
       width = 30, height = 20, units = "cm", dpi = 200)


pl_lay <- "
AAABBBB
AAACCCC
AAADDDD
"
(map_emp_rpe) + 
  (map_rand_rpe) + 
  (map_imp_tax_root_rpe) + 
  (map_imp_tax_rand_rpe) + 
  plot_layout(guides = "keep", design = pl_lay, axis_titles = "keep", 
  ) +
  plot_annotation(
    title = 'Relative Phylogenetic Endemism (RPE)'
    #caption = 'Disclaimer: None of these plots are insightful'
  )
ggsave(filename = paste0("figures/map_rpe_comp-treats", frt_out), 
       width = 30, height = 20, units = "cm", dpi = 200)

# grid proportional significance
pl_lay <- "
AAAA
BBBB
CCCC
DDDD
"

# PD
pd_res_type_pr + pd_res_type_rand_imp + pd_res_type_root_imp + pd_res_type_tax_rand_imp + 
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/grid_prop_pd_sign_pnl", frt_out), 
       width = 20, height = 25, units = "cm", dpi = 200)

# PE
pe_res_type_pr + pe_res_type_rand_imp + pe_res_type_root_imp + pe_res_type_tax_rand_imp + 
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/grid_prop_pe_sign_pnl", frt_out), 
       width = 20, height = 25, units = "cm", dpi = 200)

# RPD
rpd_res_type_rand_imp / rpd_res_type_root_imp / rpd_res_type_tax_rand_imp +
  plot_layout(guides = "collect", axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/grid_prop_rpd_sign_pnl", frt_out), 
       width = 20, height = 25, units = "cm", dpi = 200)

# RPE
rpe_res_type_rand_imp / rpe_res_type_root_imp / rpe_res_type_tax_rand_imp  +
  plot_layout(guides = "collect", axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/grid_prop_rpe_sign_pnl", frt_out), 
       width = 20, height = 25, units = "cm", dpi = 200)

