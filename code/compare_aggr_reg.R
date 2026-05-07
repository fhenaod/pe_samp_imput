library(tidyverse)

# load results using Nitta et al. 2022 dataset ####
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))

cols2sel <- c("site", 
  "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
  "pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean")

emp_df <- emp_canape_res %>% select(all_of(cols2sel)) %>% rename_if(is.numeric, ~paste0("emp_", .x)) ; emp_df %>% dim()

frt_out <- ".pdf" #  ".png"

# load pruned ####
rnd_canaper_fls <- list.files("output/canaper/", "canaper_rand_pr_[0-9]", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper/", "canaper_late_pr_[0-9]", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper/", "canaper_marg_pr_[0-9]", all.files = T, full.names = T)
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

prun_res_all <- rbind(rnd_canape_all_df, late_canape_all_df, marg_canape_all_df)

# load pruned  - random imputed ####
rnd_canaper_fls <- list.files("output/canaper/", "canaper_rand_pr_rand", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

late_canaper_fls <- list.files("output/canaper/", "canaper_late_pr_rand", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

marg_canaper_fls <- list.files("output/canaper/", "canaper_marg_pr_rand", all.files = T, full.names = T)
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

rand_imp_res_all <- rbind(rnd_canape_all_df, late_canape_all_df, marg_canape_all_df)

# load pruned - taxonomic imputed ####
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
        samp_frac = rand_frc[i], treat = "rand-pr", imputation = "taxonomy") %>% 
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
        samp_frac = late_frc[i], treat = "late-pr", imputation = "taxonomy") %>% 
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
        samp_frac = marg_frc[i], treat = "marg-pr", imputation = "taxonomy") %>% 
    bind_rows()
}
marg_canape_all_df <- data.table::rbindlist(marg_canape_df_ls) ; dim(marg_canape_all_df)

tax_imp_res_all <- rbind(rnd_canape_all_df, late_canape_all_df, marg_canape_all_df)

# all data  results raw #####
all_canape_res <- rbind(prun_res_all, rand_imp_res_all, tax_imp_res_all)

saveRDS(all_canape_res, "output/canaper_all_data_res.rds")
all_canape_res <- readRDS("output/canaper_all_data_res.rds")

# mean and sd per metric
canape_stats <- all_canape_res %>% 
  group_by(site, treat, samp_frac, imputation) %>%
  summarize(m_pd = mean(pd_obs), sd_pd = sd(pd_obs),
            m_pe = mean(pe_obs), sd_pe = sd(pe_obs),
            m_rpd = mean(rpd_obs), sd_rpd = sd(rpd_obs),
            m_rpe = mean(rpe_obs), sd_rpe = sd(rpe_obs),
            m_pd_m = mean(pd_rand_mean), sd_pd_m = sd(pd_rand_mean),
            m_pe_m = mean(pe_rand_mean), sd_pe_m = sd(pe_rand_mean),
            m_rpd_m = mean(rpd_rand_mean), sd_rpd_m = sd(rpd_rand_mean),
            m_rpe_m = mean(rpe_rand_mean), sd_rpe_m = sd(rpe_rand_mean)
            ) %>% 
  mutate(treat = case_match(treat, 
                            "late-pr" ~ "Lateral", 
                            "marg-pr" ~ "Marginal", 
                            "rand-pr" ~ "Random")) %>% 
  left_join(., emp_df, by = join_by(site))

head(canape_stats) %>% as.tibble()

# metrics to estimate and viz
estim2sel <- c("m_pd", "m_pe", "m_rpd", "m_rpe", 
               "m_pd_m", "m_pe_m", "m_rpd_m", "m_rpe_m",
               "m_pd_obs_z", "m_pe_obs_z", "m_rpd_obs_z", "m_rpe_obs_z") 
metir2sel <- c("pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
               "pd_rand_mean", "pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean",
               "pd_obs_z", "pe_obs_z", "rpd_obs_z", "rpe_obs_z") 

pars_ls_dfs <- list()
for(j in 1:length(cols2sel)){
  pd_obs_reg <- canape_stats %>% 
    nest(data = -imputation) %>% 
    mutate(fit = purrr::map(data, 
                            ~lm(eval(parse(text = paste0(".$", estim2sel[j]))) ~ eval(parse(text = paste0(".$emp_", metir2sel[j],"*.$treat")))), 
                            data = .),
           #purrr::map(data, ~lm(.$pd_obs ~ emp_df$pe_obs), data = .),
           glance = purrr::map(fit, broom::glance),
           tidied = purrr::map(fit, broom::tidy, conf.level = .95)
           #,augmented = purrr::map(fit, broom::augment)
    )
}

data.frame(pd_obs_reg %>% unnest(glance) %>% select(r.squared), 
           pd_obs_reg %>% unnest(tidied) %>% filter(term != "(Intercept)") %>% 
             select(estimate),
           metric = cols2sel[j]) 

# regression panels ####
## PD ####
canape_stats %>% 
  filter(samp_frac %in% c(0.05, 0.15, 0.3, 0.5)
         #, imputation != "random"
         ) %>%
  select(site, treat, samp_frac, imputation, contains("_pd")) %>%
  #gather( key, value, -site, -treat, -samp_frac, -imputation) %>% head()
  ggplot(aes(x = emp_pd_obs, y = m_pd, color = treat)) +
  geom_abline(intercept = 0, slope = 1, col = "coral", lty = 2, lwd = 1) +
  geom_pointrange(aes(ymin = m_pd - sd_pd, ymax =  m_pd + sd_pd), 
                  size = .1, position = position_dodge(0.05), alpha = 1/10) +
  geom_smooth(method = "lm", se  = T) +
  theme_classic(base_size = 12)  + 
  scale_colour_viridis_d(option = "viridis", begin = .95, end = .05 , alpha = .85) +
  theme(strip.background = element_blank(), legend.position = "top") +
  labs(subtitle = "Phylogenetic Diversity (PD)", 
       x = "PD observed full", y = "Mean PD", color = "Pruning\nmethod") +
  facet_wrap(~imputation*samp_frac, nrow = 3, scales = "free")
ggsave(filename = paste0("figures/mean_pd_reg_pnl", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

## PE ####
canape_stats %>% 
  filter(samp_frac %in% c(0.05, 0.15, 0.3, 0.5)
         #, imputation != "random"
  ) %>%
  select(site, treat, samp_frac, imputation, contains("_pe")) %>%
  #gather( key, value, -site, -treat, -samp_frac, -imputation) %>% head()
  ggplot(aes(x = emp_pe_obs, y = m_pe, color = treat)) +
  geom_abline(intercept = 0, slope = 1, col = "coral", lty = 2, lwd = 1) +
  geom_pointrange(aes(ymin = m_pe - sd_pe, ymax =  m_pe + sd_pe), 
                  size = .1, position = position_dodge(0.05), alpha = 1/10) +
  geom_smooth(method = "lm", se  = T) +
  theme_test(base_size = 12)  + 
  scale_colour_viridis_d(option = "viridis", begin = .95, end = .05 , alpha = .75) +
  theme(strip.background = element_blank(), legend.position = "top") +
  labs(subtitle = "Phylogenetic Endemism (PE)", 
       x = "PE observed full", y = "Mean PE", color = "Pruning\nmethod") +
  facet_wrap(~imputation*samp_frac, nrow = 3, scales = "free")
ggsave(filename = paste0("figures/mean_pe_reg_pnl", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

## RPD ####
canape_stats %>% 
  filter(samp_frac %in% c(0.05, 0.15, 0.3, 0.5)
         #, imputation != "random"
  ) %>%
  select(site, treat, samp_frac, imputation, contains("_rpd")) %>% 
  #gather( key, value, -site, -treat, -samp_frac, -imputation) %>% head()
  ggplot(aes(x = emp_rpd_obs, y = m_rpd, color = treat)) +
  geom_abline(intercept = 0, slope = 1, col = "coral", lty = 2, lwd = 1) +
  geom_pointrange(aes(ymin = m_rpd - sd_rpd, ymax =  m_rpd + sd_rpd), 
                  size = .1, position = position_dodge(0.05), alpha = 1/10) +
  geom_smooth(method = "lm", se  = T) +
  theme_test(base_size = 12)  + 
  scale_colour_viridis_d(option = "viridis", begin = .95, end = .05 , alpha = .75) +
  theme(strip.background = element_blank(), legend.position = "top") +
  labs(subtitle = "Relative Phylogenetic Diversity (RPD)", 
       x = "RPD observed full", y = "Mean RPD", color = "Pruning\nmethod") +
  facet_wrap(~imputation*samp_frac, nrow = 3, scales = "free")
ggsave(filename = paste0("figures/mean_rpd_reg_pnl", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

## RPE ####
canape_stats %>% 
  filter(samp_frac %in% c(0.05, 0.15, 0.3, 0.5)
         #, imputation != "random"
  ) %>%
  select(site, treat, samp_frac, imputation, contains("_rpe")) %>% 
#gather( key, value, -site, -treat, -samp_frac, -imputation) %>% head()
ggplot(aes(x = emp_rpe_obs, y = m_rpe, color = treat)) +
  geom_abline(intercept = 0, slope = 1, col = "coral", lty = 2, lwd = 1) +
  geom_pointrange(aes(ymin = m_rpe - sd_rpe, ymax =  m_rpe + sd_rpe), 
                  size = .1, position = position_dodge(0.05), alpha = 1/10) +
  geom_smooth(method = "lm", se  = T) +
  theme_test(base_size = 12)  + 
  scale_colour_viridis_d(option = "viridis", begin = .95, end = .05 , alpha = .75) +
  theme(strip.background = element_blank(), legend.position = "top") +
  labs(subtitle = "Relative Phylogenetic Endemism (RPE)", 
       x = "RPE observed full", y = "Mean RPE", color = "Pruning\nmethod") +
  facet_wrap(~imputation*samp_frac, nrow = 3, scales = "free")
ggsave(filename = paste0("figures/mean_rpe_reg_pnl", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 300)

# maps #####
# get japan's layout
world <- map_data("world")
japan <- world[world$region == "Japan", ]

spat_df <- read.csv("data/nitta-et-al/results/japan_ferns_biodiv.csv", header = T) %>%
  select(grids, lat, long) %>% mutate(grids = as.character(grids))

# PD difference
canape_stats %>% 
  filter(samp_frac %in% c(0.05, 0.15, 0.3, 0.5)
         #, imputation != "random"
  ) %>%
  select(site, treat, samp_frac, imputation, contains("_pd")) %>% 
  mutate(m_pd_diff = m_pd - emp_pd_obs, 
         treat = paste0(treat, " - ", samp_frac)) %>% ungroup() %>% 
  left_join(., spat_df, 
            by = join_by(site == grids)) %>%
  select(site, lat, long, m_pd_diff, treat, imputation) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = m_pd_diff, fill = m_pd_diff), shape = 20, size = 1/10) + 
  coord_equal() + theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.45, "cm"),
        legend.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "", subtitle = "Phylogenetic Diversity (PD)", 
       color = "Difference", fill = "Difference", 
       x = "Longitude", y = "Latitude") +
  scale_color_gradientn(colours = (RColorBrewer::brewer.pal(5, "Spectral"))) +
  scale_fill_gradientn( colours = (RColorBrewer::brewer.pal(5, "Spectral"))) +
  facet_wrap(~imputation * treat, ncol = 9)
ggsave(filename = paste0("figures/map_pd_diff_pnl", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 350)

# RPD difference
canape_stats %>% 
  filter(samp_frac %in% c(0.05, 0.15, 0.3, 0.5)
         #, imputation != "random"
  ) %>%
  select(site, treat, samp_frac, imputation, contains("_rpd")) %>% 
  mutate(m_rpd_diff = m_rpd - emp_rpd_obs, 
         treat = paste0(treat, " - ", samp_frac)) %>% ungroup() %>% 
  left_join(., spat_df, 
            by = join_by(site == grids)) %>% 
  select(site, lat, long, m_rpd_diff, treat, imputation) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "darkgray", fill = "lightgray") + 
  geom_point(aes(col = m_rpd_diff, fill = m_rpd_diff), shape = 20, size = 1/50) + 
  coord_equal() + theme_bw() +
  theme(legend.position = "top", 
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(fill = "gray80"),
        legend.key.size = unit(.45, "cm"),
        legend.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.placement  = "outside",
        #axis.title.x     = element_blank(), axis.title.y     = element_blank()
  ) + 
  labs(title = "", subtitle = "Phylogenetic Relative Diversity (RPD)", 
       color = "Difference", fill = "Difference", 
       x = "Longitude", y = "Latitude") +
  scale_color_gradientn(colours = (RColorBrewer::brewer.pal(5, "Spectral"))) +
  scale_fill_gradientn( colours = (RColorBrewer::brewer.pal(5, "Spectral"))) +
  facet_wrap(~imputation * treat, ncol = 9)
ggsave(filename = paste0("figures/map_rpd_diff_pnl", frt_out), 
       width = 35, height = 25, units = "cm", dpi = 350)



