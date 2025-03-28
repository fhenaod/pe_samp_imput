library(tidyverse)

# load results using Nitta et al. 2022 dataset ####
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))[[1]]
emp_canape_sig <- readRDS(paste0("output/canaper/canaper_empirical.rds"))[[2]]

cols2sel <- c(#"site", 
              "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              "pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean")

emp_df <- emp_canape_res %>% select(all_of(cols2sel))

# load results from canaper 
# pruned ####
# random 
rnd_canaper_fls <- list.files("output/canaper/", "canaper_rand", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(rand_frc)){
    rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
    rnd_canape_df <- lapply(rnd_canape_ls, "[[", 1) %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls)))) %>% 
      bind_rows()
    
    print(paste0("Iteration: ", i)) #; rnd_canape_df %>% group_by(sim_n) %>% tally()
    
    pars_ls_dfs <- list()
    for(j in 1:length(cols2sel)){
      pd_obs_reg <- rnd_canape_df %>% 
        nest(data = -sim_n) %>% 
        mutate(fit = purrr::map(data, 
                                ~lm(eval(parse(text = paste0(".$", cols2sel[j]))) ~ eval(parse(text = paste0("emp_df$", cols2sel[j])))), 
                                data = .),
               #purrr::map(data, ~lm(.$pd_obs ~ emp_df$pe_obs), data = .),
               glance = purrr::map(fit, broom::glance),
               tidied = purrr::map(fit, broom::tidy, conf.level = .95)
               #,augmented = purrr::map(fit, broom::augment)
        )
      
      pars_ls_dfs[[j]] <- data.frame(pd_obs_reg %>% unnest(glance) %>% select(sim_n, r.squared), 
                                     pd_obs_reg %>% unnest(tidied) %>% filter(term != "(Intercept)") %>% 
                                       select(estimate),
                                     metric = cols2sel[j],
                                     pruning = "random",
                                     samp_frac = rand_frc[i]) 
    }
    
    pars_ls_dfs_res[[i]] <- pars_ls_dfs %>% bind_rows() 
  }
  rand_res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
  rand_res_pars_df %>% saveRDS("output/pars_rand_prun.rds") ; gc()
} else {
  rand_res_pars_df <- readRDS("output/pars_rand_prun.rds")
}

# lateral
late_canaper_fls <- list.files("output/canaper/", "canaper_late", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(late_canaper_fls)){
    late_canape_ls <- readRDS(late_canaper_fls[i])
    late_canape_df <- lapply(late_canape_ls, "[[", 1) %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls)))) %>% 
      bind_rows()
    
    if_else(is_empty(late_canape_df), "WARNING!!! Empty data frame", "")
    
    print(paste0("Iteration: ", i))
    
    pars_ls_dfs <- list()
    for(j in 1:length(cols2sel)){
      pd_obs_reg <- late_canape_df %>% 
        nest(data = -sim_n) %>% 
        mutate(fit = purrr::map(data, 
                                ~lm(eval(parse(text = paste0(".$", cols2sel[j]))) ~ eval(parse(text = paste0("emp_df$", cols2sel[j])))), 
                                data = .),
               #purrr::map(data, ~lm(.$pd_obs ~ emp_df$pe_obs), data = .),
               glance = purrr::map(fit, broom::glance),
               tidied = purrr::map(fit, broom::tidy, conf.level = .95)
               #,augmented = purrr::map(fit, broom::augment)
        )
      
      pars_ls_dfs[[j]] <- data.frame(pd_obs_reg %>% unnest(glance) %>% select(sim_n, r.squared), 
                                     pd_obs_reg %>% unnest(tidied) %>% filter(term != "(Intercept)") %>% 
                                       select(estimate),
                                     metric = cols2sel[j],
                                     pruning = "lateral",
                                     samp_frac = late_frc[i]) 
    }
    pars_ls_dfs_res[[i]] <- pars_ls_dfs %>% bind_rows() 
  }
  late_res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
  late_res_pars_df %>% saveRDS("output/pars_late_prun.rds") ; gc()
} else {
  late_res_pars_df <- readRDS("output/pars_late_prun.rds")
}

# marginal
marg_canaper_fls <- list.files("output/canaper/", "canaper_marg", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(marg_canaper_fls)){
    marg_canape_ls <- readRDS(marg_canaper_fls[i])
    marg_canape_df <- lapply(marg_canape_ls, "[[", 1) %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canape_ls)))) %>% 
      bind_rows()
    
    if_else(is_empty(marg_canape_df), "WARNING!!! Empty data frame", "")
    
    print(paste0("Iteration: ", i))
    
    pars_ls_dfs <- list()
    for(j in 1:length(cols2sel)){
      pd_obs_reg <- marg_canape_df %>% 
        nest(data = -sim_n) %>% 
        mutate(fit = purrr::map(data, 
                                ~lm(eval(parse(text = paste0(".$", cols2sel[j]))) ~ eval(parse(text = paste0("emp_df$", cols2sel[j])))), 
                                data = .),
               #purrr::map(data, ~lm(.$pd_obs ~ emp_df$pe_obs), data = .),
               glance = purrr::map(fit, broom::glance),
               tidied = purrr::map(fit, broom::tidy, conf.level = .95)
               #,augmented = purrr::map(fit, broom::augment)
        )
      
      pars_ls_dfs[[j]] <- data.frame(pd_obs_reg %>% unnest(glance) %>% select(sim_n, r.squared), 
                                     pd_obs_reg %>% unnest(tidied) %>% filter(term != "(Intercept)") %>% 
                                       select(estimate),
                                     metric = cols2sel[j],
                                     pruning = "marginal",
                                     samp_frac = marg_frc[i]) 
    }
    pars_ls_dfs_res[[i]] <- pars_ls_dfs %>% bind_rows() 
  }
  marg_res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
  marg_res_pars_df %>% saveRDS("output/pars_marg_prun.rds") ; gc()
} else {
  marg_res_pars_df <- readRDS("output/pars_marg_prun.rds")
}

# figs #####
# by pruning method
# random 
rand_res_pars_df %>% select(-estimate) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  ggplot(aes(x = as.factor(samp_frac), y = value)) +
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + geom_smooth() +
  labs(subtitle = "Random pruning", x = "Proportion of tips pruned", y = "R-squared") +
  theme_classic(base_size = 15) + theme(strip.background = element_blank()) +
  #scale_x_continuous(breaks = rand_frc, labels = as.character(rand_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 4, scales = "fixed")
ggsave(filename = paste0("figures/r2_rand_pr.svg"), 
       width = 25, height = 20, units = "cm", dpi = 200)

rand_res_pars_df %>% select(-r.squared) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  ggplot(aes(x = as.factor(samp_frac), y = value)) +
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + #geom_smooth() +
  labs(subtitle = "Random pruning", x = "Proportion of tips pruned", y = "Slope") +
  theme_classic(base_size = 15) + theme(strip.background = element_blank()) +
  #scale_x_continuous(breaks = rand_frc, labels = as.character(rand_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 4, scales = "fixed")
ggsave(filename = paste0("figures/slope_rand_pr.svg"), 
       width = 25, height = 20, units = "cm", dpi = 200)

# lateral
late_res_pars_df %>% select(-estimate) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  ggplot(aes(x = as.factor(samp_frac), y = value)) +
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + geom_smooth() +
  labs(subtitle = "Lateral pruning", x = "Proportion of tips pruned", y = "R-squared") +
  theme_classic(base_size = 15) + theme(strip.background = element_blank()) +
  #scale_x_continuous(breaks = late_frc, labels = as.character(late_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 4, scales = "fixed")
ggsave(filename = paste0("figures/r2_late_pr.svg"), 
       width = 25, height = 20, units = "cm", dpi = 200)

late_res_pars_df %>% select(-r.squared) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  ggplot(aes(x = as.factor(samp_frac), y = value)) +
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + geom_smooth() +
  labs(subtitle = "Lateral pruning", x = "Proportion of tips pruned", y = "Slope") +
  theme_classic(base_size = 15) + theme(strip.background = element_blank()) +
  #scale_x_continuous(breaks = late_frc, labels = as.character(late_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 4, scales = "fixed")
ggsave(filename = paste0("figures/slope_late_pr.svg"), 
       width = 25, height = 20, units = "cm", dpi = 200)

# marginal
marg_res_pars_df %>% select(-estimate) %>% 
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  ggplot(aes(x = as.factor(samp_frac), y = value)) +
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + geom_smooth() +
  labs(subtitle = "Marginal pruning", x = "Proportion of tips pruned", y = "R-squared") +
  theme_classic(base_size = 15)  + theme(strip.background = element_blank()) +
  #scale_x_continuous(breaks = marg_frc, labels = as.character(marg_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 4, scales = "fixed")
ggsave(filename = paste0("figures/r2_marg_pr.svg"), 
       width = 25, height = 20, units = "cm", dpi = 200)

marg_res_pars_df %>% select(-r.squared) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  ggplot(aes(x = as.factor(samp_frac), y = value)) +
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + geom_smooth() +
  labs(subtitle = "Marginal pruning", x = "Proportion of tips pruned", y = "Slope") +
  theme_classic(base_size = 15)  + theme(strip.background = element_blank()) +
  #scale_x_continuous(breaks = marg_frc, labels = as.character(marg_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 4, scales = "fixed")
ggsave(filename = paste0("figures/slope_marg_pr.svg"), 
       width = 25, height = 20, units = "cm", dpi = 200)

# metric comparison by observed or randomized
all_prun_res_df <- rbind(rand_res_pars_df, late_res_pars_df, marg_res_pars_df)

# observed
all_prun_res_df %>% filter(grepl("obs", metric)) %>% 
  select(-estimate) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>%
  ggplot(aes(x = as.factor(samp_frac), y = value)) + 
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + geom_smooth() + 
  labs(subtitle = "Observed", x = "Proportion of tips pruned", y = "R-squared") +
  theme_classic(base_size = 15)  + theme(strip.background = element_blank()) +
  facet_wrap(~metric*pruning, ncol = 3, scales = "fixed" #, labeller = labeller("PD")
             )
ggsave(filename = paste0("figures/r2_obs.svg"), 
       width = 20, height = 20, units = "cm", dpi = 200)

all_prun_res_df %>% filter(grepl("obs", metric)) %>% 
  select(-r.squared) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>%
  ggplot(aes(x = as.factor(samp_frac), y = value)) + 
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + geom_smooth() + 
  labs(subtitle = "Observed", x = "Proportion of tips pruned", y = "Slope") +
  theme_classic(base_size = 15) + theme(strip.background = element_blank()) +
  facet_wrap(~metric*pruning, ncol = 3, scales = "fixed" #, labeller = labeller("PD")
  )
ggsave(filename = paste0("figures/slope_obs.svg"), 
       width = 20, height = 20, units = "cm", dpi = 200)

# randomized observed
all_prun_res_df %>% filter(grepl("rand", metric)) %>% 
  select(-estimate) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>%
  ggplot(aes(x = as.factor(samp_frac), y = value)) + 
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + geom_smooth() + 
  labs(subtitle = "Randomized", x = "Proportion of tips pruned", y = "R-squared") +
  theme_classic(base_size = 15) + theme(strip.background = element_blank()) +
  facet_wrap(~metric*pruning, ncol = 3, scales = "free" #, labeller = labeller("PD")
  )
ggsave(filename = paste0("figures/r2_rand.svg"), 
       width = 20, height = 20, units = "cm", dpi = 200)

all_prun_res_df %>% filter(grepl("rand", metric)) %>% 
  select(-r.squared) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>%
  ggplot(aes(x = as.factor(samp_frac), y = value)) + 
  geom_boxplot(outlier.size = .85, varwidth = T) +
  geom_point(cex = .65) + geom_smooth() + 
  labs(subtitle = "Randomized", x = "Proportion of tips pruned", y = "Slope") +
  theme_classic(base_size = 15) + theme(strip.background = element_blank()) +
  facet_wrap(~metric*pruning, ncol = 3, scales = "free" #, labeller = labeller("PD")
  )
ggsave(filename = paste0("figures/slope_rand.svg"), 
       width = 20, height = 20, units = "cm", dpi = 200)

# imputed tips ####


# viz regs idea ####
library(ggpubr)
p1 <- emp_df %>% 
  ggscatter(x = "pd_obs", y = "pd_obs") +
  geom_point(cex = .5, col = "black") +
  geom_smooth() +
  stat_regline_equation(label.y = .54) +
  stat_cor(label.y = .52) +
  labs(x = "PD observed full", y = "PD observed full") +
  theme_classic()

p2 <- emp_df %>% 
  rowwise() %>% 
  mutate(pd_obs_add = c(pd_obs+runif(1, min = 0, max = 1))) %>% 
  select(pd_obs, pd_obs_add) %>% #head()
  ggscatter(x = "pd_obs", y = "pd_obs_add") +
  geom_point(cex = .5, col = "black") +
  geom_smooth() +
  stat_regline_equation(label.y = 1.45,) +
  stat_cor(label.y = 1.4) + 
  labs(x = "PD observed full", y = "PD observed pruned") +
  theme_classic()

ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
ggsave(filename = paste0("figures/regs_comp2.png"), 
       width = 25, height = 15, units = "cm", dpi = 200)

