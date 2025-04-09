library(tidyverse)

# load results using Nitta et al. 2022 dataset ####
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))

cols2sel <- c("site", 
              "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              "pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean")

emp_df <- emp_canape_res %>% select(all_of(cols2sel))

# load results from canaper, run regressions, and compare r-square and slope
# pruned ####
# random 
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_[0-9]", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

run_reg <- T
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(rand_frc)){
    rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
    rnd_canape_df <- rnd_canape_ls %>% 
      lapply(., "[", cols2sel) %>%
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls)))) %>% 
      bind_rows() %>% 
      left_join(., emp_df %>% rename_if(is.numeric, ~paste0("emp_", .x)), by = join_by(site))
    
    print(paste0("Iteration: ", i)) #; rnd_canape_df %>% group_by(sim_n) %>% tally()
    
    pars_ls_dfs <- list()
    for(j in 2:length(cols2sel)){
      pd_obs_reg <- rnd_canape_df %>% 
        nest(data = -sim_n) %>% 
        mutate(fit = purrr::map(data, 
                                ~lm(eval(parse(text = paste0(".$", cols2sel[j]))) ~ eval(parse(text = paste0(".$emp_", cols2sel[j])))), 
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
late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_[0-9]", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(late_canaper_fls)){
    late_canape_ls <- readRDS(late_canaper_fls[i])
    late_canape_df <- late_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls)))) %>% 
      bind_rows() %>% 
      left_join(., emp_df %>% rename_if(is.numeric, ~paste0("emp_", .x)), by = join_by(site))
    
    if_else(is_empty(late_canape_df), "WARNING!!! Empty data frame", "")
    
    print(paste0("Iteration: ", i))
    
    pars_ls_dfs <- list()
    for(j in 2:length(cols2sel)){
      pd_obs_reg <- late_canape_df %>% 
        nest(data = -sim_n) %>% 
        mutate(fit = purrr::map(data, 
                                ~lm(eval(parse(text = paste0(".$", cols2sel[j]))) ~ eval(parse(text = paste0(".$emp_", cols2sel[j])))), 
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
marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_[0-9]", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(marg_canaper_fls)){
    marg_canape_ls <- readRDS(marg_canaper_fls[i])
    marg_canape_df <- marg_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canape_ls)))) %>% 
      bind_rows() %>% 
      left_join(., emp_df %>% rename_if(is.numeric, ~paste0("emp_", .x)), by = join_by(site))
    
    if_else(is_empty(marg_canape_df), "WARNING!!! Empty data frame", "")
    
    print(paste0("Iteration: ", i))
    
    pars_ls_dfs <- list()
    for(j in 2:length(cols2sel)){
      pd_obs_reg <- marg_canape_df %>% 
        nest(data = -sim_n) %>% 
        mutate(fit = purrr::map(data, 
                                ~lm(eval(parse(text = paste0(".$", cols2sel[j]))) ~ eval(parse(text = paste0(".$emp_", cols2sel[j])))), 
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

# metric comparison by observed or randomized
all_prun_res_df <- rbind(rand_res_pars_df, late_res_pars_df, marg_res_pars_df)

# imputed ####
## random imputation ####
rnd_canaper_fls <- list.files("output/canaper/", "canaper_rand_pr_rand", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(rand_frc)){
    rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
    rnd_canape_df <- rnd_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls)))) %>% 
      bind_rows()
    
    print(paste0("Iteration: ", i)) #; rnd_canape_df %>% group_by(sim_n) %>% tally()
    
    pars_ls_dfs <- list()
    for(j in 2:length(cols2sel)){ # start in 2 to avoid site column
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
                                     imputation = "random",
                                     samp_frac = rand_frc[i]) 
    }
    
    pars_ls_dfs_res[[i]] <- pars_ls_dfs %>% bind_rows() 
  }
  rand_res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
  rand_res_pars_df %>% saveRDS("output/pars_rand_pr_rand_imp.rds") ; gc()
} else {
  rand_res_pars_df <- readRDS("output/pars_rand_pr_rand_imp.rds")
}

# lateral
late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_rand", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(late_canaper_fls)){
    late_canape_ls <- readRDS(late_canaper_fls[i])
    late_canape_df <- late_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls)))) %>% 
      bind_rows()
    
    if_else(is_empty(late_canape_df), "WARNING!!! Empty data frame", "")
    
    print(paste0("Iteration: ", i))
    
    pars_ls_dfs <- list()
    for(j in 2:length(cols2sel)){
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
                                     imputation = "random",
                                     samp_frac = late_frc[i]) 
    }
    pars_ls_dfs_res[[i]] <- pars_ls_dfs %>% bind_rows() 
  }
  late_res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
  late_res_pars_df %>% saveRDS("output/pars_late_pr_rand_imp.rds") ; gc()
} else {
  late_res_pars_df <-  readRDS("output/pars_late_pr_rand_imp.rds")
}

# marginal
marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_rand", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(marg_canaper_fls)){
    marg_canape_ls <- readRDS(marg_canaper_fls[i])
    marg_canape_df <- marg_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canape_ls)))) %>% 
      bind_rows()
    
    if_else(is_empty(marg_canape_df), "WARNING!!! Empty data frame", "")
    
    print(paste0("Iteration: ", i))
    
    pars_ls_dfs <- list()
    for(j in 2:length(cols2sel)){
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
                                     imputation = "random",
                                     samp_frac = marg_frc[i]) 
    }
    pars_ls_dfs_res[[i]] <- pars_ls_dfs %>% bind_rows() 
  }
  marg_res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
  marg_res_pars_df %>% saveRDS("output/pars_marg_pr_rand_imp.rds") ; gc()
} else {
  marg_res_pars_df <- readRDS("output/pars_marg_pr_rand_imp.rds")
}

# metric comparison by observed or randomized
all_prun_rand_imp_res_df <- rbind(rand_res_pars_df, late_res_pars_df, marg_res_pars_df)

## imputed taxonomically ####
# random 
rnd_canaper_fls <- list.files("output/canaper", "canaper_rand_pr_tax", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(rand_frc)){
    rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
    rnd_canape_df <- rnd_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls)))) %>% 
      bind_rows()
    
    cat("Iteration: ", i) #; rnd_canape_df %>% group_by(sim_n) %>% tally()
    
    pars_ls_dfs <- list()
    for(j in 2:length(cols2sel)){ # start in 2 to avoid site column
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
                                     imputation = "taxonomy",
                                     samp_frac = rand_frc[i]) 
    }
    
    pars_ls_dfs_res[[i]] <- pars_ls_dfs %>% bind_rows() 
  }
  rand_res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
  rand_res_pars_df %>% saveRDS("output/pars_rand_pr_tax_imp.rds") ; gc()
} else {
  rand_res_pars_df <- readRDS("output/pars_rand_pr_tax_imp.rds")
}

# lateral
late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_tax", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(late_canaper_fls)){
    late_canape_ls <- readRDS(late_canaper_fls[i])
    late_canape_df <- late_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls)))) %>% 
      bind_rows()
    
    if_else(is_empty(late_canape_df), "WARNING!!! Empty data frame", "")
    
    print(paste0("Iteration: ", i))
    
    pars_ls_dfs <- list()
    for(j in 2:length(cols2sel)){
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
                                     imputation = "",
                                     samp_frac = late_frc[i]) 
    }
    pars_ls_dfs_res[[i]] <- pars_ls_dfs %>% bind_rows() 
  }
  late_res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
  late_res_pars_df %>% saveRDS("output/pars_late_pr_tax_imp.rds") ; gc()
} else {
  late_res_pars_df <-  readRDS("output/pars_late_pr_tax_imp.rds")
}

# marginal
marg_canaper_fls <- list.files("output/canaper", "canaper_marg_pr_tax", all.files = T, full.names = T)
marg_frc <- readr::parse_number(marg_canaper_fls)

run_reg <- F
if(run_reg){
  pars_ls_dfs_res <- list()
  for(i in 1:length(marg_canaper_fls)){
    marg_canape_ls <- readRDS(marg_canaper_fls[i])
    marg_canape_df <- marg_canape_ls %>% 
      lapply(., "[", cols2sel) %>% 
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(marg_canape_ls)))) %>% 
      bind_rows()
    
    if_else(is_empty(marg_canape_df), "WARNING!!! Empty data frame", "")
    
    print(paste0("Iteration: ", i))
    
    pars_ls_dfs <- list()
    for(j in 2:length(cols2sel)){
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
                                     imputation = "taxonomy",
                                     samp_frac = marg_frc[i]) 
    }
    pars_ls_dfs_res[[i]] <- pars_ls_dfs %>% bind_rows() 
  }
  marg_res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
  marg_res_pars_df %>% saveRDS("output/pars_marg_pr_tax_imp.rds") ; gc()
} else {
  marg_res_pars_df <- readRDS("output/pars_marg_pr_tax_imp.rds")
}

# metric comparison by observed or randomized
all_prun_tax_imp_res_df <- rbind(rand_res_pars_df, late_res_pars_df, marg_res_pars_df)

# figs #####
#  all one plot by estimate
# pruned 
all_prun_res_df %>% select(-estimate) %>% 
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  mutate(pruning = str_to_title(pruning)) %>%
  ggplot(aes(x = (samp_frac), y = (value), col = pruning)) +
  #geom_hline(yintercept = log(.5), lty = 2, col = "darkgray") +
  geom_point(cex = .65, position = position_jitterdodge(jitter.width = .05, dodge.width = .05)) +
  geom_smooth(method = "lm") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) +
  labs(subtitle = "Pruned", x = "Proportion of tips pruned", y = "R-squared", color = "Pruning\nmethod") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "top") +
  #scale_x_continuous(breaks = marg_frc, labels = as.character(marg_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 2, scales = "free")
ggsave(filename = paste0("figures/r2_pruned.png"), 
       width = 20, height = 25, units = "cm", dpi = 200)

all_prun_res_df %>% select(-r.squared) %>% 
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  mutate(pruning = str_to_title(pruning)) %>%
  ggplot(aes(x = (samp_frac), y = value, col = pruning)) +
  geom_hline(yintercept = 1, lty = 2, col = "darkgray") +
  geom_point(cex = .65, position = position_jitterdodge(jitter.width = .05, dodge.width = .05)) +
  geom_smooth(method = "lm") +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) +
  #scale_colour_viridis_d(option = "plasma", begin = .75, end = .05 , alpha = .75) +
  labs(subtitle = "Pruned", x = "Proportion of tips pruned", y = "Slope", color = "Pruning\nmethod") +
  theme_classic(base_size = 12)  + 
  theme(strip.background = element_blank(), legend.position = "top") +
  #scale_x_continuous(breaks = marg_frc, labels = as.character(marg_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 2, scales = "free")
ggsave(filename = paste0("figures/slope_pruned.png"), 
       width = 20, height = 25, units = "cm", dpi = 200)

# imputed random
all_prun_rand_imp_res_df %>% select(-estimate) %>% 
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning, -imputation) %>% 
  mutate(pruning = str_to_title(pruning)) %>%
  ggplot(aes(x = (samp_frac), y = (value), col = pruning)) +
  #geom_hline(yintercept = log(.5), lty = 2, col = "darkgray") +
  geom_point(cex = .65, position = position_jitterdodge(jitter.width = .05, dodge.width = .05)) +
  geom_smooth(method = "lm") +
  scale_color_manual(values = MetBrewer::met.brewer("Renoir", 3)) + 
  labs(subtitle = "Random imputation", x = "Proportion of tips pruned", y = "R-squared", color = "Pruning\nmethod") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "top") +
  #scale_x_continuous(breaks = marg_frc, labels = as.character(marg_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 2, scales = "free")
ggsave(filename = paste0("figures/r2_pruned_rand_imp.png"), 
       width = 20, height = 25, units = "cm", dpi = 200)

all_prun_rand_imp_res_df %>% select(-r.squared) %>% 
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning, -imputation) %>% 
  mutate(pruning = str_to_title(pruning)) %>%
  ggplot(aes(x = (samp_frac), y = value, col = pruning)) +
  geom_hline(yintercept = 1, lty = 2, col = "darkgray") +
  geom_point(cex = .65, position = position_jitterdodge(jitter.width = .05, dodge.width = .05)) +
  geom_smooth(method = "lm") +
  scale_color_manual(values = MetBrewer::met.brewer("Renoir", 3)) + 
  labs(subtitle = "Random imputation", x = "Proportion of tips pruned", y = "Slope", color = "Pruning\nmethod") +
  theme_classic(base_size = 12)  + 
  theme(strip.background = element_blank(), legend.position = "top") +
  #scale_x_continuous(breaks = marg_frc, labels = as.character(marg_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 2, scales = "free")
ggsave(filename = paste0("figures/slope_pruned_rand_imp.png"), 
       width = 20, height = 25, units = "cm", dpi = 200)

# imputed taxonomic
all_prun_tax_imp_res_df %>% select(-estimate) %>% 
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning, -imputation) %>% 
  mutate(pruning = str_to_title(pruning)) %>%
  ggplot(aes(x = (samp_frac), y = (value), col = pruning)) +
  #geom_hline(yintercept = log(.5), lty = 2, col = "darkgray") +
  geom_point(cex = .65, position = position_jitterdodge(jitter.width = .05, dodge.width = .05)) +
  geom_smooth(method = "lm") +
  scale_colour_viridis_d(option = "viridis", begin = .9, end = .05, alpha = .95) +
  labs(subtitle = "Taxonomic imputation", x = "Proportion of tips pruned", y = "R-squared", color = "Pruning\nmethod") +
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), legend.position = "top") +
  #scale_x_continuous(breaks = marg_frc, labels = as.character(marg_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 2, scales = "free")
ggsave(filename = paste0("figures/r2_pruned_tax_imp.png"), 
       width = 20, height = 25, units = "cm", dpi = 200)

all_prun_tax_imp_res_df %>% select(-r.squared) %>% 
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning, -imputation) %>% 
  mutate(pruning = str_to_title(pruning)) %>%
  ggplot(aes(x = (samp_frac), y = value, col = pruning)) +
  geom_hline(yintercept = 1, lty = 2, col = "darkgray") +
  geom_point(cex = .65, position = position_jitterdodge(jitter.width = .05, dodge.width = .05)) +
  geom_smooth(method = "lm") +
  scale_colour_viridis_d(option = "viridis", begin = .9, end = .05, alpha = .95) +
  #scale_colour_viridis_d(option = "plasma", begin = .75, end = .05 , alpha = .75) +
  labs(subtitle = "Taxonomic imputation", x = "Proportion of tips pruned", y = "Slope", color = "Pruning\nmethod") +
  theme_classic(base_size = 12)  + 
  theme(strip.background = element_blank(), legend.position = "top") +
  #scale_x_continuous(breaks = marg_frc, labels = as.character(marg_frc)) +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 2, scales = "free")
ggsave(filename = paste0("figures/slope_pruned_tax_imp.png"), 
       width = 20, height = 25, units = "cm", dpi = 200)

# viz regs  ####
library(ggpubr)
library(patchwork)
late_canaper_fls <- list.files("output/canaper", "canaper_late_pr_tax", all.files = T, full.names = T)
canape_df2 <- readRDS(late_canaper_fls[4])[[15]]
canape_df3 <- readRDS(late_canaper_fls[4])[[32]]
canape_df4 <- readRDS(late_canaper_fls[4])[[68]]
canape_df5 <- readRDS(late_canaper_fls[4])[[97]]

p1 <- emp_df %>% 
  ggscatter(x = "pd_obs", y = "pd_obs", color = "black", fill = "black", size = 1.5, alpha = .65) +
  geom_smooth(lwd = .5) +
  stat_regline_equation(label.y = .54) +
  stat_cor(label.y = .52) +
  labs(subtitle = "Full dataset", x = "PD observed full", y = "PD observed full") +
  theme_classic()

p2 <- data.frame(pd_obs_prun = canape_df2$pd_obs, 
                 pd_obs_emp = emp_df$pd_obs) %>% 
  ggscatter(x = "pd_obs_emp", y = "pd_obs_prun", color = "black", fill = "black", size = .75, alpha = .85) +
  geom_abline(intercept = 0, slope = 1, col = "coral", lty = 2, lwd = 1) +
  geom_smooth(lwd = .5) +
  stat_regline_equation(label.y = 1.45, cex = 3.5) +
  stat_cor(label.y = 1.35, cex = 3.5) + 
  labs(subtitle = "Simulation 15", x = "", y = "PD observed pruned") +
  theme_classic()

p3 <- data.frame(pd_obs_prun = canape_df3$pd_obs, 
                 pd_obs_emp = emp_df$pd_obs) %>% 
  ggscatter(x = "pd_obs_emp", y = "pd_obs_prun", color = "black", fill = "black", size = .75, alpha = .85) +
  geom_abline(intercept = 0, slope = 1, col = "coral", lty = 2, lwd = 1) +
  geom_smooth(lwd = .5) +
  stat_regline_equation(label.y = 1.45, cex = 3.5) +
  stat_cor(label.y = 1.35, cex = 3.5) + 
  labs(subtitle = "Simulation 32", x = "PD observed full", y = "PD observed pruned") +
  theme_classic()

p4 <- data.frame(pd_obs_prun = canape_df4$pd_obs, 
                 pd_obs_emp = emp_df$pd_obs) %>% 
  ggscatter(x = "pd_obs_emp", y = "pd_obs_prun", color = "black", fill = "black", size = .75, alpha = .85) +
  geom_abline(intercept = 0, slope = 1, col = "coral", lty = 2, lwd = 1) +
  geom_smooth(lwd = .5) +
  stat_regline_equation(label.y = 1.45, cex = 3.5) +
  stat_cor(label.y = 1.35, cex = 3.5) + 
  labs(subtitle = "Simulation 68", x = "PD observed full", y = "PD observed pruned") +
  theme_classic()

p5 <- data.frame(pd_obs_prun = canape_df5$pd_obs, 
                 pd_obs_emp = emp_df$pd_obs) %>% 
  ggscatter(x = "pd_obs_emp", y = "pd_obs_prun", color = "black", fill = "black", size = .75, alpha = .85) +
  geom_abline(intercept = 0, slope = 1, col = "coral", lty = 2, lwd = 1) +
  geom_smooth(lwd = .5) +
  stat_regline_equation(label.y = 1.45, cex = 3.5) +
  stat_cor(label.y = 1.35, cex = 3.5) + 
  labs(subtitle = "Simulation 97", x = "PD observed full", y = "PD observed pruned") +
  theme_classic()

p1 + p2 + p3 + p4 + p5 + 
plot_layout(design = "AAAABBCC
                      AAAADDEE")
#ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
ggsave(filename = paste0("figures/regs_comp.png"), 
       width = 25, height = 25, units = "cm", dpi = 300)
 
