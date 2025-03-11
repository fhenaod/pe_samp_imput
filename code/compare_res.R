library(tidyverse)

# load results using Nitta et al. 2022 dataset ####
emp_canape_res <- readRDS(paste0("output/empirical_canaper.rds"))[[1]]
emp_canape_sig <- readRDS(paste0("output/empirical_canaper.rds"))[[2]]

cols2sel <- c(#"site", 
              "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              "pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean")

emp_df <- emp_canape_res %>% select(all_of(cols2sel))

# load results from canaper 
# pruned ####
# random 
rnd_canaper_fls <- list.files("output/", "canaper_rand", all.files = T, full.names = T)
rand_frc <- readr::parse_number(rnd_canaper_fls)

pars_ls_dfs_res <- list()
for(i in 1:length(rand_frc)){
  rnd_canape_ls <- readRDS(rnd_canaper_fls[i])
  rnd_canape_df <- lapply(rnd_canape_ls, "[[", 1) %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls)))) %>% 
    bind_rows()
  
  print(paste0("Iteration: ", i)); rnd_canape_df %>% group_by(sim_n) %>% tally()
  
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
res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
res_pars_df %>% saveRDS("output/pars_rand_df.rds") ; gc()

# lateral
late_canaper_fls <- list.files("output/", "canaper_late", all.files = T, full.names = T)
late_frc <- readr::parse_number(late_canaper_fls)

pars_ls_dfs_res <- list()
for(i in 1:length(late_canaper_fls)){
  late_canape_ls <- readRDS(late_canaper_fls[i])
  late_canape_df <- lapply(late_canape_ls, "[[", 1) %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls)))) %>% 
    bind_rows()
  
  if_else(is_empty(late_canape_df), "WARNING!!! Empty data frame", "")
  
  print(paste0("Iteration: ", i)); 
  late_canape_df %>% group_by(sim_n) %>% tally()
  emp_df %>% dim()
  
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
res_pars_df <- pars_ls_dfs_res %>% bind_rows() 
res_pars_df %>% saveRDS("output/pars_late_df.rds") ; gc()

# marginal
canaper_fls <- list.files("output/", "canaper_", all.files = T, full.names = T)
canaper_fls

## figs #####
# random 
res_pars_df %>% select(-estimate) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  ggplot(aes(x = samp_frac, y = value)) +
  geom_point() + geom_smooth() +
  labs(subtitle = "Random pruning", x = "Percentage of tips pruned", y = "R-squared") +
  theme_classic() +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 4, scales = "free")

res_pars_df %>% select(-r.squared) %>%
  gather(key, value, -sim_n, -samp_frac, -metric, -pruning) %>% 
  ggplot(aes(x = samp_frac, y = value)) +
  geom_point() + geom_smooth() +
  labs(subtitle = "Random pruning", x = "Percentage of tips pruned", y = "Estimate") +
  theme_classic() +
  facet_wrap(~factor(metric, c(cols2sel)), ncol = 4, scales = "free")

# imputed ####


