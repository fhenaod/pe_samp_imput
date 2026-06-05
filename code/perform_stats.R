library(tidyverse)
library(patchwork)

# load functions
source("code/pe_samp_funs.R")

# endemicity categories and significance ####
# columns to select for analyses
cols2sel <- c("site", "endem_type",
              "pd_signif", "pe_signif"
              #, "rpd_signif", "rpe_signif"
              )

## load results using Nitta et al. 2022 dataset ####
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))
emp_df <- emp_canape_res %>% select(cols2sel) ; emp_df %>% head()

spat_df <- read.csv("data/nitta-et-al/results/japan_ferns_biodiv.csv", header = T) %>% 
  select(grids, lat, long) %>% mutate(site = as.character(grids))

# check that all canape non-sig grid cells are ok
emp_canape_res %>% 
  select("site", "endem_type", "pe_obs_p_upper", "pe_alt_obs_p_upper") %>% 
  filter(endem_type != "not significant") %>% 
  filter(pe_obs_p_upper < 0.95 & pe_alt_obs_p_upper < 0.95)

## load results data - analyses and make tables ####
out_fold <- "fern"
canaper_res_df <- readRDS(paste0("output/canaper/canaper_all_data_res_", out_fold,".rds"))

# join with empirical df
canaper_res_all_df <- left_join(canaper_res_df, 
                                emp_df %>% dplyr::rename(emp_endem_type = endem_type,
                                                         emp_pd_signif = pd_signif, 
                                                         emp_pe_signif = pe_signif), 
                                by = join_by(site)
                                )

# endemicity categories
canaper_res_all_df %>% group_by(endem_type) %>% tally

# count by significant categories
canaper_res_all_df %>% group_by(pd_signif) %>% tally
canaper_res_all_df %>% group_by(pe_signif) %>% tally

# NAs count 
canaper_res_all_df %>% 
  group_by(treat, samp_frac, imputation) %>% 
  select(endem_type) %>% 
  tally(is.na(endem_type)) %>% filter(n > 0)

# CANAPE NA grids per pruning method and treatment
canaper_res_all_df %>% 
  group_by(treat, samp_frac, imputation) %>% 
  select(endem_type, emp_endem_type) %>%
  mutate(uncat = ifelse(is.na(endem_type), "is_na", "not_na")) %>% 
  summarise(prop_l = mean(uncat == "is_na")) %>%  
  mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed")) %>% 
  ggplot(aes(x = samp_frac, y = (prop_l), color = treat, shape = treat)) + 
    geom_point() + geom_line() + 
  labs(title = "CANAPE", 
       x = "Sampling proportion", y = "Proportion of NAs") + 
  facet_wrap(~is_imp, scales = "free") + 
  theme_classic() +
  theme(strip.background = element_blank())
    
ggsave(filename = paste0("figures/", out_fold,"/canaper_na_grids.pdf"), 
       width = 30, height = 20, units = "cm", dpi = 200)

## Estimate: Precision, recall and f-metric #####
library(yardstick)

canaper_cat_df <- canaper_res_all_df %>% 
  mutate(endem_type = gsub("super", "mixed", endem_type),
         emp_endem_type = gsub("super", "mixed", emp_endem_type)) %>%
  #filter(endem_type !=is.na(endem_type)) %>% 
  mutate(emp_endem_type = factor(emp_endem_type), endem_type = factor(endem_type),
         emp_pd_signif = signif2cats(emp_pd_signif) %>% factor, pd_signif = signif2cats(pd_signif) %>% factor,
         emp_pe_signif = signif2cats(emp_pe_signif) %>% factor, pe_signif = signif2cats(pe_signif) %>% factor)

# confusion matrix per sims 
conf_mats <- canaper_cat_df %>% 
  group_by(sim_n, treat, samp_frac, imputation) %>%
  conf_mat(truth = emp_endem_type, estimate = endem_type, na_rm = T) 

conf_mats$conf_mat[[1]] %>% 
  ggplot2::autoplot(type = "heatmap") + # "mosaic"
  geom_label(aes(label = Freq), fill = "white")

# predictive stats
# precision, recall, f-score, Matthews correlation coeff
canaper_stat <- canaper_cat_df %>% 
  group_by(sim_n, treat, samp_frac, imputation) %>% 
  summarize(precision = precision_vec(truth = emp_endem_type, estimate = endem_type),
            recall = recall_vec(truth = emp_endem_type, estimate = endem_type),
            f_mes = f_meas_vec(truth = emp_endem_type, estimate = endem_type),
            m_ccr = mcc_vec(truth = emp_endem_type, estimate = endem_type), # endemicity type
            
            pd_mcc = mcc_vec(truth = emp_pd_signif, estimate = pd_signif), # significance
            pe_mcc = mcc_vec(truth = emp_pe_signif, estimate = pe_signif)
            )

canaper_mean_stats <- canaper_stat %>% 
  ungroup() %>% 
  group_by(samp_frac, treat, imputation) %>% 
  summarise(mean_prec = mean(precision, na.rm = T),
            sd_prec = sd(precision, na.rm = T),
            mean_reca = mean(recall, na.rm = T),
            sd_reca = sd(recall, na.rm = T),
            mean_fmes = mean(f_mes, na.rm = T),
            sd_fmes = sd(f_mes, na.rm = T),
            
            mean_mcc = mean(m_ccr, na.rm = T),
            sd_mcc = sd(m_ccr, na.rm = T),
            
            m_pd_mcc = mean(pd_mcc, na.rm = T),
            sd_pd_mcc = sd(pd_mcc, na.rm = T),
            m_pe_mcc = mean(pe_mcc, na.rm = T),
            sd_pe_mcc = sd(pe_mcc, na.rm = T)
            ) 

# add grouping columns
canaper_mean_stats <- canaper_mean_stats %>%
  mutate(prun_type = case_when(str_detect(treat, "late-pr") ~ "Lateral - Phylogenetic", 
                               str_detect(treat, "late-spa") ~ "Lateral - Spatial",
                               str_detect(treat, "marg-pr") ~ "Bilateral - Phylogenetic", 
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
         prun_mhtd = case_when(str_detect(treat, "late-pr") ~ "Phylogenetic", 
                               str_detect(treat, "late-spa") ~ "Spatial",
                               str_detect(treat, "marg-pr") ~ "Phylogenetic", 
                               str_detect(treat, "marg-spa") ~ "Spatial",
                               str_detect(treat, "rand-pr") ~ "Random", 
                               TRUE ~ "other" # Default value
                               )
         ) %>% 
  mutate(is_imp = ifelse(imputation == "none", "No imputed", "Imputed"),
         is_rand = ifelse(treat == "rand-pr", "Random pruning", "Biased pruning"),
         cat2map = paste0(treat, "-" , imputation)) %>% 
  mutate(prun_type = factor(prun_type, levels = c("Random", "Trait", "Spatial"), labels = c("Random", "Phylogenetic", "Spatial")),
         treat = factor(treat, levels = c("rand-pr", "late-pr", "marg-pr", "late-spat_pr", "marg-spat_pr"), 
                        labels = c("Random", "Lateral - Phylogenetic", "Bilateral - Phylogenetic", "Lateral - Spatial", "Bilateral - Spatial")),
         imputation = factor(imputation, levels = c("none", "random", 'taxonomy - node', "taxonomy - random"), 
                             labels = c("No imputation", "Random", "At crown node", "Bellow crown node")) )
  
## figs #####
dir.create(paste0("figures/", out_fold), recursive = T)

forms2print <- c(".png", ".pdf", ".svg") 
dpi2plot <- 200

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
  
  ### MCC ####
  # in one fig
  canaper_mean_stats %>% 
    #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
    ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    theme_classic(base_size = 12) + 
    coord_cartesian(ylim = c(0.5, 0.98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_mcc_sd", fig_out), 
         width = 20, height = 20, units = "cm", dpi = dpi2plot)
  
  # without random imputation in one fig
  canaper_mean_stats %>% 
    filter(imputation != "Random") %>% 
    #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
    ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    theme_classic(base_size = 12) + 
    coord_cartesian(ylim = c(0.5, 0.98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_mcc_sd-rand_imp", fig_out), 
         width = 20, height = 20, units = "cm", dpi = dpi2plot)
  
  # by pruning type
  canaper_mean_stats %>% 
    #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
    ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    theme_classic(base_size = 12) + 
    coord_cartesian(ylim = c(0.5, 0.98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), ncol = 3) + 
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_mcc_sd_prun3type", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # by pruning type without random imputation
  canaper_mean_stats %>% 
    filter(imputation != "Random") %>% 
    #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
    ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    theme_classic(base_size = 12) + 
    coord_cartesian(ylim = c(0.5, 0.98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), ncol = 3) + 
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_mcc_sd_prun3type-rand_imp", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # by pruning method
  canape_p <- canaper_mean_stats %>% 
    #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
    ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    theme_classic(base_size = 12) + 
    coord_cartesian(ylim = c(0.5, 0.98)) +
    labs(title = "CANAPE", 
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
  ggsave(plot = canape_p, filename = paste0("figures/", out_fold, "/canaper_mcc_sd_prun_mhtd", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # by pruning method without random imputation
  canaper_mean_stats %>% 
    filter(imputation != "Random") %>% 
    #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
    ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    theme_classic(base_size = 12) + 
    coord_cartesian(ylim = c(0.5, 0.98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_mcc_sd_prun_mthd-rand_imp", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # by if imputed
  canaper_mean_stats %>% 
    #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
    ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    theme_classic(base_size = 12) + 
    coord_cartesian(ylim = c(0.5, 0.98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(is_imp, c("No imputed", "Imputed")), labeller = label_value, ncol = 3) +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15),strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_mcc_sd_if-imput", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # by imputing method
  canaper_mean_stats %>% 
    #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
    ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    theme_classic(base_size = 12) + 
    coord_cartesian(ylim = c(0.5, 0.98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(imputation), labeller = label_value, ncol = 2, scales = "fixed") +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
    ggsave(filename = paste0("figures/", out_fold, "/canaper_mcc_sd_imputation", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
    
    # panel canaper, PD and PE significance levels
    pd_sig_p <- canaper_mean_stats %>% 
      #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
      ggplot(aes(x = samp_frac, 
                 y = m_pd_mcc, ymin = m_pd_mcc - sd_pd_mcc, ymax = m_pd_mcc + sd_pd_mcc,
                 color = imputation, shape = treat)) +
      geom_line(linewidth = 0.9, alpha = 0.65) + 
      geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
      scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
      scale_shape_manual(values = c("Random" = 4, 
                                    "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                    "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
      coord_cartesian(ylim = c(0.05, 0.98)) +
      labs(title = "PD significance", 
        x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
        color = "Imputation", shape = "Pruning\nmethod") +
      facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
      guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
             shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
      theme_classic(base_size = 12) +
      theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")
    
    pe_sig_p <- canaper_mean_stats %>% 
      #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
      ggplot(aes(x = samp_frac, 
                 y = m_pe_mcc, 
                 ymin = m_pe_mcc - sd_pe_mcc, ymax = m_pe_mcc + sd_pe_mcc,
                 color = imputation, shape = treat)) +
      geom_line(linewidth = 0.9, alpha = 0.65) + 
      geom_pointrange(position = position_dodge2(width = .02, padding = 0.5), linewidth = 0.3, fatten = 3) + 
      scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
      scale_shape_manual(values = c("Random" = 4, 
                                    "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                    "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
      coord_cartesian(ylim = c(0.45, 0.98)) +
      labs(title = "PE significance", 
           x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
           color = "Imputation", shape = "Pruning\nmethod") +
      facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
      guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
             shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
      theme_classic(base_size = 12) +
      theme(strip.text = element_blank(), strip.background = element_blank(), legend.position = "right")
    
    pl_long <- "
  AAAA
  AAAA
  BBBB
  BBBB
  CCCC
  CCCC
  "
    
    canape_p + pd_sig_p + pe_sig_p +  
      plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new")
    ggsave(filename = paste0("figures/", out_fold, "/mcc_sd_signif_pnl", fig_out), 
           width = 30, height = 40, units = "cm", dpi = 200)
    
    # MCC sd trends 
    canaper_mean_stats %>% 
      filter(imputation != "random") %>% 
      mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
             cat2map = paste0(treat, "-" , imputation)) %>% 
      select(samp_frac, treat, imputation, contains("sd_")) %>% 
      select(samp_frac, treat, imputation, contains("_mcc")) %>% 
      gather(key, value, -samp_frac, -treat, -imputation) %>% 
      ggplot(aes(x = samp_frac, y = value, 
                 color = imputation, shape = treat)) +
      geom_line() + geom_point() +
      theme_classic(base_size = 12) + 
      theme(strip.background = element_blank(), legend.position = "right") +
      labs(x = "Sampling fraction", y = "MCC's SD") +
      scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
      facet_wrap(~ key, scales = "fixed", ncol = 3)
    
    ggsave(filename = paste0("figures/", out_fold, "/canaper_mcc_metric_sd_pnl", fig_out), 
           width = 30, height = 40, units = "cm", dpi = 200)
  
  ## F measure ####
    # mean and sd in one fig
    canaper_mean_stats %>% 
      ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
                 color = imputation, shape = treat)) +
      geom_line(linewidth = 0.9, alpha = 0.65) + 
      geom_pointrange(position = position_dodge2(width = .02, padding = 0.1), linewidth = 0.3, fatten = 3) + 
      scale_shape_manual(values = c("Random" = 4, 
                                    "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                    "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
      scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
      ylim(c(0.35, .98)) +
      labs(#title = "Endemism categorization", 
        x = "Proportion of tips removed", y = "Mean F-score", 
        color = "Imputation", shape = "Pruning\nmethod") + 
      guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
             shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
      theme_classic(base_size = 12) +
      theme(strip.background = element_blank(), legend.position = "right")
    
    ggsave(filename = paste0("figures/", out_fold, "/canaper_fmeas_sd", fig_out), 
           width = 20, height = 20, units = "cm", dpi = dpi2plot)
    
  # mean and sd in one fig
  canaper_mean_stats %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.1), linewidth = 0.3, fatten = 3) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    ylim(c(0.35, .98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean F-score", 
      color = "Imputation", shape = "Pruning\nmethod") + 
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_fmeas_sd", fig_out), 
         width = 20, height = 20, units = "cm", dpi = dpi2plot)  
  
  # only mean values 
  canaper_mean_stats %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, color = imputation, shape = treat)) +
      geom_line(linewidth = 0.9, alpha = 0.65) + 
      geom_point(size = 2.5) + 
      scale_shape_manual(values = c("Random" = 4, 
                                    "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                    "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
      scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    theme_classic(base_size = 12) + 
    ylim(c(0.35, .98)) +
      labs(#title = "Endemism type", 
           x = "Proportion of tips removed", y = "Mean F-score", 
           color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(is_imp, c("No imputed", "Imputed")), ncol = 2, scales = "fixed") +
      theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")

  ggsave(filename = paste0("figures/", out_fold, "/canaper_fmeas_if-imput", fig_out), 
         width = 25, height = 15, units = "cm", dpi = dpi2plot)

  # mean and sd by treatment
  canaper_mean_stats %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
               color = imputation, shape = treat)) + 
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.1), linewidth = 0.3, fatten = 3) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    theme_classic(base_size = 12) + 
    ylim(c(0.35, .98)) +
    labs(#title = "Endemism type", 
         x = "Proportion of tips removed", y = "Mean F-score", 
         color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(is_imp, c("No imputed", "Imputed")), ncol = 2, scales = "fixed") +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_fmeas_sd_if-imput", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # with sd by pruning
  canaper_mean_stats %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.1), linewidth = 0.3, fatten = 3) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    ylim(c(0.35, .98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean F-score", 
      color = "Imputation", shape = "Pruning\nmethod") + 
      facet_wrap(~factor(is_rand, c("Random pruning", "Biased pruning")), ncol = 2, scales = "free_y") + 
      guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
             shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
      theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
    
  ggsave(filename = paste0("figures/", out_fold, "/canaper_fmeas_sd_pruning", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # with sd without random imputation by if imputed
  canaper_mean_stats %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.1), linewidth = 0.3, fatten = 3) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    ylim(c(0.35, .98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean F-score", 
      color = "Imputation", shape = "Pruning\nmethod") + 
    facet_wrap(~factor(is_imp, c("No imputed", "Imputed")), ncol = 2, scales = "free_y") + 
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_fmeas_sd-rand_imp_if-imput", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # with sd by pruning method 
  canaper_mean_stats %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.1), linewidth = 0.3, fatten = 3) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    ylim(c(0.35, .98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean F-score", 
      color = "Imputation", shape = "Pruning\nmethod") + 
    facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic" ,"Spatial")), ncol = 3, scales = "free_y") + 
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_fmeas_sd_pruning_mthd", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # with sd separated by pruning type (Random, Lateral, Bilateral)
  canaper_mean_stats %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.1), linewidth = 0.3, fatten = 3) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    ylim(c(0.35, .98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean F-score", 
      color = "Imputation", shape = "Pruning\nmethod") + 
    facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), ncol = 3) +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_fmeas_sd_prun3type", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
  
  # with sd but separated by pruning type (Random, Lateral, Bilateral), without random imputation
  canaper_mean_stats %>% 
    filter(imputation != "random") %>% 
    ggplot(aes(x = samp_frac, y = mean_fmes, ymin = mean_fmes - sd_fmes, ymax = mean_fmes + sd_fmes,
               color = imputation, shape = treat)) +
    geom_line(linewidth = 0.9, alpha = 0.65) + 
    geom_pointrange(position = position_dodge2(width = .02, padding = 0.1), linewidth = 0.3, fatten = 3) + 
    scale_shape_manual(values = c("Random" = 4, 
                                  "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                  "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
    scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
    ylim(c(0.35, .98)) +
    labs(#title = "Endemism categorization", 
      x = "Proportion of tips removed", y = "Mean F-score", 
      color = "Imputation", shape = "Pruning\nmethod") + 
    facet_wrap(~factor(prun3type, c("Random", "Bilateral", "Lateral")), ncol = 3) +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right")
  
  ggsave(filename = paste0("figures/", out_fold, "/canaper_fmeas_sd_prun3type-rand_imp", fig_out), 
         width = 30, height = 15, units = "cm", dpi = dpi2plot)
}

## animated #####
library(gganimate)

# by pruning method
mcc_meas_pt <- canaper_mean_stats %>% 
  filter(imputation != "Random") %>% 
  #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
  ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
             color = imputation, shape = treat)) +
  geom_line(aes(group = interaction(prun_mhtd, treat, imputation)), linewidth = 0.9, alpha = 0.65) +
  geom_errorbar(aes(ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc, group = interaction(prun_mhtd, treat, imputation, samp_frac)),
                width = 0, linewidth = 0.3, alpha = 0.7) +
  geom_point(aes(group = interaction(prun_mhtd, treat, imputation, samp_frac)), size = 3) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  coord_cartesian(ylim = c(0.5, 0.98)) +
  labs(subtitle = "CANAPE", 
    x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
    color = "Imputation", shape = "Pruning\nmethod") +
  facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(size = 15), strip.background = element_blank(), legend.position = "right") +
  transition_reveal(samp_frac) + shadow_mark(alpha = 0.3, size = 0.5) 

animate(mcc_meas_pt, width = 40, height = 20, units = "cm", 
        res = 200, nframes = 60,  fps = 10, duration = 20, 
        start_pause = 2, end_pause = 60)

anim_save(filename = paste0("figures/", out_fold, "/canaper_mcc_anim_pnl.gif"))

# in one fig
f_meas_pt2 <- canaper_mean_stats %>% 
  filter(imputation != "Random") %>% 
  #mutate(across(starts_with("sd"), ~replace_na(.x, 0))) %>%
  ggplot(aes(x = samp_frac, y = mean_mcc, ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc,
             color = imputation, shape = treat)) +
  geom_line(aes(group = interaction(prun_mhtd, treat, imputation)), linewidth = 0.9, alpha = 0.65) +
  geom_errorbar(aes(ymin = mean_mcc - sd_mcc, ymax = mean_mcc + sd_mcc, group = interaction(prun_mhtd, treat, imputation, samp_frac)),
    width = 0, linewidth = 0.3, alpha = 0.7) +
  geom_point(aes(group = interaction(prun_mhtd, treat, imputation, samp_frac)), size = 3) +
  scale_colour_viridis_d(option = "inferno", begin = .8, end = .05, alpha = .95) + 
  scale_shape_manual(values = c("Random" = 4, 
                                "Lateral - Phylogenetic" = 20,  "Bilateral - Phylogenetic" = 18, 
                                "Lateral - Spatial" = 15, "Bilateral - Spatial" = 17)) +
  coord_cartesian(ylim = c(0.5, 0.98)) +
  labs(subtitle = "CANAPE", 
    x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
    color = "Imputation", shape = "Pruning\nmethod") +
  #facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
  guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
         shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right") +
  transition_reveal(samp_frac) + shadow_mark(alpha = 0.3, size = 0.5)  + 
  ease_aes('sine-in-out')
animate(f_meas_pt2, width = 25, height = 25, units = "cm", 
        res = 200, nframes = 60,  fps = 10, duration = 10, 
        start_pause = 2, end_pause = 60)
anim_save(filename = paste0("figures/", out_fold, "/canaper_mcc_anim.gif"))

## descriptive values ####
canaper_mean_stats %>% 
  dplyr::select(samp_frac, treat, imputation, mean_mcc) %>% # mean
  gather(key, value, -samp_frac, -treat, -imputation) %>% 
  filter(#treat == "rand-pr" #, 
         imputation != "none" , samp_frac == 0.1
  ) %>% view
  group_by(key, samp_frac, treat,imputation) %>% 
  summarise(min = min(value, na.rm = T) %>% round(2), 
            max = max(value, na.rm = T) %>% round(2)) %>% 
  filter(samp_frac == 0.5) %>% 
    view

# figure p-values by endemicity type ####
cols2sel <- c("site", "endem_type",
              "pe_signif", "rpe_signif",
              "pe_obs_p_upper", "pe_alt_obs_p_upper")
  
## full dataset ####
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))
emp_df <- emp_canape_res %>% select(cols2sel) ; emp_df %>% head()

x1end <- 50 ; x2start <- 75

emp_pl <- emp_df %>% 
    #filter(endem_type != "not significant" #, pe_signif == "not significant") %>% 
    gather(key, value, -site, -endem_type, -pe_signif, -rpe_signif) %>% 
    ggplot(aes(x = value, color = endem_type)) +
    geom_freqpoly() + 
    ylim(0, x2start) +
    #geom_histogram(alpha = .5) +
    #geom_density(alpha = .5) + 
    #theme_minimal() + 
    geom_vline(xintercept = 0.95, colour = "steelblue", lwd = .5, lty = 2) +
    labs(title = "Full tree", x = "p-value", y = "grid cell count", color = "CANAPE") +
    scale_color_viridis_d(option = "plasma", begin = 0, end = .8) + 
    facet_wrap(~ key, scales = "fixed", ncol = 1)

## treated datasets #### 
cols2sel <- c("site", "endem_type", 
              "pe_signif", #"rpe_signif", 
              "pe_obs_p_upper", "pe_alt_obs_p_upper")

out_fold <- "fern"
canaper_res_df <- readRDS(paste0("output/canaper/canaper_all_data_res_", out_fold,".rds"))

late_prun_pl <- canaper_res_df %>%
    select("site", "endem_type", 
           "imputation", "samp_frac", "treat",
           "pe_signif", 
           "pe_obs_p_upper", "pe_alt_obs_p_upper") %>% 
    filter(imputation == "none", 
           samp_frac == 0.5,
           treat == "late-spat_pr") %>% 
    filter(!is.na(endem_type)) %>% 
    droplevels() %>% 
    #filter(endem_type != "not significant"#, pe_signif == "not significant") %>% 
    gather(key, value, -site, -endem_type, 
           -imputation, -samp_frac, -treat, 
           -pe_signif) %>% 
    ggplot(aes(x = value, color = endem_type)) +
    #geom_histogram(alpha = .5) +
    #geom_density(alpha = .5) + 
    geom_freqpoly() +
    ylim(0, x2start) +
    geom_vline(xintercept = 0.95, colour = "steelblue", lwd = .5, lty = 2) +
    #theme_minimal() + 
    labs(title = "Lateral - Spatial pruned (50%)", x = "p-value", y = "grid cell count", color = "CANAPE") +
    scale_color_viridis_d(option = "plasma", begin = 0, end = .8) + 
    facet_wrap(~ key, scales = "fixed", ncol = 1)
  
rand_prun_pl <- canaper_res_df %>%
    select("site", "endem_type", "sim_n",
           "imputation", "samp_frac", "treat",
           "pe_signif", 
           "pe_obs_p_upper", "pe_alt_obs_p_upper") %>%
    filter(imputation == "none", 
           samp_frac == 0.5,
           treat == "rand-pr", #sim_n == "sim_1"
    ) %>% 
    filter(!is.na(endem_type)) %>% 
    #filter(endem_type != "not significant"#, pe_signif == "not significant") %>% 
    gather(key, value, -site, -endem_type, 
           -imputation, -samp_frac, -treat, -sim_n,
           -pe_signif) %>% 
    ggplot(aes(x = value, ..count../ 100,
               color = endem_type)) +
    #geom_histogram(alpha = .5) +
    #geom_density(alpha = .5) + 
    geom_freqpoly() +
    ylim(0, x2start) +
    geom_vline(xintercept = 0.95, colour = "steelblue", lwd = .5, lty = 2) +
    #theme_minimal() + 
    labs(title = "Random pruned (50%)", x = "p-value", y = "grid cell count", color = "CANAPE") +
    scale_color_viridis_d(option = "plasma", begin = 0, end = .8) + 
    facet_wrap(~ key, scales = "fixed", ncol = 1)
  
## panel plot 
emp_pl + late_prun_pl + rand_prun_pl + 
    plot_layout(guides = "collect", axis_titles = "collect", tag_level = "new")
ggsave(paste0("figures/", out_fold,"/canape_pe_pval_pnl", ".pdf"), width = 25, height = 20, units = "cm", dpi = 200)

