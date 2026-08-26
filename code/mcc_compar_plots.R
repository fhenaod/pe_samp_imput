library(tidyverse)
library(yardstick)
library(patchwork)

# load functions
source("code/pe_samp_funs.R")

# load results ####
cols2sel <- c("site", "endem_type",
              "pd_signif", "pe_signif"
              #, "rpd_signif", "rpe_signif"
              )

# from Nitta et al. 2022 dataset
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))
emp_df <- emp_canape_res %>% select(cols2sel) ; emp_df %>% head()

# names files to load
names2load <- c("fern", "old", "cnt", "young", "rec", "exp")

# names for subtitles
nams4til <- c("fern", "300 my ± 25", "217 my ± 25", "100 my ± 25" , "20 my ± 25", "Exponential")

forms2print <- c(".png", ".pdf", ".svg") 

# horizontal arrangement ####
# loop over files creating a mcc panel plot
mcc_plt_ls_h <- list()
for(h in seq_along(names2load)) {
  out_fold <- names2load[h]
  
  canaper_res_df <- readRDS(paste0("output/canaper/canaper_all_data_res_", out_fold, ".rds"))
  print(paste0(toupper(out_fold), " loaded"))
  
  # join with empirical df
  canaper_res_all_df <- left_join(canaper_res_df, 
                                  emp_df %>% rename(emp_endem_type = endem_type,
                                                    emp_pd_signif = pd_signif,
                                                    emp_pe_signif = pe_signif), 
                                  by = join_by(site))
  
  
  # prepare dataframe
  canaper_cat_df <- canaper_res_all_df %>% 
    mutate(endem_type = gsub("super", "mixed", endem_type),
           emp_endem_type = gsub("super", "mixed", emp_endem_type)) %>%
    #filter(endem_type !=is.na(endem_type)) %>% 
    mutate(emp_endem_type = factor(emp_endem_type), endem_type = factor(endem_type),
           emp_pd_signif = signif2cats(emp_pd_signif) %>% factor, pd_signif = signif2cats(pd_signif) %>% factor,
           emp_pe_signif = signif2cats(emp_pe_signif) %>% factor, pe_signif = signif2cats(pe_signif) %>% factor)
  
  #  Matthews correlation coeff
  canaper_stat <- canaper_cat_df %>% 
    group_by(sim_n, treat, samp_frac, imputation) %>% 
    summarize(#precision = precision_vec(truth = emp_endem_type, estimate = endem_type),
      #recall = recall_vec(truth = emp_endem_type, estimate = endem_type),
      #f_mes = f_meas_vec(truth = emp_endem_type, estimate = endem_type),
      m_ccr = mcc_vec(truth = emp_endem_type, estimate = endem_type)#,
      
      #pd_mcc = mcc_vec(truth = emp_pd_signif, estimate = pd_signif),
      #pe_mcc = mcc_vec(truth = emp_pe_signif, estimate = pe_signif)
    )
  
  # get summary stats
  canaper_mean_stats <- canaper_stat %>% 
    ungroup() %>% 
    group_by(samp_frac, treat, imputation) %>% 
    summarise(#mean_prec = mean(precision, na.rm = T),
      #sd_prec = sd(precision, na.rm = T),
      #mean_reca = mean(recall, na.rm = T),
      #sd_reca = sd(recall, na.rm = T),
      #mean_fmes = mean(f_mes, na.rm = T),
      #sd_fmes = sd(f_mes, na.rm = T),
      mean_mcc = mean(m_ccr, na.rm = T),
      sd_mcc = sd(m_ccr, na.rm = T)#,
      
      #m_pd_mcc = mean(pd_mcc, na.rm = T),
      #sd_pd_mcc = sd(pd_mcc, na.rm = T),
      #m_pe_mcc = mean(pe_mcc, na.rm = T),
      #sd_pe_mcc = sd(pe_mcc, na.rm = T)
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
    mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
           is_rand = ifelse(treat == "rand-pr", "Random pruning", "Biased pruning"),
           cat2map = paste0(treat, "-" , imputation)) %>% 
    mutate(prun_type = factor(prun_type, levels = c("Random", "Trait", "Spatial"), labels = c("Random", "Phylogenetic", "Spatial")),
           treat = factor(treat, levels = c("rand-pr", "late-pr", "marg-pr", "late-spat_pr", "marg-spat_pr"), 
                          labels = c("Random", "Lateral - Phylogenetic", "Bilateral - Phylogenetic", "Lateral - Spatial", "Bilateral - Spatial")),
           imputation = factor(imputation, levels = c("none", "random", 'taxonomy - node', "taxonomy - random"), 
                               labels = c("No imputation", "Random", "At crown node", "Below crown node")) )
  
  # plot mcc by pruning method
  mcc_plt_ls_h[[h]] <- canaper_mean_stats %>% 
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
    coord_cartesian(ylim = c(0.3, 0.98)) +
    labs(title = nams4til[h],
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, ncol = 3) +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), legend.position = "right")
}
names(mcc_plt_ls_h) <- names2load
length(mcc_plt_ls_h)

# panel: empirical and all sim trees - horizontal
pl_long <- "
  AAAA
  BBBB
  CCCC
  DDDD
  EEEE
  "

(mcc_plt_ls_h$old) + 
(mcc_plt_ls_h$cnt) +
(mcc_plt_ls_h$young) +
(mcc_plt_ls_h$rec) +
(mcc_plt_ls_h$exp) + 
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new") +
  plot_annotation(title = "CANAPE", theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0)))

for(i in 1:length(forms2print)){
  fig_out <- forms2print[i]
ggsave(filename = paste0("figures/mcc_sim_tress_all_pnl_h", fig_out), 
       width = 30-10, height = 45-10, units = "cm", dpi = 200)
}

# vertical arrangement #####
# loop over files creating a mcc panel plot
mcc_plt_ls_v <- list()
for(h in seq_along(names2load)) {
  out_fold <- names2load[h]
  
  canaper_res_df <- readRDS(paste0("output/canaper/canaper_all_data_res_", out_fold, ".rds"))
  print(paste0(toupper(out_fold), " loaded"))
  
  # join with empirical df
  canaper_res_all_df <- left_join(canaper_res_df, 
                                  emp_df %>% rename(emp_endem_type = endem_type,
                                                    emp_pd_signif = pd_signif,
                                                    emp_pe_signif = pe_signif), 
                                  by = join_by(site))
  
  # prepare dataframe
  canaper_cat_df <- canaper_res_all_df %>% 
    mutate(endem_type = gsub("super", "mixed", endem_type),
           emp_endem_type = gsub("super", "mixed", emp_endem_type)) %>%
    #filter(endem_type !=is.na(endem_type)) %>% 
    mutate(emp_endem_type = factor(emp_endem_type), endem_type = factor(endem_type),
           emp_pd_signif = signif2cats(emp_pd_signif) %>% factor, pd_signif = signif2cats(pd_signif) %>% factor,
           emp_pe_signif = signif2cats(emp_pe_signif) %>% factor, pe_signif = signif2cats(pe_signif) %>% factor)
  
  #  Matthews correlation coeff
  canaper_stat <- canaper_cat_df %>% 
    group_by(sim_n, treat, samp_frac, imputation) %>% 
    summarize(#precision = precision_vec(truth = emp_endem_type, estimate = endem_type),
      #recall = recall_vec(truth = emp_endem_type, estimate = endem_type),
      #f_mes = f_meas_vec(truth = emp_endem_type, estimate = endem_type),
      m_ccr = mcc_vec(truth = emp_endem_type, estimate = endem_type)#,
      
      #pd_mcc = mcc_vec(truth = emp_pd_signif, estimate = pd_signif),
      #pe_mcc = mcc_vec(truth = emp_pe_signif, estimate = pe_signif)
    )
  
  # get summary stats
  canaper_mean_stats <- canaper_stat %>% 
    ungroup() %>% 
    group_by(samp_frac, treat, imputation) %>% 
    summarise(#mean_prec = mean(precision, na.rm = T),
      #sd_prec = sd(precision, na.rm = T),
      #mean_reca = mean(recall, na.rm = T),
      #sd_reca = sd(recall, na.rm = T),
      #mean_fmes = mean(f_mes, na.rm = T),
      #sd_fmes = sd(f_mes, na.rm = T),
      mean_mcc = mean(m_ccr, na.rm = T),
      sd_mcc = sd(m_ccr, na.rm = T)#,
      
      #m_pd_mcc = mean(pd_mcc, na.rm = T),
      #sd_pd_mcc = sd(pd_mcc, na.rm = T),
      #m_pe_mcc = mean(pe_mcc, na.rm = T),
      #sd_pe_mcc = sd(pe_mcc, na.rm = T)
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
    mutate(is_imp = ifelse(imputation == "none", "Only pruned", "Imputed"),
           is_rand = ifelse(treat == "rand-pr", "Random pruning", "Biased pruning"),
           cat2map = paste0(treat, "-" , imputation)) %>% 
    mutate(prun_type = factor(prun_type, levels = c("Random", "Trait", "Spatial"), labels = c("Random", "Phylogenetic", "Spatial")),
           treat = factor(treat, levels = c("rand-pr", "late-pr", "marg-pr", "late-spat_pr", "marg-spat_pr"), 
                          labels = c("Random", "Lateral - Phylogenetic", "Bilateral - Phylogenetic", "Lateral - Spatial", "Bilateral - Spatial")),
           imputation = factor(imputation, levels = c("none", "random", 'taxonomy - node', "taxonomy - random"), 
                               labels = c("No imputation", "Random", "At crown node", "Below crown node")) )
  
  # plot mcc by pruning method
  mcc_plt_ls_v[[h]] <- canaper_mean_stats %>% 
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
    coord_cartesian(ylim = c(0.5, 0.98)) +
    labs(#title = "CANAPE", 
      x = "Proportion of tips removed", y = "Mean Matthews correlation coefficient", 
      color = "Imputation", shape = "Pruning\nmethod") +
    facet_wrap(~factor(prun_mhtd, c("Random", "Phylogenetic", "Spatial")), labeller = label_value, nrow = 3) +
    guides(colour = guide_legend(override.aes = list(shape = "\u2013", linetype =  0, size = 3)),
           shape = guide_legend(override.aes = list(shape = c(4, 20, 18, 15, 17), linetype =  0))) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), legend.position = "right")
}
names(mcc_plt_ls_v) <- names2load
length(mcc_plt_ls_v)

# panel: empirical and all sim trees - vertical
pl_long <- "
  ABCDE
  ABCDE
  ABCDE
  ABCDE
  ABCDE
  "

(mcc_plt_ls_v$old) + 
  (mcc_plt_ls_v$cnt) +
  (mcc_plt_ls_v$young) +
  (mcc_plt_ls_v$rec) +
  (mcc_plt_ls_v$exp) + 
  plot_layout(guides = "collect", design = pl_long, axis_titles = "collect", tag_level = "new") +
  plot_annotation(title = "CANAPE", theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0)))

for(i in 1:length(forms2print)){
fig_out <- forms2print[i]
ggsave(filename = paste0("figures/mcc_sim_tress_all_pnl_v", fig_out), 
       width = 45, height = 30, units = "cm", dpi = 200)
}
