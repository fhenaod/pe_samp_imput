library(tidyverse)

# load and get canaper tables ####
# empirical full data
cols2sel <- c("site", 
              "pd_obs", "pe_obs", "rpd_obs", "rpe_obs", 
              "pd_obs_z", "pe_obs_z", "rpd_obs_z", "rpe_obs_z"
              #"pd_rand_mean","pe_rand_mean", "rpd_rand_mean", "rpe_rand_mean"
              )

canaper_emp_raw <- readRDS("output/canaper/canaper_empirical.rds")
canaper_emp <- canaper_emp_raw %>% 
  select(any_of(cols2sel)) %>% 
  mutate(treat = "empirical", samp_frac = 1, imputation = "none") %>% 
  mutate(rw_nam = paste0(site,"_", samp_frac, "_" , treat, "_", imputation)) %>% 
  column_to_rownames(.data = ., var = "rw_nam")
head(canaper_emp)

# results data - analyses
canaper_res_all <- readRDS("output/canaper_all_data_res.rds")

# get mean across 100 sims per samp fraction, pruning and imputation methods
canaper_res_df <- canaper_res_all %>% 
  mutate(treat = gsub("_", "-", treat)) %>%
  group_by(site, samp_frac, treat, imputation) %>%
  summarise(across(pd_obs:rpe_obs_z, mean, .names = "{.col}")) %>% 
  mutate(rw_nam = paste0(site,"_", samp_frac, "_" , treat, "_", imputation)) %>% 
  column_to_rownames(var = "rw_nam") %>% select(names(canaper_emp))
head(canaper_res_df)  

# pca by sampling fraction and colored by treatment ####
library(factoextra) ; library(ggpubr)

samps_fr <- canaper_res_df %>% group_by(samp_frac) %>% tally() %>% pull(samp_frac) ; samps_fr <- samps_fr[-2]
canaper_res_df %>% group_by(treat) %>% tally()
canaper_res_df %>% group_by(imputation) %>% tally()
impt <- list(c("none"), c("random", "taxonomy - node", "taxonomy - random"))

# loop over sampling fractions and biplot
biplot_ls <- list()
for(j in 1:length(impt)){
  pca_res_ls <- list() ; pca_plt_ls <- list()
  for(i in 1:length(samps_fr)){
  dat2pca <- rbind(canaper_emp %>% select(-c(site, treat, samp_frac, imputation)), 
                   canaper_res_df %>% ungroup() %>% 
                     filter(samp_frac == samps_fr[[i]], 
                            imputation %in% impt[[j]]) %>% 
                       filter(!is.na(rpd_obs)) %>% select(-c(site, samp_frac, treat, imputation))
                   )
  
  pca_res_ls[[i]] <- dat2pca %>% prcomp(., scale = T)
  
  meta_cat <- pca_res_ls[[i]]$x %>% data.frame() %>%
    rownames_to_column(.data = ., var = "rnw") %>% 
    separate_wider_delim(cols = rnw, 
                         names = c("site", "samp_frac", "treat", "imputation"), 
                         delim = "_") %>% 
    mutate(cat2map = paste0(treat, "-" , imputation)) %>% pull(cat2map)
  
  pca_plt_ls[[i]] <- pca_res_ls[[i]] %>% 
    fviz_pca_biplot(., geom.ind = "point", #c("point", "text", "arrows"),
                    #col.var = "cos2", 
                    col.ind = meta_cat, 
                    label = "none", #label = "var", labelsize = .35, 
                    pointshape = 21, pointsize = .6, arrowsize = 0, 
                    alpha = .35,
                    #gradient.cols = c("steelblue", "orange", "red"),
                    addEllipses = T, ellipse.level = 0.95, 
                    #gradient.cols = "npg", 
                    #select.ind = list(contrib = 30),
                    title = paste0("Removed proportion: ", samps_fr[[i]]),
                    subtitle = ifelse(length(impt[[j]]) > 1, "Imputation", "No imputation") , 
                    legend.title = "", 
                    repel = T) + 
    scale_color_brewer(palette = "Set1") +
    theme_minimal(base_size = 11) + 
    theme(legend.position = "top", 
          legend.title = element_blank())
  } 
  biplot_ls[[j]] <- pca_plt_ls
}


# one biplot panel
ggarrange(plotlist = do.call(c, biplot_ls), nrow = 2, ncol = 5, common.legend = T)
ggsave(filename = paste0("figures/pca_res.pdf"), 
       width = 40, height = 20, units = "cm", dpi = 200)

# separated panels
ggarrange(plotlist = biplot_ls[[1]], nrow = 1, common.legend = T)
ggarrange(plotlist = biplot_ls[[2]], nrow = 1, common.legend = T)
