library(tidyverse)
library(canaper)
library(yardstick)

cols2sel <- c("site", "endem_type", "pe_obs_p_upper", "pe_alt_obs_p_upper")

# load results using Nitta et al. 2022 dataset ####
emp_canape_res <- readRDS(paste0("output/canaper/canaper_empirical.rds"))
emp_df <- emp_canape_res %>% select(cols2sel) ; emp_df %>% head()

spat_df <- read.csv("data/nitta-et-al/results/japan_ferns_biodiv.csv", header = T) %>% 
  select(grids, #endem_type, 
         lat, long) %>% mutate(grids = as.character(grids))


empirical_spat_df <- left_join(emp_df, spat_df, by = join_by(site == grids)) %>% 
  mutate(endem_type = gsub("super", "mixed", endem_type))

# get japan's layout
world <- map_data("world")
japan <- world[world$region == "Japan", ]

frt_out <- ".png" #".pdf"

# full dataset ####
emp_canaper <- empirical_spat_df %>% 
  filter(!is.na(lat)) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "gray80", fill = "gray90") + 
  geom_point(aes(col = endem_type, fill = endem_type), shape = 15, size = 2) + 
  # cpr_endem_cols_4 is a CVD-friendly color palette in canaper
  scale_fill_manual( values = cpr_endem_cols_4) +
  scale_color_manual(values = cpr_endem_cols_4) +
  coord_fixed() + coord_equal() + 
  ggthemes::theme_map(base_size = 12) +
  theme_minimal() +
  theme(legend.position = c(0.15, 0.85),
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(color = "gray80", fill = "gray90"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside", # panel.background = element_rect(fill = "lightblue1", color = "NA")
        axis.text = element_text(size = 12),
        axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) +
  labs(title = "Empirical - Full dataset", 
    color = "CANAPE", fill = "CANAPE",
    x = "Longitude", y = "Latitude")
emp_canaper
ggsave(filename = paste0("figures/map_canaper", frt_out), 
       width = 20, height = 25, units = "cm", dpi = 200)

# results data - analyses 
canaper_res_all <- readRDS("output/canaper_all_data_res.rds")

# 50% randomly pruned ####
rnd_canaper_fls <- readRDS("output/canaper/fern/canaper_rand_pr_0.5.rds")
rand_frc <- 0.5

rnd_canape_ls <- rnd_canaper_fls
rnd_canape_df_ls <- list()
rnd_canape_df_ls <- rnd_canape_ls %>% 
    lapply(., "[", cols2sel) %>% 
    Map(cbind, ., sim_n = paste0("sim_", rep(1:length(rnd_canape_ls))), 
        samp_frac = rand_frc, treat = "rand-pr", imputation = "none") %>% 
    bind_rows()
rnd_canape_all_df <- (rnd_canape_df_ls) %>% 
  mutate(endem_type = gsub("super", "mixed", endem_type)); dim(rnd_canape_all_df)

# 50% lateral pruned ####
late_spat_canaper_fls <- readRDS("output/canaper/fern/canaper_late_spat_pr_0.5.rds")
late_frc <- 0.5

late_canape_ls <- late_spat_canaper_fls
late_canape_df_ls <- list()
  if(class(late_canape_ls) == "list"){
    late_canape_df_ls <- late_canape_ls %>% 
      lapply(., "[", cols2sel) %>%  
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(late_canape_ls))), 
          samp_frac = late_frc, treat = "late-spat_pr", imputation = "none") %>% 
      bind_rows()
  } else {
    late_canape_df_ls <- late_canape_ls %>% 
      dplyr::select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = late_frc, treat = "late-spat_pr", imputation = "none")
  }
late_spat_canape_all_df <- late_canape_df_ls %>% 
  mutate(endem_type = gsub("super", "mixed", endem_type)); dim(late_spat_canape_all_df)

# 50% lateral pruned and randtip imputed #####
canaper_randtip <- readRDS("output/canaper/canaper_late_pr_randtip_0.5.rds")
late_frc <- 0.5

load_randtip <- F
if(load_randtip){
  late_canape_ls <- canaper_randtip
  late_canape_df_ls <- list()
  if(class(canaper_randtip) == "list"){
    late_canape_df_ls <- canaper_randtip %>% 
      lapply(., "[", cols2sel) %>%  
      Map(cbind, ., sim_n = paste0("sim_", rep(1:length(canaper_randtip))), 
          samp_frac = late_frc, treat = "late-spat_pr", imputation = "randtip") %>% 
      bind_rows()
  } else {
    late_canape_df_ls <- canaper_randtip %>% 
      dplyr::select(all_of(cols2sel)) %>% 
      mutate(sim_n = "none", samp_frac = late_frc, treat = "late-spat_pr", imputation = "randtip") 
  }
  late_spat_canape_randtip_df <- (late_canape_df_ls) ; dim(late_spat_canape_all_df)
  late_spat_canape_randtip_df %>% saveRDS("output/canaper/canaper_late_pr_randtip_0.5_res.rds")
} else {
  late_spat_canape_randtip_df <- readRDS("output/canaper/canaper_late_pr_randtip_0.5_res.rds")
}
head(late_spat_canape_randtip_df); tail(late_spat_canape_randtip_df)

threshold <- 0.7

canape_cons_df <- late_spat_canape_randtip_df %>% 
  mutate(endem_type = gsub("super", "mixed", endem_type)) %>% 
  group_by(site, endem_type) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  select(-n) %>% 
  spread(endem_type, freq, fill = 0) %>% as.data.frame %>% 
  rowwise() %>% 
  mutate(endem_consensus = {
    vals <- unlist(pick(mixed, neo, `not significant`, paleo))
    hits <- names(vals)[vals >= threshold]
    if (length(hits) == 0) 'no consensus' else paste(hits, collapse = ", ")
  }) #%>% mutate(endem_consensus = as.factor(endem_consensus))

canape_cons_df %>% group_by(endem_consensus) %>% 
  summarise(n = n()) %>% 
  mutate(freq_100 = (n / sum(n))*100)

canape_cons_df <- left_join(canape_cons_df, spat_df, 
                            by = join_by(site == grids)) #%>% mutate(endem_consensus = ct_relevel(endem_consensus, 'paleo', 'mixed', 'super', 'not significant', 'no consensus'))

# map plots ####
# set pallete with 5 colors 
cpr_endem_cols_5 <- c(cpr_endem_cols_2, "#A9838F") # "#D996BB" "#CC79A7" "#E5B2CE" "#E6C2D6"
cpr_endem_cols_5[['not significant']] <- "transparent"
names(cpr_endem_cols_5)[6] <- c('no consensus')

# remove super from legend
cpr_endem_cols_5 <- cpr_endem_cols_5[names(cpr_endem_cols_5) != "super"] 
scales::show_col(cpr_endem_cols_5)

dot_size <- 1

## full dataset #####
# CANAPE
emp_canaper5 <- empirical_spat_df %>% 
  filter(!is.na(lat)) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "gray80", fill = "gray90") + 
  geom_point(aes(col = endem_type, fill = endem_type), shape = 15, size = dot_size, show.legend = T) + 
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  coord_fixed() + coord_equal() + 
  ggthemes::theme_map(base_size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom", #legend.position = c(0.15, 0.85),
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(color = "gray80", fill = "gray90"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside", # panel.background = element_rect(fill = "lightblue1", color = "NA")
        axis.text = element_text(size = 12),
        axis.title.x     = element_blank(), axis.title.y = element_blank(),
  ) +
  guides(color = guide_legend(override.aes = list(size = 4.5))) +
  labs(title = "Full dataset", subtitle = "",
       color = "CANAPE", fill = "CANAPE",
       x = "Longitude", y = "Latitude")

# PE significance
emp_pe_sig_p <- empirical_spat_df %>% 
  filter(!is.na(lat),
         endem_type != "not significant") %>%
  mutate(pe_obs_sig = case_when(pe_obs_p_upper >= .95 ~ "*", TRUE ~ ""),
         pe_alt_sig = case_when(pe_alt_obs_p_upper >= .95 ~ "*", TRUE ~ "")) %>%  
  gather(key, value, -site, -endem_type, -pe_obs_p_upper, -pe_alt_obs_p_upper, -lat, -long) %>%
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "gray80", fill = "gray90") + 
  geom_point(aes(col = endem_type, fill = endem_type), shape = 15, size = dot_size, alpha = .2, show.legend = T) + 
  geom_text(aes(label = value), size = 3, color = "black") +
  facet_wrap(~key, nrow = 2) +
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  coord_fixed() + coord_equal() + 
  guides(color = guide_legend(override.aes = list(size = 4.5))) +
  labs(title = "Full dataset", subtitle = "",
       color = "CANAPE", fill = "CANAPE",
       x = "Longitude", y = "Latitude") +
  ggthemes::theme_map(base_size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom", #legend.position = c(0.15, 0.85),
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(color = "gray80", fill = "gray90"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside", # panel.background = element_rect(fill = "lightblue1", color = "NA")
        axis.text = element_text(size = 12),
        axis.title.x     = element_blank(), axis.title.y = element_blank(),
  )

## pruned ####
# random
# CANAPE
rand_canape_prun_df <- left_join(rnd_canape_all_df, spat_df, by = join_by(site == grids))

rand_prun_canaper <- rand_canape_prun_df %>% 
  mutate(endem_type = gsub("super", "mixed", endem_type)) %>% 
  filter(sim_n == "sim_1") %>% 
  filter(!is.na(lat)) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "gray80", fill = "gray90") + 
  geom_point(aes(col = endem_type, fill = endem_type), shape = 15, size = dot_size, show.legend = T) + 
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  coord_fixed() + coord_equal() + 
  ggthemes::theme_map(base_size = 10) +
  theme_minimal() +
  theme(legend.position = "none",
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(color = "gray80", fill = "gray90"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",# panel.background = element_rect(fill = "lightblue1", color = "NA")
        #axis.ticks = element_text(color = "gray9")
        axis.text = element_blank(), #element_text(size = 12),
        axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(title = paste0("Randomly pruned (50%)"), 
       subtitle = "", 
       color = "CANAPE", fill = "CANAPE",
       x = "Longitude", y = "Latitude")

# PE significance
rand_pe_sig_p <- rand_canape_prun_df %>% 
  filter(sim_n == "sim_1") %>% 
  select(-sim_n, -samp_frac, -treat, -imputation) %>% 
  filter(!is.na(lat),
         endem_type != "not significant") %>%
  mutate(pe_obs_sig = case_when(pe_obs_p_upper >= .95 ~ "*", TRUE ~ ""),
         pe_alt_sig = case_when(pe_alt_obs_p_upper >= .95 ~ "*", TRUE ~ "")) %>% 
  gather(key, value, -site, -endem_type, -pe_obs_p_upper, -pe_alt_obs_p_upper, -lat, -long) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "gray80", fill = "gray90") + 
  geom_point(aes(col = endem_type, fill = endem_type), shape = 15, size = dot_size, alpha = .2, show.legend = T) + 
  geom_text(aes(label = value), size = 2.5, color = "black") +
  facet_wrap(~key, nrow = 2) +
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  coord_fixed() + coord_equal() + 
  guides(color = guide_legend(override.aes = list(size = 4.5))) +
  labs(title = "Random pruned (50%) ", subtitle = "",
       color = "CANAPE", fill = "CANAPE",
       x = "Longitude", y = "Latitude") +
  ggthemes::theme_map(base_size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom", #legend.position = c(0.15, 0.85),
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(color = "gray80", fill = "gray90"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside", # panel.background = element_rect(fill = "lightblue1", color = "NA")
        axis.text = element_text(size = 12),
        axis.title.x     = element_blank(), axis.title.y = element_blank(),
  )

# lateral
# CANAPE
late_canape_prun_df <- left_join(late_spat_canape_all_df, spat_df, by = join_by(site == grids))

late_prun_canaper <- late_canape_prun_df %>% 
  filter(!is.na(lat)) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "gray80", fill = "gray90") + 
  geom_point(aes(col = endem_type, fill = endem_type), shape = 15, size = dot_size, show.legend = T) + 
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  coord_fixed() + coord_equal() + 
  ggthemes::theme_map(base_size = 10) +
  theme_minimal() +
  theme(legend.position = "none",
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(color = "gray80", fill = "gray90"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",# panel.background = element_rect(fill = "lightblue1", color = "NA")
        #axis.ticks = element_text(color = "gray9")
        axis.text = element_blank(), #element_text(size = 12),
        axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(title = paste0("Lateral spatial pruned (50%)"), 
       subtitle = "", 
       color = "CANAPE", fill = "CANAPE",
       x = "Longitude", y = "Latitude")

# PE significance
lat_pe_sig_p <- late_canape_prun_df %>% 
  select(-sim_n, -samp_frac, -treat, -imputation) %>% 
  filter(!is.na(lat),
         endem_type != "not significant") %>%
  mutate(pe_obs_sig = case_when(pe_obs_p_upper >= .95 ~ "*", TRUE ~ ""),
         pe_alt_sig = case_when(pe_alt_obs_p_upper >= .95 ~ "*", TRUE ~ "")) %>% 
  gather(key, value, -site, -endem_type, -pe_obs_p_upper, -pe_alt_obs_p_upper,-lat, -long) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "gray80", fill = "gray90") + 
  geom_point(aes(col = endem_type, fill = endem_type), shape = 15, size = dot_size, alpha = .2, show.legend = T) + 
  geom_text(aes(label = value), size = 2.5, color = "black") +
  facet_wrap(~key, nrow = 2) +
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  coord_fixed() + coord_equal() + 
  guides(color = guide_legend(override.aes = list(size = 4.5))) +
  labs(title = "Lateral - Spatial pruned (50%) ", subtitle = "",
       color = "CANAPE", fill = "CANAPE",
       x = "Longitude", y = "Latitude") +
  ggthemes::theme_map(base_size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom", #legend.position = c(0.15, 0.85),
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(color = "gray80", fill = "gray90"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside", # panel.background = element_rect(fill = "lightblue1", color = "NA")
        axis.text = element_text(size = 12),
        axis.title.x     = element_blank(), axis.title.y = element_blank(),
  )

# confusion matrix
lat_spat_cf_mat <- 
left_join(late_canape_prun_df,
          empirical_spat_df %>% select(site, emp_endem_type = endem_type), by = join_by(site)
          ) %>% 
  mutate(endem_type = endem_type %>% factor, emp_endem_type = emp_endem_type %>% factor) %>% 
  #group_by(endem_type) %>% 
  conf_mat(truth = emp_endem_type, estimate = endem_type, na_rm = T)

lat_spat_cf_mat %>% summary

## imputed  consensus #####
conse_canaper <- canape_cons_df %>% 
  filter(!is.na(lat)) %>% 
  ggplot(aes(x = long, y = lat)) +
  geom_polygon(data = japan, aes(x = long, y = lat, group = group), 
               color = "gray80", fill = "gray90") + 
  geom_point(aes(col = endem_consensus, fill = endem_consensus), shape = 15, size = dot_size, show.legend = T) + 
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  coord_fixed() + coord_equal() + 
  ggthemes::theme_map(base_size = 10) +
  theme_minimal() +
  theme(legend.position = "none",
        #legend.frame = element_rect(colour = "lightgray", fill = "lightgray"),
        legend.key = element_rect(color = "gray80", fill = "gray90"),
        legend.key.size = unit(.5, "cm"),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement  = "outside",# panel.background = element_rect(fill = "lightblue1", color = "NA")
        #axis.ticks = element_text(color = "gray9")
        axis.text = element_blank(), #element_text(size = 12),
        axis.title.x     = element_blank(), axis.title.y = element_blank()
  ) + 
  labs(title = paste0("Taxonomic imputed on a lateral spatial pruned (50%)"), 
       subtitle = paste0("Consensus at ", threshold*100, "%"), 
       color = "CANAPE", fill = "CANAPE",
       x = "Longitude", y = "Latitude")

## confusion matrix
# does not work because uneven levels
consens_cf_mat <- 
  left_join(canape_cons_df,
            empirical_spat_df %>% select(site, emp_endem_type = endem_type), 
            by = join_by(site)
            ) %>% 
  mutate(endem_type = endem_consensus %>% factor, emp_endem_type = emp_endem_type %>% factor) %>% 
  #group_by(endem_type) %>% 
  conf_mat(truth = emp_endem_type, estimate = endem_type, na_rm = T) 

# canape histograms ####
cpr_endem_cols_5[['not significant']] <- "gray60"

# full empirical tree
emp_hist <- empirical_spat_df %>% 
  filter(!is.na(endem_type)) %>% 
  group_by(endem_type) %>% tally %>% 
  ggplot(aes(x = "", y = n, fill = endem_type)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  labs(x = "", y = "") +
  scale_fill_manual(values = cpr_endem_cols_5) +
  theme_minimal() + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.text = element_blank(), 
        panel.grid = element_blank())

# pruned spatial tree
prun_hist <- late_canape_prun_df %>% 
  filter(!is.na(endem_type)) %>% 
  group_by(endem_type) %>% tally %>% 
  ggplot(aes(x = "", y = n, fill = endem_type)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  labs(x = "", y = "") +
  scale_fill_manual(values = cpr_endem_cols_5) +
  theme_minimal() + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.text = element_blank(), 
        panel.grid = element_blank())

# multiple imputed trees
randt_hist <- canape_cons_df %>% 
  filter(!is.na(endem_consensus)) %>% 
  group_by(endem_consensus) %>% tally %>% 
  ggplot(aes(x = "", y = n, fill = endem_consensus)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  labs(x = "", y = "") +
  scale_fill_manual(values = cpr_endem_cols_5) +
  theme_minimal() + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.text = element_blank(), 
        panel.grid = element_blank())

# panel map plots ####
# CANAPE
pl_long <- "
  AAACCC
  AAACCC
  BBBDDD
  BBBDDD
  "

emp_canaper5 + rand_prun_canaper + late_prun_canaper + conse_canaper +
  plot_layout(guides = "keep", design = pl_long, tag_level = "keep") + 
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.margin = margin(1, 5, 5, 5))
ggsave(filename = paste0("figures/map_compara_canaper", frt_out), 
       width = 25, height = 25, units = "cm", dpi = 200)

# PE significance
emp_pe_sig_p + lat_pe_sig_p + rand_pe_sig_p +
  plot_layout(guides = "collect", axis_titles = "collect", tag_level = "new")  &
  theme(legend.position = 'bottom')
ggsave(paste0("figures/", out_fold,"/canape_pe_signif_pnl", ".pdf"), width = 30, height = 25, units = "cm", dpi = 300)

# alluvial plots ####
library(ggalluvial)

# set pallet with 5 colors 
cpr_endem_cols_5 <- c(cpr_endem_cols_2, "#A9838F") # "#D996BB" "#CC79A7" "#E5B2CE" "#E6C2D6"
cpr_endem_cols_5[['not significant']] <- "gray70"
names(cpr_endem_cols_5)[6] <- c('no consensus')
# remove super from legend
cpr_endem_cols_5 <- cpr_endem_cols_5[names(cpr_endem_cols_5) != "super"] 
#scales::show_col(cpr_endem_cols_5)
names(cpr_endem_cols_5) <- str_to_title(names(cpr_endem_cols_5))

# load all results
canaper_res_df <- readRDS(paste0("output/canaper/canaper_all_data_res_fern.rds"))

# join with empirical df
canaper_res_all_df <- left_join(canaper_res_df, 
                                emp_df %>% rename(emp_endem_type = endem_type), 
                                by = join_by(site))

# randomly pruned 50%
tab <- canaper_res_all_df %>% 
  filter(imputation == "none", 
         samp_frac == 0.5,
         treat == "rand-pr") %>% 
  filter(!is.na(endem_type)) %>% 
  select(rand_pr = (endem_type), full = emp_endem_type) %>% 
  table

df <- as.data.frame(tab) %>% 
  mutate(
    full = str_to_title(as.character(full)),
    rand_pr = str_to_title(as.character(rand_pr)),
    full = factor(full, levels = c("Neo", "Mixed", "Super", "Paleo", "Not Significant")),
    rand_pr = factor(rand_pr, levels = c("Neo", "Mixed", "Super", "Paleo", "Not Significant"))
  )

random_aluv_pt <- 
  ggplot(df, aes(axis1 = full, axis2 = rand_pr, y = Freq)) +
  geom_alluvium(aes(fill = full), width = 1/4, alpha = 0.6) +
  geom_stratum(aes(fill = after_stat(stratum), color = after_stat(stratum)), width = 1/4) +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Full", "Pruned"), expand = c(0.1, 0.05)) +
  #labs(x = NULL, y = "Count", fill = "late_spat_pr") +
  theme_minimal(base_line_size = 0, base_rect_size = 0) + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title = element_blank(), legend.position = "none") +
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") + 
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  labs(title = paste0("Randomly pruned (50%)"))

# lateral spatially pruned 50%
tab <- canaper_res_all_df %>% 
  filter(imputation == "none", 
         samp_frac == 0.5,
         treat == "late-spat_pr") %>% 
  filter(!is.na(endem_type)) %>% 
  select(late_spat_pr = (endem_type), full = emp_endem_type) %>% 
  table

df <- as.data.frame(tab) %>% 
  mutate(
    full = str_to_title(as.character(full)),
    late_spat_pr = str_to_title(as.character(late_spat_pr)),
    full = factor(full, levels = c("Neo", "Mixed", "Super", "Paleo", "Not Significant")),
    late_spat_pr = factor(late_spat_pr, levels = c("Neo", "Mixed", "Super", "Paleo", "Not Significant"))
  )

late_aluv_pt <- 
  ggplot(df, aes(axis1 = full, axis2 = late_spat_pr, y = Freq)) +
  geom_alluvium(aes(fill = full), width = 1/4, alpha = 0.6) +
  geom_stratum(aes(fill = after_stat(stratum), color = after_stat(stratum)), width = 1/4) +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Full", "Pruned"), expand = c(0.1, 0.05)) +
  #labs(x = NULL, y = "Count", fill = "late_spat_pr") +
  theme_minimal(base_line_size = 0, base_rect_size = 0) + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title = element_blank(), legend.position = "none") +
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") + 
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  labs(title = paste0("Lateral spatial pruned (50%)"))

# randtip
df <- left_join(canape_cons_df %>% select(site, imputed = endem_consensus), 
                emp_df %>% select(site, full = endem_type), 
                by = join_by(site == site)) %>% select(-site) %>% table %>% 
  as.data.frame() %>% 
  mutate(
    full = str_to_title(as.character(full)),
    imputed = str_to_title(as.character(imputed)),
    full = factor(full, levels = c("Neo", "Mixed", "Super", "Paleo", "Not Significant")),
    imputed = factor(imputed, levels = c("Neo", "Mixed", "Super", "Paleo", "No Consensus", "Not Significant"))
  )

imp_aluv_pt <- 
  ggplot(df, aes(axis1 = full, axis2 = imputed, y = Freq)) +
  geom_alluvium(aes(fill = full, 
                    alpha = after_stat(stratum) %in% c("Neo", "Super", "Paleo", "No Consensus")), width = 1/4) +
  geom_stratum(aes(fill = after_stat(stratum), color = after_stat(stratum)), width = 1/4) +
  scale_alpha_manual(values = c(`FALSE` = 0.55, `TRUE` = .75), guide = "none") +
  #geom_text(stat = "stratum", aes(label = after_stat(count)), size = 3) +
  scale_x_discrete(limits = c("Full", "Imputed"), expand = c(0.1, 0.05)) +
  #labs(x = NULL, y = "Count", fill = "late_spat_pr") +
  theme_minimal(base_line_size = 0, base_rect_size = 0) + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title = element_blank(), legend.position = "right") +
  scale_fill_manual( values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") + 
  scale_color_manual(values = cpr_endem_cols_5, limits = names(cpr_endem_cols_5), breaks = names(cpr_endem_cols_5), drop = F, name = "CANAPE") +
  labs(title = paste0("Taxonomic imputed on a lateral spatial pruned (50%)"), 
       subtitle = paste0("Consensus at ", threshold*100, "%"))

random_aluv_pt + late_aluv_pt + imp_aluv_pt + 
  plot_layout(guides = "collect", tag_level = "keep") + 
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.margin = margin(1, 5, 5, 5))
ggsave(filename = paste0("figures/aluvial_compara_canaper", frt_out), 
       width = 30, height = 15, units = "cm", dpi = 200)

# maps + alluvial plots
pl_long <- "
  AAABBBCCC
  AAABBBCCC
  AAABBBCCC
  DDDEEEFFF
  DDDEEEFFF
  "

(rand_prun_canaper) + 
  (late_prun_canaper) + 
  (conse_canaper) +
  (random_aluv_pt + theme(plot.title = element_blank())) +
  (late_aluv_pt + theme(plot.title = element_blank())) + 
  (imp_aluv_pt + theme(plot.title = element_blank(), plot.subtitle = element_blank())) + 
  plot_layout(guides = "collect", design = pl_long, tag_level = "keep") + 
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.margin = margin(1, 5, 5, 5))
ggsave(filename = paste0("figures/map_aluv_canaper", frt_out), 
       width = 35, height = 20, units = "cm", dpi = 200)

