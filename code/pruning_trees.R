library(tidyverse)
library(geiger)

#load tree
fern_tre <- read.tree("output/nitta_node-named.tre") 
is.binary(fern_tre) ; is.ultrametric(fern_tre)
targ_tre <- fern_tre

dir.create("output/pruned_trees/")

# random pruning ####
samp_prop <- c(.05, .15,.3, .5) ; n_reps <- 100

prn_trs <- list()
for(j in 1:length(samp_prop)){
tips2drop <- list()
  for(i in 1:n_reps){
    tips2drop[[i]] <- sample(1:Ntip(targ_tre), 
                             size = floor(Ntip(targ_tre)*samp_prop[j]), 
                             replace = F)
    prn_trs <- lapply(X = tips2drop, function(X) drop.tip(phy = targ_tre, tip = X))
  }
prn_trs %>% saveRDS(paste0("output/pruned_trees/rand_pruned_trs_", samp_prop[j], ".rds"))
}

# sim character ####
n_reps <- 100
run_char <- F
if(run_char){
  sim_char <- sim.char(targ_tre, par = .2, nsim = n_reps, model = "BM")
  sim_char %>% saveRDS("output/sim_char_bm.rds")
} else {
  sim_char <- readRDS("output/sim_char_bm.rds")
  class(sim_char); dim(sim_char)
}

## lateral pruning ####
# classifies continuous values into five percentiles, and samples/prunes the first four in increments
lat_samp_prop <- c(.05, .15, .3, .5) 
lat_funct=function(x){cut(x, 
                          breaks = c(-Inf, quantile(x, lat_samp_prop), +Inf), 
                          labels = c("below.05", ".05to.15", ".15to.3", ".3to.5", "above.5"))
  }

lat_cat2filt <- c("below.05", ".05to.15", ".15to.3", ".3to.5", "above.5")

lat_tips2drop <- list()
for(j in 1:length(lat_samp_prop)){
  lat_tips2drop <- data.frame(sim_char) %>% 
    rownames_to_column(var = "tip_names") %>% 
    mutate(across(where(is.numeric), lat_funct)) %>%
    gather(key = sim_n, value = cat2perc, -tip_names) %>% 
    filter(cat2perc %in% lat_cat2filt[1:j]) %>% 
    group_split(sim_n, .keep = F) %>% lapply(., '[[', "tip_names")
  
  lat_trs <- lapply(X = lat_tips2drop, function(X) drop.tip(phy = targ_tre, tip = X))
  lat_trs %>% saveRDS(paste0("output/pruned_trees/lateral_pruned_trs_", lat_samp_prop[j], ".rds"))
}

## marginal pruning ####
marg_samp_prop <- seq(from = .025, to = 0.975, by = .025)
cat2filt <- paste0(rep("perc_", length(marg_samp_prop)+1), c(marg_samp_prop, 1))

# classifies continuous values into percentiles
margin_funct=function(x){cut(x, 
                             breaks = c(-Inf, quantile(x, marg_samp_prop), +Inf),
                             labels = cat2filt)
  }

# marginal percentiles to remove
marg_cat2filt_ls <- list(c(head(cat2filt, 1), tail(cat2filt, 1)),
                    c(head(cat2filt, 2), tail(cat2filt, 2)),
                    c(head(cat2filt, 4), tail(cat2filt, 4)),
                    c(head(cat2filt, 6), tail(cat2filt, 6)),
                    c(head(cat2filt, 10), tail(cat2filt, 10))
                    )

marg_trs <- list()
for(i in 1:length(marg_cat2filt_ls)){
  marg_tips2drop <- data.frame(sim_char) %>% rownames_to_column(var = "tip_names") %>% 
    mutate(across(where(is.numeric), margin_funct)) %>% 
    gather(key = sim_n, value = cat2perc, -tip_names) %>% #head()
    filter(cat2perc %in% marg_cat2filt_ls[[i]]) %>% 
    group_split(sim_n, .keep = F) %>% lapply(., '[[', "tip_names")
  
  marg_trs <- lapply(X = marg_tips2drop, function(X) drop.tip(phy = targ_tre, tip = X))
  marg_trs %>% 
    saveRDS(paste0("output/pruned_trees/marginal_pruned_trs_", c(length(marg_cat2filt_ls[[i]])*.025), ".rds"))
}

# check  pruned tree sizes ####
rand_fls <- list.files("output/pruned_trees", "rand_pruned", all.files = T, full.names = T)
for(i in 1:length(rand_fls)){
  targ_tr_ls <- readRDS(rand_fls[i])
  sapply(targ_tr_ls, Ntip) %>% table() %>% print
}

late_fls <- list.files("output/pruned_trees", "lateral_pruned", all.files = T, full.names = T)
for(i in 1:length(late_fls)){
  targ_tr_ls <- readRDS(late_fls[i])
  sapply(targ_tr_ls, Ntip) %>% 
    max %>% 
    #table() %>% 
    print
}

marg_fls <- list.files("output/pruned_trees", "marginal_pruned", all.files = T, full.names = T)
for(i in 1:length(marg_fls)){
  targ_tr_ls <- readRDS(marg_fls[i])
  sapply(targ_tr_ls, Ntip) %>% 
    max %>% 
    #table() %>% 
    print
}

# viz BM sims  ####
library(ggtree)
ggtree(targ_tre, aes(color=trait), continuous = 'colour', yscale = "trait") + 
  scale_color_viridis_c() + theme_minimal()

# plot five random simulated characters
data.frame(sim_char) %>% 
  select(sample(1:100, 5)) %>% 
  #rownames_to_column(var = "tip_names") %>%
  gheatmap(ggtree(targ_tre), ., offset = 8, width = .15, 
         colnames = F, legend_title = "trait\nvalue", ) +
  scale_x_ggtree() + 
  scale_y_continuous(expand = c(0, 0.2)) + 
  theme(legend.position = "none")
ggsave("figures/tree_simchar_s5.png", width = 25, height = 25, units = "cm", dpi = 200)

# plot all random simulated characters
data.frame(sim_char) %>% 
  gheatmap(ggtree(targ_tre), ., offset = 8, width = 0.6, 
           colnames = F, legend_title = "trait\nvalue", ) +
  scale_x_ggtree() + 
  scale_y_continuous(expand = c(0, 0.2)) + 
  theme(legend.position = "none")
ggsave("figures/tree_simchar.png", width = 25, height = 25, units = "cm", dpi = 200)

# viz pruning data in hist ####
library(ggpubr)
# lateral
v2samp <- sample(2:100, 3, replace = F)

data2plt <- data.frame(sim_char) %>% rownames_to_column(var = "tip_names") %>% 
  select(1, all_of(v2samp)) %>% 
  gather(key = sim_n, val, -tip_names) %>% 
  group_by(sim_n) %>% 
  mutate(perc = lat_funct(val)) 

lat2print <- list()
for(i in 1:c(length(lat_cat2filt)-1)){
  lat2print[[i]] <- data2plt %>% 
    mutate(is_pr = ifelse(perc %in% lat_cat2filt[1:i], "pruned", "preserved")) %>%
    ggplot(aes(x = val, color = is_pr, fill = is_pr)) +
    geom_histogram() +
    theme_classic() + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    labs(subtitle = paste0("Lateral pruning (", lat_samp_prop[i]*100, "%)"),
         x = "Trait value", y = "Count") +
    scale_color_grey(name = "State", start = .55, end = .15) + 
    scale_fill_grey( name = "State", start = .75, end = .15) +
    geom_vline(data = data2plt %>% plyr::ddply("sim_n", summarise, grp.mean = mean(val)), 
               aes(xintercept = grp.mean), color = "coral",
               linetype = "dashed")+
    facet_wrap(~sim_n, scales = "free")
  
}
ggarrange(plotlist = lat2print, ncol = 1, common.legend = T)
ggsave("figures/lateral_prun.png", 
       width = 20, height = 25, units = "cm", dpi = 200, bg = "white")
 
# marginal
data2plt <- data.frame(sim_char) %>% 
  rownames_to_column(var = "tip_names") %>% 
  select(1, all_of(v2samp)) %>% 
  gather(key = sim_n, val, -tip_names) %>% 
  group_by(sim_n) %>% 
  mutate(perc = margin_funct(val)) 

marg2print <- list()
for(i in 1:length(marg_cat2filt_ls)){
marg2print[[i]] <- data2plt %>% 
    mutate(is_pr = ifelse(perc %in% marg_cat2filt_ls[[i]], "pruned", "preserved")) %>%
    ggplot(aes(x = val, color = is_pr, fill = is_pr)) +
    geom_histogram() +
    theme_classic() + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    labs(subtitle = paste0("Marginal pruning (", c(length(marg_cat2filt_ls[[i]])*.025)*100, "%)"),
         x = "Trait value", y = "Count") +
    scale_color_grey(name = "State", start = .55, end = .15) + 
    scale_fill_grey( name = "State", start = .75, end = .15) +
    geom_vline(data = data2plt %>% plyr::ddply("sim_n", summarise, grp.mean = mean(val)), 
               aes(xintercept = grp.mean), color = "coral",
               linetype = "dashed") +
    facet_wrap(~sim_n, scales = "free", ncol = 3)
}
ggarrange(plotlist = marg2print, ncol = 1, common.legend = T)
ggsave("figures/marginal_prun.png", 
       width = 20, height = 25, units = "cm", dpi = 200, bg = "white")

### with densities
data2plt %>% 
  ggplot(aes(x = val, color = perc, fill = perc)) +
  #geom_histogram(fill  = "red",
  #               alpha = .5,
  #               position = "identity") +
  geom_histogram(aes(y = ..density..), col = "darkgray", position = "identity",  alpha = .5) +
  geom_density(alpha = .85) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(title = "Lateral pruning", subtitle = "",
       x = "Trait value", y = "Density") +
  scale_color_grey(name = "Percentile") + scale_fill_grey(name = "Percentile") +
  geom_vline(data = data2plt %>% plyr::ddply("sim_n", summarise, grp.mean = mean(val)), 
             aes(xintercept = grp.mean), color = "coral",
             linetype = "dashed") +
  facet_wrap(~sim_n, scales = "free")
