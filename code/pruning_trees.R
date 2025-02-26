library(tidyverse)
library(geiger)

#load tree
cali_tre <- read.nexus("data/California_clades_fully_dated.nex")

# random pruning ####
samp_prop <- c(.05, .15,.3, .5) ; n_reps <- 100

prn_trs <- list()
for(j in 1:length(samp_prop)){
tips2drop <- list()
  for(i in 1:n_reps){
    tips2drop[[i]] <- sample(1:Ntip(cali_tre), 
                             size = floor(Ntip(cali_tre)*samp_prop[j]), 
                             replace = F)
    prn_trs <- lapply(X = tips2drop, function(X) drop.tip(phy = cali_tre, tip = X))
  } 
prn_trs %>% saveRDS(paste0("output/rand_pruned_trs_", samp_prop[j], ".rds"))
}

# sim character ####
n_reps <- 100
sim_char <- sim.char(cali_tre, par = .2, nsim = n_reps, model = "BM")
sim_char %>% saveRDS("output/sim_char_bm.rds")

## lateral pruning ####
# classifies continuous values into percentiles
lat_samp_prop <- c(.05, .15,.3, .5) 
lat_funct=function(x){cut(x, 
                          breaks = c(-Inf, quantile(x, lat_samp_prop), +Inf), 
                          labels = c("below.05", ".05to.15", ".15to.3", ".3to.5", "above.5"))
  }

lat_cat2filt <- c("below.05", ".05to.15", ".15to.3", ".3to.5", "above.5")

lat_tips2drop <- list()
for(j in 1:c(length(lat_cat2filt)-1)){
  lat_tips2drop <- data.frame(sim_char) %>% rownames_to_column(var = "tip_names") %>% 
    mutate(across(where(is.numeric), lat_funct)) %>% 
    gather(key = sim_n, value = cat2perc, -tip_names) %>% 
    filter(cat2perc == lat_cat2filt[j]) %>% 
    group_split(sim_n, .keep = F) %>% lapply(., '[[', "tip_names")
  
  lat_trs <- lapply(X = lat_tips2drop, function(X) drop.tip(phy = cali_tre, tip = X))
  lat_trs %>% saveRDS(paste0("output/lateral_pruned_trs_", lat_samp_prop[j], ".rds"))
}

## marginal pruning ####
marg_samp_prop <- seq(from = .025, to = 0.975, by = .025)
cat2filt <- paste0(rep("perc_", length(marg_samp_prop)+1), c(marg_samp_prop, 1))

# classifies continuous values into percentiles
margin_funct=function(x){cut(x, 
                             breaks = c(-Inf, quantile(x, marg_samp_prop), +Inf),
                             labels = cat2filt
                          )
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
  
  marg_trs <- lapply(X = marg_tips2drop, function(X) drop.tip(phy = cali_tre, tip = X))
  marg_trs %>% 
    saveRDS(paste0("output/marginal_pruned_trs_", c(length(marg_cat2filt_ls[[i]])*.025), ".rds"))
}

