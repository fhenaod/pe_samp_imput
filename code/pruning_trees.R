library(tidyverse)
library(geiger)

source("code/pe_samp_funs.R")

# load data ####
# load community matrix
fern_com <- read.csv("data/nitta-et-al/japan_ferns_comm_full.csv", header = T) %>% 
  column_to_rownames("grids")

# load empirical fern tree
fern_tre <- read.tree("output/nitta_node-named.tre") 
is.binary(fern_tre) ; is.ultrametric(fern_tre)

dir.create("output/pruned_trees/")

# quantitative character ####
# sim or load
n_reps <- 100
run_char <- F
if(run_char){
  sim_char <- sim.char(targ_tre, par = .2, nsim = n_reps, model = "BM")
  sim_char %>% saveRDS("output/sim_char_bm.rds")
} else {
  sim_char <- readRDS("output/sim_char_bm.rds")
  class(sim_char); dim(sim_char)
}

to_prune <- c("fern", "young", "old", "cnt", "exp")
out_fold <- to_prune[1]

# *Fern tree* #####
out_fold <- "fern"
dir.create(paste0("output/pruned_trees/", out_fold))
targ_tre <- fern_tre

# random pruning ####
samp_prop <- c(.05, .15, .3, .5) ; n_reps <- 100

prn_trs <- list()
for(j in 1:length(samp_prop)){
tips2drop <- list()
  for(i in 1:n_reps){
    tips2drop[[i]] <- sample(1:Ntip(targ_tre), 
                             size = floor(Ntip(targ_tre)*samp_prop[j]), 
                             replace = F)
    prn_trs <- lapply(X = tips2drop, function(X) drop.tip(phy = targ_tre, tip = X))
  }
prn_trs %>% 
  saveRDS(paste0("output/pruned_trees/", out_fold, "/rand_pruned_trs_", samp_prop[j], ".rds"))
}

## lateral pruning ####
# classifies continuous values into five percentiles, and samples/prunes the first four in increments
lat_samp_prop <- c(.05, .15, .3, .5) 
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
  lat_trs %>% 
    saveRDS(paste0("output/pruned_trees/", out_fold, "/lateral_pruned_trs_", lat_samp_prop[j], ".rds"))
}

## bilateral pruning ####
marg_samp_prop <- seq(from = .025, to = 0.975, by = .025)
cat2filt <- paste0(rep("perc_", length(marg_samp_prop)+1), c(marg_samp_prop, 1))

# bilateral percentiles to remove
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
    saveRDS(paste0("output/pruned_trees/", out_fold, "/marginal_pruned_trs_", c(length(marg_cat2filt_ls[[i]])*.025), ".rds"))
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

# *Sim trees* #####
# load simulated trees
yng_trees <- readRDS("output/sim_trees_100_age.rds")
old_trees <- readRDS("output/sim_trees_300_age.rds")
cnt_trees <- readRDS("output/sim_trees_217_age.rds")
exp_trees <- readRDS("output/sim_trees_exp_brlen.rds")

trs_ls <- list(yng_trees, old_trees, cnt_trees, exp_trees)

names(trs_ls) <- c("young", "old", "cnt", "exp")

# loop over all simulated trees and pruning proportions
for(l in 1:length(trs_ls)){
  out_fold <- names(trs_ls)[[l]]
  targ_tre <- trs_ls[[l]]
  
  dir.create(paste0("output/pruned_trees/", out_fold))
  
# random pruning ####
samp_prop <- c(.05, .15, .3, .5) ; n_reps <- 100

prn_trs <- list()
for(j in 1:length(samp_prop)){
  tips2drop <- list()
  for(i in 1:n_reps){
    tips2drop[[i]] <- sample(1:Ntip(targ_tre[[i]]), 
                             size = floor(Ntip(targ_tre[[i]])*samp_prop[j]), 
                             replace = F)
    prn_trs <- mapply(FUN = drop.tip, phy = targ_tre, tip = tips2drop, SIMPLIFY = F)
  }
  prn_trs %>% 
    saveRDS(paste0("output/pruned_trees/", out_fold, "/rand_pruned_trs_", samp_prop[j], ".rds"))
}

## lateral pruning ####
# classifies continuous values into five percentiles, and samples/prunes the first four in increments
lat_samp_prop <- c(.05, .15, .3, .5) 
lat_cat2filt <- c("below.05", ".05to.15", ".15to.3", ".3to.5", "above.5")

lat_tips2drop <- list()
for(j in 1:length(lat_samp_prop)){
  lat_tips2drop <- data.frame(sim_char) %>% 
    rownames_to_column(var = "tip_names") %>% 
    mutate(across(where(is.numeric), lat_funct)) %>%
    gather(key = sim_n, value = cat2perc, -tip_names) %>% 
    filter(cat2perc %in% lat_cat2filt[1:j]) %>% 
    group_split(sim_n, .keep = F) %>% lapply(., '[[', "tip_names")
  
  lat_trs <- mapply(FUN = drop.tip, phy = targ_tre, tip = lat_tips2drop, SIMPLIFY = F)
  lat_trs %>% 
    saveRDS(paste0("output/pruned_trees/", out_fold, "/lateral_pruned_trs_", lat_samp_prop[j], ".rds"))
}

## bilateral pruning ####
marg_samp_prop <- seq(from = .025, to = 0.975, by = .025)
cat2filt <- paste0(rep("perc_", length(marg_samp_prop)+1), c(marg_samp_prop, 1))

# bilateral percentiles to remove
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
  
  marg_trs <- mapply(FUN = drop.tip, phy = targ_tre, tip = marg_tips2drop, SIMPLIFY = F)
  marg_trs %>% 
    saveRDS(paste0("output/pruned_trees/", out_fold, "/marginal_pruned_trs_", c(length(marg_cat2filt_ls[[i]])*.025), ".rds"))
}
}

# check  pruned tree sizes ####
for(j in 1:length(trs_ls)){
  print(toupper(names(trs_ls)[[j]]))
  
  print("Random pruning")
  rand_fls <- list.files(paste0("output/pruned_trees/", out_fold), pattern = "rand_pruned", all.files = T, full.names = T)
for(i in 1:length(rand_fls)){
  targ_tr_ls <- readRDS(rand_fls[i])
  sapply(targ_tr_ls, Ntip) %>% 
    table() %>% print
}

  print("Lateral pruning")
late_fls <- list.files(paste0("output/pruned_trees/", out_fold), pattern = "lateral_pruned", all.files = T, full.names = T)
for(i in 1:length(late_fls)){
  targ_tr_ls <- readRDS(late_fls[i])
  sapply(targ_tr_ls, Ntip) %>% 
    #max %>% 
    table() %>% print
}

  print("Marginal pruning")
marg_fls <- list.files(paste0("output/pruned_trees/", out_fold), pattern = "marginal_pruned", all.files = T, full.names = T)
for(i in 1:length(marg_fls)){
  targ_tr_ls <- readRDS(marg_fls[i])
  sapply(targ_tr_ls, Ntip) %>% 
    #max %>% 
    table() %>% print
}
}
# viz tree with sim BM characters ####
library(ggtree)

# plot five random simulated characters
data.frame(sim_char) %>% 
  select(sample(1:100, 5)) %>% 
  #rownames_to_column(var = "tip_names") %>%
  gheatmap(ggtree(targ_tre), ., offset = 8, width = .15, 
         colnames = F, legend_title = "trait\nvalue") +
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

