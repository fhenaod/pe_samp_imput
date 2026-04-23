library(tidyverse)
library(ggtree)
library(patchwork)

source("code/pe_samp_funs.R")

# load data ####
# load community matrix
fern_com <- read.csv("data/nitta-et-al/japan_ferns_comm_full.csv", header = T) %>% 
  column_to_rownames("grids")

# spatial data
spp_dist_df <- fern_com %>% colSums() %>% data.frame(n_grid = .) %>% 
  rownames_to_column("tip_names") %>% 
  mutate(prop_grid = n_grid/dim(fern_com)[1],
         ln_prop_grd = log(prop_grid))

# load tree
fern_tre <- read.tree("output/nitta_node-named.tre") 
ape::is.binary(fern_tre) ; ape::is.ultrametric(fern_tre)
targ_tre <- fern_tre

sim_char <- readRDS("output/sim_char_bm.rds")

# set tree type ####
typ_tre <-  "rectangular" # "fan" #"roundrect" #

# full tree
fern_full <- fern_tre %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) + 
  geom_nodepoint(size = .15, alpha = .55) +
  labs(title = "Full tree")

## tree with geological timeline ####
library(deeptime)
fern_tre %>% 
  ggtree(layout = "fan", ladderize = T, size = .2) + 
  gggeo_scale(
    time_scale = "periods", # Or use "eras", "epochs", etc.
    pos = "bottom", # Position of the timescale
    offset = 0.5  # Space between tree and timescale
  )

revts(ggtree(fern_tre, layout = "rectangular", ladderize = T, size = .2))  + 
  coord_geo(xlim = c(-440, 0), ylim = c(-2, Ntip(fern_tre)), 
            pos = "bottom", size = 5, fontface = "bold", #alpha = 0.95,
            #dat = list("stages", "epochs", "periods"),
            neg = T, abbrv = T) + 
  scale_x_continuous(breaks = seq(-440, 0, 50), labels = abs(seq(-440, 0, 50))) +
  theme_tree2() +
  theme(axis.text.x = element_text(size = 12, face = "bold")) 
ggsave(filename = paste0("figures/fern_tree_time.pdf"), 
       width = 15, height = 27, units = "cm", dpi = 200)

# Pruned phylogenetic trees ####
## spatial ####
late_spat_fls <- dir("output/pruned_trees/", pattern = "lateral_spat_pruned", full.names = T)
bila_spat_fls <- dir("output/pruned_trees/", pattern = "marginal_spat_pruned", full.names = T)

# lateral
# pruned trees
late_spat_05_pr_tr <- readRDS(late_spat_fls[[1]]) 
late_spat_15_pr_tr <- readRDS(late_spat_fls[[2]])
late_spat_30_pr_tr <- readRDS(late_spat_fls[[3]])
late_spat_50_pr_tr <- readRDS(late_spat_fls[[4]])

# 5%
late_spat_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, late_spat_05_pr_tr$tip.label), 
                         keept  = intersect(fern_tre$tip.label, late_spat_05_pr_tr$tip.label))
lat_spat_pr05_p <- groupOTU(fern_full, late_spat_prn_ls, 'tips') +
  aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Lateral pruning\n(5%)")

# 15%
late_spat_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, late_spat_15_pr_tr$tip.label), 
                         keept  = intersect(fern_tre$tip.label, late_spat_15_pr_tr$tip.label))
lat_spat_pr15_p <- groupOTU(fern_full, late_spat_prn_ls, 'tips') +
  aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Lateral pruning\n(15%)")

# 30% 
late_spat_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, late_spat_30_pr_tr$tip.label), 
                         keept  = intersect(fern_tre$tip.label, late_spat_30_pr_tr$tip.label))
lat_spat_pr30_p <- groupOTU(fern_full, late_spat_prn_ls, 'tips') +
  aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Lateral pruning\n(30%)")

# 50% 
late_spat_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, late_spat_50_pr_tr$tip.label), 
                         keept  = intersect(fern_tre$tip.label, late_spat_50_pr_tr$tip.label))
lat_spat_pr50_p <- groupOTU(fern_full, late_spat_prn_ls, 'tips') +
  aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Lateral pruning\n(50%)")

# bilateral
# pruned trees
bil_spat_10_pr_tr <- readRDS(bila_spat_fls[[1]]) 
bil_spat_20_pr_tr <- readRDS(bila_spat_fls[[2]])
bil_spat_30_pr_tr <- readRDS(bila_spat_fls[[3]])
bil_spat_50_pr_tr <- readRDS(bila_spat_fls[[4]])

# 10%
bila_spat_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, bil_spat_10_pr_tr$tip.label), 
                         keept  = intersect(fern_tre$tip.label, bil_spat_10_pr_tr$tip.label))
bila_spat_pr10_p <- groupOTU(fern_full, bila_spat_prn_ls, 'tips') +
  aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Bilateral pruning\n(10%)")

# 20%
bila_spat_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, bil_spat_20_pr_tr$tip.label), 
                         keept  = intersect(fern_tre$tip.label, bil_spat_20_pr_tr$tip.label))
bila_spat_pr20_p <- groupOTU(fern_full, bila_spat_prn_ls, 'tips') +
  aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Bilateral pruning\n(20%)")

# 30% 
bila_spat_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, bil_spat_30_pr_tr$tip.label), 
                         keept  = intersect(fern_tre$tip.label, bil_spat_30_pr_tr$tip.label))
bila_spat_pr30_p <- groupOTU(fern_full, bila_spat_prn_ls, 'tips') +
  aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Lateral pruning\n(30%)")

# 50% 
bila_spat_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, bil_spat_50_pr_tr$tip.label), 
                         keept  = intersect(fern_tre$tip.label, bil_spat_50_pr_tr$tip.label))
bila_spat_pr50_p <- groupOTU(fern_full, bila_spat_prn_ls, 'tips') +
  aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Bilateral pruning\n(50%)")

## quantitative trait #####
rand_fls <- list.files("output/pruned_trees", "rand_pruned",     all.files = T, full.names = T)
late_fls <- list.files("output/pruned_trees", "lateral_pruned",  all.files = T, full.names = T)
marg_fls <- list.files("output/pruned_trees", "marginal_pruned", all.files = T, full.names = T)

rand_frc <- readr::parse_number(rand_fls)
late_frc <- readr::parse_number(late_fls)
marg_frc <- readr::parse_number(marg_fls)

# pruned trees
rand_prun_tree <- readRDS(rand_fls[[3]])[[10]] 
late_prun_tree <- readRDS(late_fls[[3]])[[10]]
marg_prun_tree <- readRDS(marg_fls[[4]])[[10]]

# taxonomic at node imputation
rand_tax_imp_tree <- readRDS(paste0("output/imputed_trees/rand_pr_tax_imp_root_",  rand_frc[3], ".rds"))[[10]]
late_tax_imp_tree <- readRDS(paste0("output/imputed_trees/late_pr_tax_imp_root_",  late_frc[3], ".rds"))[[10]]
marg_tax_imp_tree <- readRDS(paste0("output/imputed_trees/marg_pr_tax_imp_root_",  marg_frc[4], ".rds"))[[10]]

# taxonomic at random within
rand_tax_imp_rand_tree <- readRDS(paste0("output/imputed_trees/rand_pr_tax_imp_random_",  rand_frc[3], ".rds"))[[10]]
late_tax_imp_rand_tree <- readRDS(paste0("output/imputed_trees/late_pr_tax_imp_random_",  late_frc[3], ".rds"))[[10]]
marg_tax_imp_rand_tree <- readRDS(paste0("output/imputed_trees/marg_pr_tax_imp_random_",  marg_frc[4], ".rds"))[[10]]

### Randomly pruned ####
rand_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, rand_prun_tree$tip.label), 
                    keept = intersect(fern_tre$tip.label, rand_prun_tree$tip.label))

t1.1 <- fern_tre %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .1, alpha = .55)  +
  labs(title = "Random pruning (30%)") #+ geom_tiplab(aes(subset = label %in% grp1$pruned), color = 'firebrick', size = 1.2)

ran_prn <- groupOTU(t1.1, rand_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#f3a322", "black"))

### Randomly pruned and taxonomically imputed (at node) ####
rand_prn_ls <- list(imputed = setdiff(fern_tre$tip.label, rand_prun_tree$tip.label), 
                    keept = intersect(fern_tre$tip.label, rand_prun_tree$tip.label))
t2.1 <- rand_tax_imp_tree %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .1, alpha = .55)  +
  labs(title = "Random pruning (30%)", subtitle = "Taxonomic imputed (node)") #+ geom_tiplab(aes(subset = label %in% grp2$imputed), color = 'firebrick', size = 1.2) #+ scale_x_reverse()

ran_prn_nd_imp <- groupOTU(t2.1, rand_prn_ls, 'tips') + aes(color = tips) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = c("red", "black")) 

### Randomly pruned and taxonomically imputed (randomly within) ####
t3.1 <- rand_tax_imp_rand_tree %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .1, alpha = .55)  +
  labs(title = "Random pruning (30%)", subtitle = "Taxonomic imputed (subclade)") #+ geom_tiplab(aes(subset = label %in% grp4$imputed), color = 'firebrick', size = 1.2)

ran_prn_bin <- groupOTU(t3.1, rand_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "black"))

if(typ_tre == "fan"){
  pl_lay <- "
  AABB
  AABB
  CCDD
  CCDD
  "
} else {
  pl_lay <- "
  AABBCCDD
  AABBCCDD
  AABBCCDD
  "
}

## random pruning and imputation panel ####
free(fern_full) + free(ran_prn) + free(ran_prn_nd_imp) + free(ran_prn_bin) + 
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/comp_trees_rand_prn_imp.pdf"), 
       width = 28, height = 15, units = "cm", dpi = 200)

### lateral pruned #### 
late_pr_.05 <- readRDS(late_fls[[1]])[[10]]
late_pr_.15 <- readRDS(late_fls[[2]])[[10]]
late_pr_.3 <- readRDS(late_fls[[3]])[[10]]
late_pr_.5 <- readRDS(late_fls[[4]])[[10]]

# lateral pruned
late_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, late_pr_.05$tip.label), 
                    keept = intersect(fern_tre$tip.label, late_pr_.05$tip.label))
lat_pr05_p <- groupOTU(fern_full, late_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Lateral pruning\n(5%)")

late_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, late_pr_.15$tip.label), 
                    keept = intersect(fern_tre$tip.label, late_pr_.15$tip.label))
lat_pr15_p <- groupOTU(fern_full, late_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Lateral pruning\n(15%)")

late_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, late_pr_.3$tip.label), 
                    keept = intersect(fern_tre$tip.label, late_pr_.3$tip.label))

lat_pr30_p <- groupOTU(fern_full, late_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Lateral pruning\n(30%)")

late_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, late_pr_.5$tip.label), 
                    keept = intersect(fern_tre$tip.label, late_pr_.5$tip.label))
lat_pr50_p <- groupOTU(fern_full, late_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Lateral pruning\n(50%)")

pl_lay <- "
  AA#BBCCDDEE
  AA#BBCCDDEE
  AA#BBCCDDEE
  "
free(fern_full) + free(lat_pr05_p) + free(lat_pr15_p) + free(lat_pr30_p) + free(lat_pr50_p) +
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/comp_trees_late_prn.pdf"), 
       width = 40, height = 13, units = "cm", dpi = 300)

### Lateraly pruned and taxonomically imputed (at node) ####
late_prn_ls <- list(imputed = setdiff(fern_tre$tip.label, late_prun_tree$tip.label), 
                    keept = intersect(fern_tre$tip.label, late_prun_tree$tip.label))

t2.2 <- late_tax_imp_tree %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .1, alpha = .55)  +
  labs(#title = "Lateral pruning (30%)", 
    title = "Taxonomic imputed (node)")

late_prn_nd_imp <- groupOTU(t2.2, late_prn_ls, 'tips') + aes(color = tips) + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "black"))

### Lateraly pruned and taxonomically imputed (randomly within) ####
t3.2 <- late_tax_imp_rand_tree %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) +
  geom_nodepoint(size = .1, alpha = .55)  +
  labs(#title = "Lateral pruning (30%)", 
    title = "Taxonomic imputed (subclade)") #+ geom_tiplab(aes(subset = label %in% grp4$imputed), color = 'firebrick', size = 1.2)

late_prn_bin <- groupOTU(t3.2, late_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "black"))

## lateral pruning and imputation panel ####
free(lat_pr30_p) + free(late_prn_nd_imp) + free(late_prn_bin) + 
  plot_layout(guides = "collect", design = "AA#BBCC
                                            AA#BBCC", axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/comp_trees_late_prn_imp.png"), 
       width = 28, height = 15, units = "cm", dpi = 300)

####  marginal pruned trees #####
marg_pr_.05 <- readRDS(marg_fls[[1]])[[10]]
marg_pr_.2 <- readRDS(marg_fls[[3]])[[10]]
marg_pr_.3 <- readRDS(marg_fls[[4]])[[10]]
marg_pr_.5 <- readRDS(marg_fls[[5]])[[10]]

marg_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, marg_pr_.05$tip.label), 
                    keept = intersect(fern_tre$tip.label, marg_pr_.05$tip.label))
marg_pr05_p <- groupOTU(fern_full, marg_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Bilateral pruning\n(5%)")

marg_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, marg_pr_.2$tip.label), 
                    keept = intersect(fern_tre$tip.label, marg_pr_.2$tip.label))
marg_pr20_p <- groupOTU(fern_full, marg_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Bilateral pruning\n(20%)")

marg_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, marg_pr_.3$tip.label), 
                    keept = intersect(fern_tre$tip.label, marg_pr_.3$tip.label))
marg_pr30_p <- groupOTU(fern_full, marg_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Bilateral pruning\n(30%)")

marg_prn_ls <- list(pruned = setdiff(fern_tre$tip.label, marg_pr_.5$tip.label), 
                    keept = intersect(fern_tre$tip.label, marg_pr_.5$tip.label))
marg_pr50_p <- groupOTU(fern_full, marg_prn_ls, 'tips') + aes(color = tips) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "#f3a322")) +
  labs(title = "Bilateral pruning\n(50%)")

## marginal pruning panel ####
pl_lay <- "
  AA#BBCCDDEE
  AA#BBCCDDEE
  AA#BBCCDDEE
  "
free(fern_full) + free(marg_pr05_p) + free(marg_pr20_p) + free(marg_pr30_p) + free(marg_pr50_p) +
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/comp_trees_bilat_prn.pdf"), 
       width = 40, height = 13, units = "cm", dpi = 300)

# Histograms #####
## geographic range ####
# lateral
lat_samp_prop <- c(.05, .15, .3, .5) 
lat_cat2filt <- c("below.05", ".05to.15", ".15to.3", ".3to.5", "above.5")

data2plt <- spp_dist_df %>% 
  select(tip_names, ln_prop_grd) %>% 
  mutate(val = ln_prop_grd, perc = lat_funct(ln_prop_grd))

lat2print_spat <- list()
for(i in 1:c(length(lat_cat2filt)-1)){
  lat2print_spat[[i]] <-  data2plt%>%
    mutate(is_pr = ifelse(perc %in% lat_cat2filt[1:i], "pruned", "preserved")) %>%
    ggplot(aes(x = val, color = is_pr, fill = is_pr)) +
    geom_histogram() +
    theme_classic() + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()
          , legend.position = "none") +
    labs(#subtitle = paste0("Lateral pruning\n(", lat_samp_prop[i]*100, "%)"),
      x = "", y = ""
      #x = "Log Geographic range", y = "Count"
    ) +
    scale_fill_manual(values = c("darkgray", "#f3a322")) +
    scale_color_manual(values = c("black", "black")) +
    #scale_color_grey(name = "State", start = .55, end = .15) + scale_fill_grey( name = "State", start = .75, end = .15) +
    geom_vline(aes(xintercept = mean(val)), color = "coral", linetype = "dashed")
}
ggpubr::ggarrange(plotlist = lat2print_spat, nrow = 1, common.legend = T)

full_spat_hist <- data2plt %>% 
  ggplot(aes(x = val)) +
  geom_histogram(color = "black", fill = "darkgray") +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") +
  labs(#subtitle = paste0("Full tree"),
    #x = "", y = "",
    x = "Log Geographic range", y = "Count") +
  geom_vline(aes(xintercept = mean(val)), color = "coral", linetype = "dashed")

## BM simulated character ####
#### lateral ####
v2samp <- sample(2:100, 1, replace = F)

data2plt <- data.frame(sim_char) %>% 
  rownames_to_column(var = "tip_names") %>% 
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
          strip.text.x = element_blank()
          , legend.position = "none") +
    labs(#subtitle = paste0("Lateral pruning\n(", lat_samp_prop[i]*100, "%)"),
      x = "", y = ""
      #x = "Trait value", y = "Count"
    ) +
    scale_fill_manual(values = c("darkgray", "#f3a322")) +
    scale_color_manual(values = c("black", "black")) +
    #scale_color_grey(name = "State", start = .55, end = .15) + scale_fill_grey( name = "State", start = .75, end = .15) +
    geom_vline(data = data2plt %>% plyr::ddply("sim_n", summarise, grp.mean = mean(val)), 
               aes(xintercept = grp.mean), color = "coral",
               linetype = "dashed") +
    facet_wrap(~sim_n, scales = "free")
  
}
ggpubr::ggarrange(plotlist = lat2print, ncol = 1, common.legend = T)
ggsave("figures/lateral_prun.png", 
       width = 20, height = 25, units = "cm", dpi = 200, bg = "white")

full_tree <- data2plt %>% 
  ggplot(aes(x = val)) +
  geom_histogram(color = "black", fill = "darkgray") +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") +
  labs(#subtitle = paste0("Full tree"),
    x = "", y = ""
    #x = "Trait value", y = "Count"
  ) +
  geom_vline(data = data2plt %>% plyr::ddply("sim_n", summarise, grp.mean = mean(val)), 
             aes(xintercept = grp.mean), color = "coral",
             linetype = "dashed")

### all hists lateral panel ####
free(full_tree) + free(lat2print[[1]]) + free(lat2print[[2]]) + free(lat2print[[3]]) + free(lat2print[[4]]) +
  plot_layout(guides = "collect", design = "A#BCDE", axis_titles = "collect", tag_level = "keep") 
ggsave(filename = paste0("figures/comp_late_dist_prn.pdf"), 
       width = 40, height = 13, units = "cm", dpi = 300)

#### marginal  ####
marg_samp_prop <- seq(from = .025, to = 0.975, by = .025)
cat2filt <- paste0(rep("perc_", length(marg_samp_prop)+1), c(marg_samp_prop, 1))
marg_cat2filt_ls <- list(c(head(cat2filt, 1), tail(cat2filt, 1)),
                         c(head(cat2filt, 2), tail(cat2filt, 2)),
                         c(head(cat2filt, 4), tail(cat2filt, 4)),
                         c(head(cat2filt, 6), tail(cat2filt, 6)),
                         c(head(cat2filt, 10), tail(cat2filt, 10))
)

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
          strip.text.x = element_blank(), 
          legend.position = "none"
    ) +
    labs(#subtitle = paste0("Marginal pruning (", c(length(marg_cat2filt_ls[[i]])*.025)*100, "%)"),
      x = "", y = ""
      #x = "Trait value", y = "Count"
    ) +
    scale_fill_manual(values = c("darkgray", "#f3a322")) + 
    scale_color_manual(values = c("black", "black")) +
    #scale_color_grey(name = "State", start = .55, end = .15) + scale_fill_grey( name = "State", start = .75, end = .15) +
    geom_vline(data = data2plt %>% plyr::ddply("sim_n", summarise, grp.mean = mean(val)), 
               aes(xintercept = grp.mean), color = "coral",
               linetype = "dashed") +
    facet_wrap(~sim_n, scales = "free", ncol = 3)
}
ggpubr::ggarrange(plotlist = marg2print, ncol = 1, common.legend = T)
ggsave("figures/marginal_prun.png", 
       width = 20, height = 25, units = "cm", dpi = 200, bg = "white")

### all hist marginal panel ####
free(full_tree) + free(marg2print[[2]]) + free(marg2print[[3]]) + free(marg2print[[4]]) + free(marg2print[[5]]) +
  plot_layout(guides = "collect", design = "A#BCDE", axis_titles = "collect", tag_level = "keep") 
ggsave(filename = paste0("figures/comp_late_dist_prn.pdf"), 
       width = 40, height = 13, units = "cm", dpi = 300)

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

# pruned trees and histograms ####
pl_lay <- "
  AA#BBCCDDEE
  AA#BBCCDDEE
  AA#BBCCDDEE
  FF#GGHHIIOO
  "
## spatial ####
# lateral
free(fern_full) + free(lat_spat_pr05_p) + free(lat_spat_pr15_p) + free(lat_spat_pr30_p) + free(lat_spat_pr50_p) +
  free(full_spat_hist) + free(lat2print_spat[[1]]) + free(lat2print_spat[[2]]) + free(lat2print_spat[[3]]) + free(lat2print_spat[[4]]) +
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/comp_late_spat_dist_prn_pnl.png"), 
       width = 40, height = 20, units = "cm", dpi = 300)

# bilateral
free(fern_full) + free(bila_spat_pr10_p) + free(bila_spat_pr20_p) + free(bila_spat_pr30_p) + free(bila_spat_pr50_p) +
  free(full_spat_hist) + 
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")

## sim bm trait ####
# lateral
free(fern_full) + free(lat_pr05_p) + free(lat_pr15_p) + free(lat_pr30_p) + free(lat_pr50_p) +
  free(full_tree) + free(lat2print[[1]]) + free(lat2print[[2]]) + free(lat2print[[3]]) + free(lat2print[[4]]) +
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/comp_late_dist_prn_pnl.png"), 
       width = 40, height = 20, units = "cm", dpi = 300)

# bilateral
free(fern_full) + free(marg_pr05_p) + free(marg_pr20_p) + free(marg_pr30_p) + free(marg_pr50_p) +
  free(full_tree) + free(marg2print[[2]]) + free(marg2print[[3]]) + free(marg2print[[4]]) + free(marg2print[[5]]) +
  plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
ggsave(filename = paste0("figures/comp_bilate_dist_prn_pnl.png"), 
       width = 40, height = 20, units = "cm", dpi = 300)

# imputed trees ####

# get  number of trees analyzed ####
# pruned trees
pr_trs_p <- dir("output/pruned_trees", full.names = T)
n_prn_ts <- c()
for(i in 1:length(pr_trs_p)){
  n_prn_ts[i] <- readRDS(pr_trs_p[i]) %>% length()
}

n_prn_ts %>% sum

# imputed trees
im_trs_p <- dir("output/imputed_trees/", full.names = T)
#im_trs_p <- im_trs_p[-which(grepl("_rand_imp_", im_trs_p),im_trs_p)]

n_imp_ts <- c()
for(i in 1:length(im_trs_p)){
  n_imp_ts[i] <- readRDS(im_trs_p[i]) %>% length()
}

n_imp_ts %>% sum()

# compare simulated trees with different node ages ####

# load simulated trees
yng_trees <- readRDS("output/sim_trees/sim_trees_100_age.rds")
old_trees <- readRDS("output/sim_trees/sim_trees_300_age.rds")
cnt_trees <- readRDS("output/sim_trees/sim_trees_217_age.rds")
exp_trees <- readRDS("output/sim_trees/sim_trees_exp_brlen.rds")

typ_tre <-  "rectangular" # "roundrect" # "fan"

# panel
pl_lay <- "
  AA#BBCCDDEE
  AA#BBCCDDEE
  AA#BBCCDDEE
  "

# fern empirical tree
fern_full <- fern_tre %>% 
  ggtree(layout = typ_tre, ladderize = T, size = .2) + 
  geom_nodepoint(size = .1, alpha = .55, col = "red") +
  labs(title = "Empirical tree") + theme_tree2()

dir.create("figures/sim_trees/")

sim_trs_ls <- list()
for(i in 1:length(old_trees)){
  young_tr_p <- yng_trees[[i]] %>% 
    ggtree(layout = typ_tre, ladderize = T, size = .1) + 
    geom_nodepoint(size = .1, alpha = .55, col = "red") +
    labs(title = "Young nodes (100 myr)") + theme_tree2()
  
  old_tr_p <- old_trees[[i]] %>% 
    ggtree(layout = typ_tre, ladderize = T, size = .1) + 
    geom_nodepoint(size = .1, alpha = .55, col = "red") +
    labs(title = "Old nodes (300 myr)") + theme_tree2()
  
  cnt_tr_p <- cnt_trees[[i]] %>% 
    ggtree(layout = typ_tre, ladderize = T, size = .1) + 
    geom_nodepoint(size = .1, alpha = .55, col = "red") +
    labs(title = "Middle nodes (217 mry)") + theme_tree2()
  
  exp_tr_p <- exp_trees[[i]] %>% 
    ggtree(layout = typ_tre, ladderize = T, size = .1) + 
    geom_nodepoint(size = .1, alpha = .55, col = "red") +
    labs(title = "Exp branch lengths") + theme_tree2()
  
  sim_trs_ls[[i]] <- free(fern_full %>% revts) + 
    free(old_tr_p %>% revts) + free(cnt_tr_p %>% revts) + free(young_tr_p %>% revts) +
    free(exp_tr_p %>% revts) +
    plot_layout(guides = "collect", design = pl_lay, axis_titles = "collect", tag_level = "keep")
  
  ggsave(filename = paste0("figures/sim_trees/sim_trees_", typ_tre,"_pnl_", i,".png"), 
         width = 45, height = 25, units = "cm", dpi = 300)
}

# animate plots 
library(magick)
# load figs
img_path <- list.files("figures/sim_trees/", full.names = T)
img_list <- lapply(img_path, image_read)

img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 10, loop = 0, optimize = T)

img_animated %>% image_write(image = ., path = "figures/sim_trees_anima.gif")
