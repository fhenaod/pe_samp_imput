library(tidyverse)
library(phytools)

source("code/pe_samp_funs.R")

# load tree 
fern_tre <- read.tree("output/nitta_node-named.tre") 
is.binary(fern_tre) ; is.ultrametric(fern_tre)
targ_tre <- fern_tre

# toy example
topo <- read.tree(text = "((A,B),(C,(D,E)));")

trees <- sim_fixed_topology_ages(topo, 
                                      nrep = 100, root_age = 20, 
                                      mean_node_age = 5, # 5 - 7 
                                      sd_node_age = 0.5, 
                                      min_edge = 1e-6, 
                                      dist = "exponential")

par(mfrow = c(2, 2))
sapply(trees, plot, cex = .5, no.margin = T)
par(mfrow = c(1,1))

# ferns of japan example
topo <- targ_tre
topo$root.time <- branching.times(topo) %>% max() 
ape::node.depth.edgelength(topo) %>% summary()

# simulate trees with fixed conditions ####
dir.create("output/sim_trees/", recursive = T)

# trees with very young nodes
v_yng_trees <- sim_fixed_topology_ages(topo, 
                                     nrep = 100, root_age = topo$root.time, 
                                     mean_node_age = 20, 
                                     sd_node_age = 25, 
                                     min_edge = 1e-6, 
                                     dist = "gamma")
v_yng_trees %>% saveRDS("output/sim_trees/sim_trees_20_age.rds")
plot(v_yng_trees[[10]], show.tip.label = F, no.margin = F)
axisPhylo()

## trees with younger nodes
yng_trees <- sim_fixed_topology_ages(topo, 
                                      nrep = 100, root_age = topo$root.time, 
                                      mean_node_age = 100, 
                                      sd_node_age = 25, 
                                      min_edge = 1e-6, 
                                      dist = "gamma")
yng_trees %>% saveRDS("output/sim_trees/sim_trees_100_age.rds")
plot(yng_trees[[10]], show.tip.label = F, no.margin = F)
axisPhylo()

## trees with older nodes
old_trees <- sim_fixed_topology_ages(topo, 
                                          nrep = 100, root_age = topo$root.time, 
                                          mean_node_age = 300, 
                                          sd_node_age = 25, 
                                          min_edge = 1e-6, 
                                          dist = "gamma")
old_trees %>% saveRDS("output/sim_trees/sim_trees_300_age.rds")

plot(old_trees[[1]], show.tip.label = F, no.margin = F)
axisPhylo()

## trees with node ages centered
cnt_trees <- simulate_fixed_topology_ages(topo, 
                                          nrep = 100, root_age = topo$root.time, 
                                          mean_node_age = round(topo$root.time/2), 
                                          sd_node_age = 25, 
                                          min_edge = 1e-6, 
                                          dist = "gamma")
cnt_trees %>% saveRDS("output/sim_trees/sim_trees_217_age.rds")

plot(cnt_trees[[1]], show.tip.label = F, no.margin = F)
axisPhylo()

# sim trees with fixed topology, branch lengths drawn from exp dist
emp_rate <- 1 / mean(topo$edge.length)
hist(topo$edge.length) ; abline(v = 1/emp_rate, col = "coral")

exp_trees <- sim_fixed_topo_exp_br(topo, 
                                   nrep = 100, 
                                   rate = emp_rate*.2) 
exp_trees <- lapply(X = exp_trees, function(X)set_tree_root_age(tree = X, root_age = topo$root.time)) 
exp_trees %>% saveRDS("output/sim_trees/sim_trees_exp_brlen.rds")

 plot(exp_trees[[90]] %>% ladderize(right = F), show.tip.label = F, no.margin = F)
axisPhylo()

# simulate and compare trees after node modification #####
# sim tree
raw_tr <- phytools::pbtree(b = .8, d = .05, n = 30, scale = 10) %>% 
  geiger::drop.extinct()
raw_tr$root.time <- branching.times(raw_tr) %>% max() # set root age

# check tree
plot(raw_tr, no.margin = F)
nodelabels(frame = "none",  col = "red", adj = c(1.5, -.5), cex = .75)
edgelabels(text = round(raw_tr$edge.length, 3), frame = "none", col = "blue", cex = 0.75, adj = c(1, 1.2))
tiplabels(frame = "none", col = "darkgreen", cex = 0.5, adj = c(-3, 1))
axisPhylo()

# build dataframe with node ages, height, depth
t_ages_df <- data.frame(n_tip = 1:c(Nnode(raw_tr) + Ntip(raw_tr)),
                        node_depth = node.depth.edgelength(raw_tr),
                        node_age = Claddis::date_nodes(raw_tr) %>% round(5), 
                        stem_len = c(raw_tr$edge.length[raw_tr$edge[, 2] <= Ntip(raw_tr)], 0, raw_tr$edge.length[raw_tr$edge[, 2] > Ntip(raw_tr)]),
                        type = c(rep("terminal", Ntip(raw_tr)), rep("internal", Nnode(raw_tr)))
                        ) 
head(t_ages_df)

# Change node ages
tr <- phyloch::shiftNodeAge(raw_tr, node = 32, shift = -5)
tr <- phyloch::shiftNodeAge(tr, node = 33, shift = -5)
tr <- phyloch::shiftNodeAge(tr, node = 42, shift = 2)
tr <- phyloch::shiftNodeAge(tr, node = 57, shift = -1)

is.ultrametric(tr)

# plot compared trees
phyloch::compare.phylo(raw_tr, tr)
nodelabels(frame = "none", col = "red", adj = c(1.5, -.5), cex = .75)
axisPhylo()

# plot besides each other
par(mfrow = c(1, 2))
plot(raw_tr, no.margin = F)
nodelabels(frame = "none", col = "red", adj = c(1.5, -.5), cex = .75)

plot(tr, no.margin = F)
nodelabels(frame = "none", col = "red", adj = c(1.5, -.5), cex = .75)
par(mfrow = c(1, 1))
