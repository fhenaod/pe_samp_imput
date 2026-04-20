# 
# functions 


# classifies continuous values into five percentiles
lat_funct=function(x){cut(x, 
                          breaks = c(-Inf, quantile(x, lat_samp_prop), +Inf), 
                          labels = c("below.05", ".05to.15", ".15to.3", ".3to.5", "above.5"))
}

# classifies continuous values into percentiles
margin_funct=function(x){cut(x, 
                             breaks = c(-Inf, quantile(x, marg_samp_prop), +Inf),
                             labels = cat2filt)
}

# classifies significance values from canaper in three catagories
signif2cats=function(x){
  ifelse(x %in% c("< 0.01", "< 0.025"), "sig_low", 
         ifelse(x %in% c("> 0.975", "> 0.99"), "sig_high", "non_sig"))
}

thres2cats=function(x){
  ifelse(x >= ., "sig_low", 
         ifelse(x %in% c("> 0.975", "> 0.99"), "sig_high", "non_sig"))
}

# retrieves a result category based on comparing the type of result from the treatment and the full dataset
res_type=function(sign_cat, sign_cat_emp){
  ifelse(sign_cat == "sig_low" & sign_cat_emp == "sig_low", "true_positive", 
         ifelse(sign_cat == "sig_high" & sign_cat_emp == "sig_high", "true_positive", 
                ifelse(sign_cat == "non_sig" & sign_cat_emp %in% c("sig_low","sig_high"), "type_2", 
                       ifelse(sign_cat == "non_sig" & sign_cat_emp == "non_sig", "true_negative", 
                              ifelse(sign_cat == "n_a" & sign_cat_emp %in% c("sig_low","sig_high"), "inconclusive",
                                     ifelse(sign_cat %in% c("sig_low","sig_high") & sign_cat_emp == "non_sig", "type_1", 
                                            ifelse(sign_cat == "sig_low" & sign_cat_emp == "sig_high" |
                                                   sign_cat == "sig_hig" & sign_cat_emp == "sig_low", "switched", "inconclusive")
                                            )
                                     )
                              )
                       )
                )
         )
  }

# names nodes in a tree based on a taxonomic table for that tree 
name_nodes=function(tree, m_tax){
  require(ape)
  if(Nnode(tree)==length(tree$node.label)){} else {print("Warning ! Unequal n-nodes and n-labels")}
  
  m_tips <- tree$tip.label %>% data.frame(tip_names =.) %>% 
    mutate(tip_names,  tips = tip_names) %>%
    separate(tips, sep = "_", into = c("genus", "spp"), extra = "merge")
  
  m <- left_join(m_tips, m_tax, by = join_by(tip_names, genus, spp)) 
  
  gens2name <- m %>% pull(genus) %>% unique
  for(i in 1:length(gens2name)){
    gen_t <- filter(m, genus == gens2name[i]) %>% pull(tip_names)
    if(length(gen_t) != 1) {
      mrca_node <- getMRCA(tree, tip = gen_t)
      tree$node.label[c(mrca_node-Ntip(tree))] <- gens2name[i]
    } else {print(paste0(toupper(gens2name[i])," skipped: only one tip to name"))}
  }
  
  fams2name <- m %>% pull(family) %>% unique
  for(i in 1:length(fams2name)){
    fam_t <- filter(m, family == fams2name[i]) %>% pull(tip_names)
    if(length(fam_t) != 1) {
      mrca_node <- getMRCA(tree, tip = fam_t)
      tree$node.label[c(mrca_node-Ntip(tree))] <- fams2name[i]
    } else {print(paste0(toupper(fams2name[i])," skipped: only one tip to name"))}
  }
  
  ords2name <- m %>% pull(order) %>% unique
  for(i in 1:length(ords2name)){
    ord_t <- filter(m, order == ords2name[i]) %>% pull(tip_names)
    if(length(ord_t) > 1) {
      mrca_node <- getMRCA(tree, tip = ord_t)
      tree$node.label[c(mrca_node-Ntip(tree))] <- ords2name[i]
    } else {print(paste0(toupper(fams2name[i]), " skipped: only one tip to name"))}
  }
  
  class2name <- m %>% pull(class) %>% unique
  for(i in 1:length(class2name)){
    cla_t <- filter(m, class == class2name[i]) %>% pull(tip_names)
    if(length(cla_t) > 1) {
      mrca_node <- getMRCA(tree, tip = cla_t)
      tree$node.label[c(mrca_node-Ntip(tree))] <- class2name[i]
    } else {print(paste0(toupper(class2name[i]), "skipped: only one tip to name"))}
  }
  return(tree)
}

# get label number from tree to impute a genus based on taxonomic table and a node-named tree
node2where=function(gen, tree){
  gen_tab <- filter(tax_tbl, genus == gen)[1,]
  
  nd_ls <- list(
    family = grepl(pattern = paste0("^", gen_tab$family, "$"), x = tree$node.label) %>% which, 
    order = grepl(pattern = paste0("^", gen_tab$order , "$"), x = tree$node.label) %>% which,
    class = grepl(pattern = paste0("^", gen_tab$class , "$"), x = tree$node.label) %>% which
    )
  
  if(!is_empty(nd_ls$family)){
    nd_where <- nd_ls$family ; print("Imputed at family level")
  } else {
    if(!is_empty(nd_ls$order)) {
      nd_where <- nd_ls$order ; print("Imputed at order level")
    } else {
      if(!is_empty(nd_ls$class)) {
        nd_where <- nd_ls$class ; print("Imputed at class level")
      } else {
        nd_where <- 2 ; cat("Imputed close to root node")
      }
    }
  }
  
  if(length(nd_where) > 1) {
    print(paste0("WARNING !! ", length(nd_where), " matched named nodes, selected the closest to the root"))
    nd_where <- nd_where[1]
  } else {}
  return(nd_where)
}

# functions to estimate Robinson–Foulds distances, adjusting for phylo-info content
# by tree 
get_tree_metrics=function(tree, sum){
  require(TreeDist, quietly = T)
  data.frame(splitwise_info = SplitwiseInfo(tree, sum = sum),
             cluster_entropy = ClusteringEntropy(tree, sum = sum),
             clustering_info = ClusteringInfo(tree, sum = sum)
             )
}

# comparing between trees
get_t2t_metrics=function(ref_tree, compar_tree, sum, sim, norm){
  require(TreeDist, quietly = T)
  data.frame(rob.fou_norm = RobinsonFoulds(ref_tree, compar_tree, similarity = sim, normalize = norm),
             rob.fou_info = RobinsonFouldsInfo(ref_tree, compar_tree, similarity = sim, normalize = norm),
             shared.phy_info = SharedPhylogeneticInfo(ref_tree, compar_tree, normalize = norm),
             phylo_info_dist = PhylogeneticInfoDistance(ref_tree, compar_tree, normalize = norm)
             )
}

# simulate trees given average node age fixed topology
# gamma may produce more asymmetric tree shapes
# shape controls the concentration of node ages, larger make them more clustered while lower make them more spread
sim_fixed_topology_ages <- function(topology, nrep, root_age,
                                         mean_node_age, sd_node_age,
                                         min_edge, 
                                         dist = c("normal", "gamma", "exponential")) {
  dist <- match.arg(dist)
  out <- vector("list", nrep)
  
  ntip <- length(topology$tip.label)
  nnodes <- topology$Nnode
  internal_nodes <- (ntip + 1):(ntip + nnodes)
  
  for (i in seq_len(nrep)) {
    tr <- topology
    ages <- numeric(ntip + nnodes)
    ages[1:ntip] <- 0
    ages[ntip + 1] <- root_age
    
    k <- length(internal_nodes) - 1
    
    if (k > 0) {
      x <- switch(dist,
                  normal = rnorm(k, mean = mean_node_age, sd = sd_node_age), 
                  gamma  = rgamma(k, shape = (mean_node_age / sd_node_age)^2, scale = sd_node_age^2 / mean_node_age),
                  exponential = rexp(k, rate = 1 / mean_node_age)
      )
      
      x <- pmin(pmax(x, min_edge), root_age - min_edge)
      x <- sort(x, decreasing = T)
      
      ages[internal_nodes[-1]] <- x
    }
    
    bl <- ages[tr$edge[, 1]] - ages[tr$edge[, 2]]
    if (any(bl <= 0)) {
      stop("WARNING!! Non-positive branch lengths") # if stops adjust mean_node_age, sd_node_age, or root_age.
    }
    tr$edge.length <- bl
    out[[i]] <- tr
  }
  out
}

# simulates trees with a fixed topology, branch lengths drawn from exponential distribution
# smaller rate retrieve longer branches (spread), while larger rate shorter branches (compact)
sim_fixed_topo_exp_br <- function(topology, nrep , rate) {
  out <- vector("list", nrep)
  for (i in seq_len(nrep)) {
    tr <- topology
    bl <- rexp(nrow(tr$edge), rate = rate)
    tr$edge.length <- pmax(bl, 1e-6)
    tr <- multi2di(tr)
    tr <- phytools::force.ultrametric(tree = tr, method = "nnls", message = F)
    out[[i]] <- tr
  }
  out
}

# scale a tree to a given root age
set_tree_root_age <- function(tree, root_age) {
  nh <- nodeHeights(tree) # current root age (height of tree)
  current_root_age <- max(nh[, 2]) # height above root for each edge end
  
  scale <- root_age / current_root_age # scale branch lengths
  tree$edge.length <- tree$edge.length * scale
  
  tree
}
