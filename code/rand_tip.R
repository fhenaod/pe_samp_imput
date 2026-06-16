library(randtip)
library(parallel)
library(canaper)

# load data ####
fern_tre <- read.tree("output/nitta_node-named.tre") 

late_fls <- list.files(paste0("output/pruned_trees/fern"), "lateral_spat",  all.files = T, full.names = T)
late_frc <- readr::parse_number(late_fls)

j = 4
targ_tr_ls <- readRDS(late_fls[[j]])
targ_tre <- targ_tr_ls

# backbone mode ####
info_bb <- build_info(species = fern_tre$tip.label, tree = targ_tre, 
                        find.ranks = T, db = "ncbi", mode = "backbone")

check_bb <- check_info(info = info_bb, tree = targ_tre, parallelize = T, ncores = 3)
input_bb <- info2input(info = info_bb, tree = targ_tre, parallelize = T, ncores = 3)

# random insertion
imp_tr_bb <- rand_tip(input = input_bb, tree = targ_tre, 
                      rand.type = "random", polyphyly.scheme = "largest",
                      use.paraphyletic = T, use.singleton = F, 
                      use.stem = T, 
                      respect.mono = T, respect.para = T,
                      prob = T, prune = T, clump.puts = T, 
                      forceultrametric = T)

tip_col <- put_tip_col(newtree = imp_tr_bb, oldtree = targ_tre, 
                       placed.col = "dark grey", 
                       put.col = "red")
plot.phylo(imp_tr_bb, tip.color = tip_col, cex = 0.25)

# polytomy insertion
# does not work !!!
imp_tr_bb2 <- rand_tip(input = input_bb, tree = targ_tre,
                       rand.type = "polytomy", polyphyly.scheme = "largest",
                       use.paraphyletic = T, use.singleton = T, 
                       use.stem = T, 
                       respect.mono = F, respect.para = F,
                       prob = T, prune = T, clump.puts = T, 
                       forceultrametric = T)

tip_col <- put_tip_col(newtree = imp_tr_bb2, oldtree = targ_tre, 
                       placed.col = "dark grey", 
                       put.col = "red")
plot.phylo(imp_tr_bb, tip.color = tip_col, cex = 0.25)

# taxon list mode #####
# load taxonomy table
tax_tbl <- read.csv("output/fern_taxonomy_ncbi.csv") %>% 
  mutate(tip_names = sub(pattern = "(\\_.*?)\\_", "\\1-", tip_names))

tbl2merg <- tax_tbl %>% distinct(genus, .keep_all = T) %>% 
  select(genus, #subfamily,
         family, order, class)

info_lt <- build_info(species = fern_tre$tip.label, #tree = targ_tre, 
                      find.ranks = F, mode = "list")

info_cols <- info_lt %>% names()
info_lt2 <- left_join(info_lt %>% select(#-subfamily,
                                         -family, -order, -class), 
                      tbl2merg, 
                      by = join_by("genus" == "genus")
                      ) %>% 
  select(info_cols) #%>% mutate(rand.type = "polytomy")

info_lt2 %>% filter(is.na(family))

check_lt <- check_info(info = info_lt2, tree = targ_tre, parallelize = T, ncores = 3)
input_lt <- info2input(info = info_lt2, tree = targ_tre, parallelize = T, ncores = 3)

# random insertion
imp_tr_lt1 <- rand_tip(input = input_lt, tree = targ_tre,  
                      rand.type = "random", polyphyly.scheme = "largest",
                      use.paraphyletic = T, use.singleton = T, 
                      use.stem = T, 
                      respect.mono = T, respect.para = T,
                      prob = T, prune = T, clump.puts = T, 
                      forceultrametric = T)
imp_tr_lt1
tip_col <- put_tip_col(newtree = imp_tr_lt1, oldtree = targ_tre, 
                       placed.col = "dark grey", 
                       put.col = "red")
plot.phylo(imp_tr_lt1, tip.color = tip_col, cex = 0.25#, type = "fan"
           )

# loop ##
n_reps <- 100

rnd_tp_tre_ls <- list()
for(i in seq(n_reps)){
  rnd_tp_tre_ls[[i]] <- rand_tip(input = input_lt, tree = targ_tre,  
              rand.type = "random", polyphyly.scheme = "largest",
              use.paraphyletic = T, use.singleton = T, 
              use.stem = T, 
              respect.mono = T, respect.para = T,
              prob = T, prune = T, clump.puts = T, 
              forceultrametric = T, verbose = F)
  print(paste0("Tree: ", i))
}
sapply(rnd_tp_tre_ls, Ntip)
sapply(rnd_tp_tre_ls, is.binary)
sapply(rnd_tp_tre_ls, is.ultrametric)

all.equal(rnd_tp_tre_ls[[1]], rnd_tp_tre_ls[[3]],
          use.edge.length = T, use.tip.label = T)

rnd_tp_tre_ls %>% saveRDS("output/imputed_trees/fern_lat_spat_pr_randtip_0.5.rds")

# polytomy insertion
# does not work !!!
imp_tr_lt2 <- rand_tip(input = input_lt, tree = targ_tre,  
                       rand.type = "polytomy", polyphyly.scheme = "largest",
                       use.paraphyletic = T, use.singleton = F, 
                       use.stem = T, 
                       respect.mono = F, respect.para = F,
                       prob = T, prune = F, clump.puts = T, 
                       forceultrametric = T)
imp_tr_lt2

tip_col <- put_tip_col(newtree = imp_tr_lt2, oldtree = targ_tre, 
                       placed.col = "dark grey", 
                       put.col = "red")
plot.phylo(imp_tr_lt2, tip.color = tip_col, cex = 0.25, type = "fan")

fern_tre$tip.label %in% imp_tr_lt1$tip.label %>% table
