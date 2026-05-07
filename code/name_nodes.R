library(tidyverse) ; library(phytools)
library(taxize)

source("code/pe_samp_funs.R")
#library(pteridocat) library(taxastand)

# load tree
fern_tre <- read.tree("data/nitta-et-al/results/japan_ferns_tree_dated.tre") %>% force.ultrametric
is.binary(fern_tre) ; is.ultrametric(fern_tre) ; is.rooted(fern_tre)
targ_tre <- fern_tre %>% ladderize()

# get taxonomy from genus list ####
x <- targ_tre$tip.label %>% 
  data.frame(tip_names =.) %>% separate(tip_names, sep = "_", into = c("genus", "spp"), extra = "merge") %>%
  pull(genus) %>% unique() %>% sort()

tax_dbase <- "ncbi" # "itis" # "both"  # 
t <- tax_name(sci = c(x), get = c("phylum", "infradivision",
                                  "class",  "infraclass",
                                  "order",  "infraorder",
                                  "subfamily", "family", "superfamily",
                                   "subtribe", "tribe", 
                                  "genus"), db = tax_dbase)

tax_table <- left_join(t, 
                       targ_tre$tip.label %>% data.frame(tip_names = .) %>% 
                         separate(tip_names, sep = "_", into = c("genus", "spp"), extra = "merge") %>% 
                         mutate(tip_names = paste0(genus, "_", spp))
                       )
tax_table %>% filter(is.na(family))
tax_table %>% filter(is.na(subfamily))

tax_table %>% write.csv("output/fern_taxonomy_ncbi.csv", row.names = F)

# name tree nodes based on a taxonomic table ####
m <- read.csv(paste0("output/fern_taxonomy_ncbi.csv"), header = T)

# check that a genus is in the taxonomic table and in the tree
m %>% filter(genus == x[2])
targ_tre$tip.label[grepl(pattern = x[2], x = targ_tre$tip.label) %>% which]

targ_tre$node.label <- rep("", length(targ_tre$node.label)) # replace bootstrap support for empty values

named_tree <- name_nodes(targ_tre, m)
Nnode(named_tree) == length(named_tree$node.label) # check n nodes matches n labels
named_tree$node.label %>% table() %>% data.frame %>% filter(Freq > 1) # no repeated rank names

# replace the second underscore (_) in tip names for a dash (-)
named_tree$tip.label <- sub(pattern = "(\\_.*?)\\_", "\\1-", named_tree$tip.label)

# save tree ####
write.tree(named_tree, "output/nitta_node-named.tre")
named_tree <- read.tree("output/nitta_node-named.tre")

# plot tree and named nodes #####
plot(named_tree, cex = .005, no.margin = T)
nodelabels(text = named_tree$node.label, #adj = c(-.1, .1), 
           cex = .75, col = "red", frame = "none")

# tree from the Fern Tree of Life group Nitta et al. 2022
library(ftolr)
fern_tol <- ft_tree(branch_len = "ultra", consensus = T, backbone = F, 
                   rooted = T, drop_og = T, label_ages = F)
is.binary(fern_tol) ; is.ultrametric(fern_tol)

fern_tol$tip.label <- sub(pattern = "(\\_.*?)\\_", "\\1-", fern_tol$tip.label)

# phylogeny sampling completeness/fraction
fern_com <- read.csv("data/nitta-et-al/japan_ferns_comm_full.csv", header = T) %>% 
  column_to_rownames("grids")

colnames(fern_com) <- gsub(pattern = "\\.", replacement = "-", colnames(fern_com))
colnames(fern_com) <- sub(pattern = "(\\_.*?)\\_", "\\1-", colnames(fern_com))

fern_com %>% colnames %in% named_tree$tip.label %>% table()
round((661/ncol(fern_com))*100, 1)

fern_com %>% colnames %in% fern_tol$tip.label %>% table()
(446/ncol(fern_com))*100
