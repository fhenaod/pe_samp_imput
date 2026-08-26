library(ggplot2)
library(patchwork)
library(grid)
library(cowplot)

# --- helper theme ---
theme_panel <- function() {
  theme_void(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# --- helper to make labeled panels ---
make_panel <- function(letter, title, subtitle, fill) {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = fill, color = NA, alpha = 0.12) +
    annotate("label", x = 0.08, y = 0.92, label = letter,
             size = 6, fontface = "bold", label.size = 0,
             fill = fill, color = "white") +
    annotate("text", x = 0.22, y = 0.92, label = title,
             hjust = 0, size = 5.2, fontface = "bold") +
    annotate("text", x = 0.06, y = 0.72, label = subtitle,
             hjust = 0, vjust = 1, size = 3.8, lineheight = 0.95) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_panel()
}

# simple text ####

# --- A: input data ---
pA <- make_panel(
  "A", "Input data",
  "Spatial fern occurrences from Japan\n673 taxa, 1239 grid cells, ~20 x 20 km\nTime-calibrated phylogeny: 663 terminals",
  "#4DAA57"
) +
  annotate("segment", x = 0.10, xend = 0.42, y = 0.45, yend = 0.45,
           linewidth = 1.1, color = "#2E7D32", arrow = arrow(length = unit(0.15, "inches"))) +
  annotate("segment", x = 0.58, xend = 0.90, y = 0.45, yend = 0.45,
           linewidth = 1.1, color = "#2E7D32", arrow = arrow(length = unit(0.15, "inches"))) +
  annotate("text", x = 0.26, y = 0.52, label = "Spatial data", size = 4, color = "#2E7D32") +
  annotate("text", x = 0.74, y = 0.52, label = "Phylogeny", size = 4, color = "#2E7D32")

# --- B: pruning ---
pB <- make_panel(
  "B", "Pruning schemes",
  "Spatial pruning: range-size tails\nPhylogenetic pruning: BM trait tails\nBaseline: MCAR random pruning",
  "#4C78A8"
) +
  annotate("text", x = 0.50, y = 0.57, label = "5%, 15%, 30%, 50%", size = 4) +
  annotate("text", x = 0.50, y = 0.23, label = "Bilateral: 10%, 30%, 50%", size = 4) +
  annotate("text", x = 0.18, y = 0.5, label = "Spatial", size = 4, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.5, label = "Phylogenetic", size = 4, fontface = "bold") +
  annotate("text", x = 0.82, y = 0.5, label = "MCAR", size = 4, fontface = "bold") +
  annotate("segment", x = 0.12, xend = 0.88, y = 0.35, yend = 0.35,
           linewidth = 1.1, color = "#4C78A8", arrow = arrow(length = unit(0.15, "inches")))

# --- C: imputation ---
pC <- make_panel(
  "C", "Imputation",
  "Deleted tips reinserted using taxonomy-based methods\n1) Crown-node polytomy\n2) Random below crown node\n3) Random edge anywhere",
  "#F28E2B"
) +
  annotate("text", x = 0.20, y = 0.5, label = "Polytomy", size = 4, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.5, label = "Below crown", size = 4, fontface = "bold") +
  annotate("text", x = 0.80, y = 0.5, label = "Random edge", size = 4, fontface = "bold") +
  annotate("segment", x = 0.15, xend = 0.85, y = 0.38, yend = 0.38,
           linewidth = 1.1, color = "#F28E2B", arrow = arrow(length = unit(0.15, "inches")))

# --- D: analyses ---
pD <- make_panel(
  "D", "Spatial phylogenetic analyses",
  "Alpha-phylodiversity: PD, PE, RPD, RPE\n\nObserved and SES metrics from 1000 null communities [SES = (obs - null) / sd] \n\nCANAPE classifications from PE and RPE: Neo | Paleo | Mixed",
  "#8E63C7"
) +
  annotate("text", x = 0.18, y = 0.52, label = "PD", size = 5, fontface = "bold") +
  annotate("text", x = 0.38, y = 0.52, label = "PE", size = 5, fontface = "bold") +
  annotate("text", x = 0.58, y = 0.52, label = "RPD", size = 5, fontface = "bold") +
  annotate("text", x = 0.78, y = 0.52, label = "RPE", size = 5, fontface = "bold") #+
  #annotate("text", x = 0.50, y = 0.28, label = "SES = (obs - null) / sd", size = 4) +annotate("text", x = 0.50, y = 0.16, label = "CANAPE: Neo | Paleo | Mixed", size = 4)

# --- E: performance ---
pE <- make_panel(
  "E", "Performance comparisons",
  "Compare pruned/imputed datasets to the full dataset\nMCC for endemicity recovery\nMPE for bias in phylodiversity metrics",
  "#6B6B6B"
) +
  annotate("text", x = 0.22, y = 0.55, label = "Confusion\nmatrix", size = 4, fontface = "bold") +
  annotate("text", x = 0.65, y = 0.58, label = "MCC\n-1 to +1", size = 4, fontface = "bold") +
  annotate("text", x = 0.65, y = 0.28, label = "MPE\nover vs under", size = 4, fontface = "bold") +
  annotate("segment", x = 0.35, xend = 0.58, y = 0.55, yend = 0.55,
           linewidth = 1.1, color = "#6B6B6B", arrow = arrow(length = unit(0.15, "inches"))) +
  annotate("segment", x = 0.35, xend = 0.58, y = 0.35, yend = 0.35,
           linewidth = 1.1, color = "#6B6B6B", arrow = arrow(length = unit(0.15, "inches")))

# --- combine ---
fig <- (pA | pB | pC) / (pD | pE) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1)) +
  plot_annotation(
    title = "'Henao-Diaz. et al. - The effects of missing taxa and phylogenetic imputation on spatial phylogenetic studies' ",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )
fig
ggsave("fern_methods_schematic.png", fig, width = 18, height = 10, dpi = 300, bg = "white")

# graphical illustrated ####
# ---------- theme ----------
theme_blank <- function() {
  theme_void(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# ---------- reusable panel ----------
make_panel <- function(letter, title, subtitle, fill) {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = fill, alpha = 0.10, color = NA) +
    annotate("label", x = 0.08, y = 0.92, label = letter,
             label.size = 0, fill = fill, color = "white",
             fontface = "bold", size = 6) +
    annotate("text", x = 0.22, y = 0.92, label = title,
             hjust = 0, fontface = "bold", size = 5.4) +
    annotate("text", x = 0.06, y = 0.76, label = subtitle,
             hjust = 0, vjust = 1, size = 3.8, lineheight = 0.96) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_blank()
}

# ---------- icon placeholders ----------
fern_icon <- function(x0, y0, scale = 0.10, col = "#2E7D32") {
  ggplot() +
    annotate("segment", x = x0, xend = x0, y = y0 - scale, yend = y0 + scale,
             linewidth = 1.2, color = col) +
    annotate("segment", x = x0, xend = x0 + scale*0.55, y = y0 + scale*0.2, yend = y0 + scale*0.55,
             linewidth = 1.1, color = col) +
    annotate("segment", x = x0, xend = x0 - scale*0.55, y = y0 + scale*0.3, yend = y0 + scale*0.65,
             linewidth = 1.1, color = col) +
    annotate("segment", x = x0, xend = x0 + scale*0.50, y = y0 - scale*0.05, yend = y0 - scale*0.35,
             linewidth = 1.1, color = col) +
    annotate("segment", x = x0, xend = x0 - scale*0.55, y = y0 - scale*0.15, yend = y0 - scale*0.45,
             linewidth = 1.1, color = col) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

tree_icon <- function() {
  ggplot() +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.18, yend = 0.80,
             linewidth = 1.3, color = "#333333") +
    annotate("segment", x = 0.50, xend = 0.34, y = 0.55, yend = 0.72,
             linewidth = 1.1, color = "#333333") +
    annotate("segment", x = 0.50, xend = 0.66, y = 0.55, yend = 0.72,
             linewidth = 1.1, color = "#333333") +
    annotate("segment", x = 0.50, xend = 0.34, y = 0.42, yend = 0.30,
             linewidth = 1.1, color = "#333333") +
    annotate("segment", x = 0.50, xend = 0.66, y = 0.42, yend = 0.30,
             linewidth = 1.1, color = "#333333") +
    annotate("point", x = c(0.34, 0.66, 0.34, 0.66), y = c(0.72, 0.72, 0.30, 0.30),
             size = 3.2, color = "#333333") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

map_placeholder <- function() {
  ggplot() +
    annotate("rect", xmin = 0.20, xmax = 0.80, ymin = 0.22, ymax = 0.78,
             fill = "white", color = "#4DAA57", linewidth = 1.1) +
    annotate("segment", x = seq(0.28, 0.72, by = 0.12), xend = seq(0.28, 0.72, by = 0.12),
             y = 0.22, yend = 0.78, color = "#4DAA57", alpha = 0.45) +
    annotate("segment", x = 0.20, xend = 0.80, y = seq(0.30, 0.70, by = 0.10), yend = seq(0.30, 0.70, by = 0.10),
             color = "#4DAA57", alpha = 0.45) +
    annotate("point", x = c(0.32, 0.50, 0.65), y = c(0.65, 0.45, 0.32),
             color = "#4DAA57", size = 2.6) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

confusion_icon <- function() {
  ggplot() +
    annotate("rect", xmin = 0.20, xmax = 0.78, ymin = 0.25, ymax = 0.75,
             fill = "white", color = "#666666", linewidth = 1.0) +
    annotate("segment", x = 0.49, xend = 0.49, y = 0.25, yend = 0.75, color = "#666666") +
    annotate("segment", x = 0.20, xend = 0.78, y = 0.50, yend = 0.50, color = "#666666") +
    annotate("text", x = 0.34, y = 0.62, label = "TP", size = 5, fontface = "bold") +
    annotate("text", x = 0.62, y = 0.62, label = "FP", size = 5, fontface = "bold") +
    annotate("text", x = 0.34, y = 0.38, label = "FN", size = 5, fontface = "bold") +
    annotate("text", x = 0.62, y = 0.38, label = "TN", size = 5, fontface = "bold") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

# ---------- panels ----------
pA <- make_panel(
  "A", "Input data",
  "Japanese fern spatial data + time-calibrated phylogeny",
  "#4DAA57"
) + map_placeholder() +
  annotate("text", x = 0.22, y = 0.28, label = "Spatial grid", size = 4, color = "#2E7D32") +
  annotate("text", x = 0.70, y = 0.28, label = "Phylogeny", size = 4, color = "#2E7D32") +
  annotation_custom(ggplotGrob(fern_icon(0.5, 0.45, 0.12, "#2E7D32")))

pB <- make_panel(
  "B", "Pruning schemes",
  "Spatial pruning, phylogenetic pruning, and MCAR random pruning",
  "#4C78A8"
) +
  annotate("segment", x = 0.15, xend = 0.85, y = 0.45, yend = 0.45,
           linewidth = 1.1, color = "#4C78A8",
           arrow = arrow(length = unit(0.15, "inches"))) +
  annotate("text", x = 0.22, y = 0.62, label = "Range-size tails", size = 4, fontface = "bold") +
  annotate("text", x = 0.52, y = 0.62, label = "BM trait tails", size = 4, fontface = "bold") +
  annotate("text", x = 0.80, y = 0.62, label = "MCAR", size = 4, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.30, label = "5%, 15%, 30%, 50%\nBilateral: 10%, 30%, 50%", size = 3.8)

pC <- make_panel(
  "C", "Imputation",
  "Deleted tips reinserted by taxonomy-based strategies",
  "#F28E2B"
) +
  annotate("text", x = 0.22, y = 0.60, label = "Crown-node\npolytomy", size = 4, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.60, label = "Below-crown\ninsertion", size = 4, fontface = "bold") +
  annotate("text", x = 0.78, y = 0.60, label = "Random-edge\ninsertion", size = 4, fontface = "bold") +
  annotate("segment", x = 0.15, xend = 0.85, y = 0.38, yend = 0.38,
           linewidth = 1.1, color = "#F28E2B",
           arrow = arrow(length = unit(0.15, "inches")))

pD <- make_panel(
  "D", "Spatial phylogenetic analyses",
  "PD, PE, RPD, RPE; SES from null communities; CANAPE classification",
  "#8E63C7"
) +
  annotate("segment", x = 0.15, xend = 0.85, y = 0.52, yend = 0.52,
           linewidth = 1.0, color = "#8E63C7") +
  annotate("text", x = 0.18, y = 0.60, label = "PD", size = 5, fontface = "bold") +
  annotate("text", x = 0.36, y = 0.60, label = "PE", size = 5, fontface = "bold") +
  annotate("text", x = 0.54, y = 0.60, label = "RPD", size = 5, fontface = "bold") +
  annotate("text", x = 0.72, y = 0.60, label = "RPE", size = 5, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.32, label = "1000 null communities\ncurveball randomization", size = 4) +
  annotate("text", x = 0.50, y = 0.18, label = "CANAPE: neoendemism / paleoendemism / mixed", size = 4)

pE <- make_panel(
  "E", "Performance comparisons",
  "MCC for endemicity recovery; MPE for bias in metric estimates",
  "#6B6B6B"
) +
  confusion_icon() +
  annotate("text", x = 0.72, y = 0.65, label = "MCC\n-1 to +1", size = 4.2, fontface = "bold") +
  annotate("text", x = 0.72, y = 0.35, label = "MPE\nover / under", size = 4.2, fontface = "bold") +
  annotate("segment", x = 0.45, xend = 0.63, y = 0.63, yend = 0.63,
           linewidth = 1.1, color = "#6B6B6B",
           arrow = arrow(length = unit(0.15, "inches"))) +
  annotate("segment", x = 0.45, xend = 0.63, y = 0.35, yend = 0.35,
           linewidth = 1.1, color = "#6B6B6B",
           arrow = arrow(length = unit(0.15, "inches")))

# ---------- compose ----------
fig <- (pA | pB | pC) / (pD | pE) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Schematic workflow for fern pruning, imputation, and evaluation",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  )
fig
ggsave("fern_methods_mockup.png", fig, width = 18, height = 10, dpi = 300, bg = "white")
