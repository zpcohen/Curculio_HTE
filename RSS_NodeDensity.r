
dt <- fread("/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/C1_13_Feb26_input.tsv", col.names = c("source","node_id","pos","size","ref"))

# ---- parse source (species A/B/C) ----
dt[, source_species := tstrsplit(source, "#", fixed = TRUE)[[1]]]
dt[, sp_ABC := fcase(
  source_species == "Ccaryae", "A",
  source_species == "Cnanu",   "B",
  source_species == "Cglandium", "C",
  default = NA_character_
)]

head(dt)
# ---- parse reference coordinate (always Cglandium) ----
dt[, c("ref_species","ref_mid","ref_chrom_label") := tstrsplit(ref, "#", fixed = TRUE)]
dt[, chrom := str_remove(ref_chrom_label, "^Chrom")]
dt[, chrom := ifelse(chrom == "" | is.na(chrom), ref_chrom_label, chrom)]

win_size <- 1e6
dt[, win := (pos) %/% win_size]   # 0-based bin index

# ---- count UNIQUE node_id per species x chrom x 1Mb bin ----
counts <- dt[, .(n_nodes = uniqueN(node_id)), by = .(sp_ABC, chrom, win)]
head(counts)
# ---- make bins contiguous (fill missing with zero counts) ----
# complete all bins observed per (species, chrom)
setorder(counts, sp_ABC, chrom, win)
# coerce and check
counts[, win := as.integer(as.character(win))]
counts[ , .N, by=.(sp_ABC, chrom)][order(-N)][1:5]      # quick sanity
counts[!is.finite(win) | is.na(win)]     

counts <- counts[, {
  full <- data.table(win = seq(min(win), max(win)))
  merge(full, .SD, by = "win", all.x = TRUE)
}, by = .(sp_ABC, chrom)]

counts[is.na(n_nodes), n_nodes := 0]
counts[, win_start := win * win_size]
counts[, win_mid_Mb := (win_start + 0.5 * win_size) / 1e6]


# ---- rolling 5×1Mb root-sum-of-squares (RSS) per species × chrom ----
# trailing window: uses bins w, w-1, w-2, w-3, w-4  (total span 5 Mb)
counts <- counts[
  , PBSa_5Mb := slide_dbl(
    n_nodes,
    ~ sqrt(sum((.x)^2)),
    .before = 4, .after = 0, .complete = TRUE
  ),
  by = .(sp_ABC, chrom)
]


##### ONE AXIS, ONE PLOT
library(data.table)
library(ggplot2)
library(scales)

#
df <- as.data.table(counts)

# === knobs ===
yvar         <- "PBSa_5Mb"             # or "PBSa_5Mb"
percentiles  <- c(.90, .95, .99)
band_alpha   <- c(.10, .15, .20)      # alpha for 90/95/99
species_order <- c("C. caryae","C. glandium","C. nanulus")  # optional
use_log      <- FALSE                 # TRUE if dynamic range is huge

# === prep ===
stopifnot(all(c("species_label","chrom","win_mid_Mb", yvar) %in% names(df)))

# ensure numeric & finite
df[, chrom_num := as.integer(as.character(chrom))]
df <- df[is.finite(get(yvar)) & is.finite(win_mid_Mb)]
if (!is.null(species_order)) df[, species_label := factor(species_label, levels = species_order)]

# chromosome lengths and cumulative offsets (Mb)
setorder(df, chrom_num, win_mid_Mb)
chr_len <- df[, .(chr_len_Mb = max(win_mid_Mb)), by = chrom_num][order(chrom_num)]
chr_len[, chrom_start_Mb := shift(cumsum(chr_len_Mb), type = "lag", fill = 0)]
chr_len[, chrom_end_Mb   := chrom_start_Mb + chr_len_Mb]

# attach cumulative x
df <- chr_len[df, on = .(chrom_num)]
df[, cum_Mb := chrom_start_Mb + win_mid_Mb]

# per-chromosome & genome means (per species)
means_chr <- df[, .(mean_chr = mean(get(yvar))), by = .(species_label, chrom_num, chrom_start_Mb, chrom_end_Mb)]
means_genome <- df[, .(mean_genome = mean(get(yvar))), by = species_label]

# genome-wide percentile bands for yvar
q <- quantile(df[[yvar]], probs = percentiles, na.rm = TRUE)
bands <- data.table(ymin = q, ymax = Inf, alpha = band_alpha)

# chromosome tick labels at midpoints
chr_ticks <- chr_len[, .(x = (chrom_start_Mb + chrom_end_Mb)/2, lab = chrom_num)]

df    <- df[!is.na(species_label)]
marks <- marks[!is.na(species_label)]

df$species_label <- factor(
  df$species_label,
  levels = c("C. glandium", "C. nanulus", "C. caryae")
)

sp_colors <- c("C. glandium"="#E69F00",
               "C. nanulus"="#CC79A7",
               "C. caryae"="#7F3C00")
scale_colour_manual(values = sp_colors)

# === plot ===
p <- ggplot() +
   #percentile backgrounds
  geom_rect(data = bands,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "grey60", alpha = bands$alpha, inherit.aes = FALSE) +
  
  # traces per species
  geom_line(data = df,
            aes(x = cum_Mb, y = .data[[yvar]], colour = species_label, group = species_label),
            linewidth = 2.1, alpha = 0.95) +
  
  # per-chromosome means
  geom_segment(data = means_chr,
               aes(x = chrom_start_Mb, xend = chrom_end_Mb,
                   y = mean_chr, yend = mean_chr, colour = species_label),
               linewidth = 0.4, alpha = 0.9) +
  
  # genome-wide means
  geom_hline(data = means_genome,
             aes(yintercept = mean_genome, colour = species_label),
             linetype = 2, linewidth = 0.45, alpha = 0.9, show.legend = FALSE) +
  # chromosome boundaries
  geom_vline(data = chr_len[chrom_start_Mb > 0],
             aes(xintercept = chrom_start_Mb),
             colour = "grey75", linewidth = 0.25) +
 scale_colour_manual(values = sp_colors) +
 scale_x_continuous(breaks = chr_ticks$x, labels = chr_ticks$lab,
                     expand = c(0, 0.002),
                     name = "Chromosomes (concatenated; labels = chromosome number)") +
  labs(y = yvar, colour = "Species",
       title = sprintf("Concatenated trace of %s (chromosomes ordered, percentile bands)", yvar),
       caption = sprintf("Shaded bands: genome-wide %s percentiles = %s",
                         yvar, paste0(100*percentiles, collapse=", "))) +
   theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black",
                                fill = NA,
                                linewidth = 0.5)
     
if (use_log) p <- p + scale_y_continuous(trans = "log1p", labels = comma)
p <- p+ theme(aspect.ratio = 1)
print(p)
