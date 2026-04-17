setDT(df)

rss_wide <- dcast(
  unique(df[, .(chrom_num, chrom, win, win_start, win_mid_Mb, cum_Mb, sp_ABC, PBSa_5Mb)]),
  chrom_num + chrom + win + win_start + win_mid_Mb + cum_Mb ~ sp_ABC,
  value.var = "PBSa_5Mb"
)

setnames(rss_wide, c("A", "B", "C"), c("rss_A", "rss_B", "rss_C"))

thr_A <- quantile(rss_wide$rss_A, 0.95, na.rm = TRUE)
thr_B <- quantile(rss_wide$rss_B, 0.95, na.rm = TRUE)
thr_C <- quantile(rss_wide$rss_C, 0.95, na.rm = TRUE)

rss_wide[, high_A := rss_A >= thr_A]
rss_wide[, high_B := rss_B >= thr_B]
rss_wide[, high_C := rss_C >= thr_C]

rss_wide[, shared_AB_not_C := high_A & high_B & !high_C]

rss_wide[, shared_score := pmin(rss_A, rss_B) - rss_C]

thr_score <- quantile(rss_wide$shared_score, 0.95, na.rm = TRUE)

rss_wide[, shared_AB_not_C := (rss_A >= thr_A) &
           (rss_B >= thr_B) &
           (shared_score >= thr_score)]
# test first:
thr_A <- quantile(rss_wide$rss_A, 0.95, na.rm = TRUE)
thr_B <- quantile(rss_wide$rss_B, 0.95, na.rm = TRUE)
thr_C_low <- quantile(rss_wide$rss_C, 0.75, na.rm = TRUE)

rss_wide[, shared_score := pmin(rss_A, rss_B) - rss_C]

rss_wide[, shared_AB_not_C := (rss_A >= thr_A) &
           (rss_B >= thr_B) &
           (rss_C < thr_C_low)]

rss_tr = rss_wide[rss_wide$shared_AB_not_C =="TRUE"]



########## PBS LIKE (failed)
# from counts df <- as.data.table(counts)

library(data.table)
library(slider)
library(ggplot2)
library(scales)

setDT(df)

# keep only columns needed for PBS-like calculation
bc_dt <- unique(df[, .(
  chrom_num, chr_len_Mb, chrom_start_Mb, chrom_end_Mb,
  chrom, win, win_start, win_mid_Mb, cum_Mb,
  sp_ABC, species_label, n_nodes
)])

# wide format: one row per chromosome-window
bc_wide <- dcast(
  bc_dt,
  chrom_num + chr_len_Mb + chrom_start_Mb + chrom_end_Mb +
    chrom + win + win_start + win_mid_Mb + cum_Mb ~ sp_ABC,
  value.var = "n_nodes",
  fill = 0
)

# confirm columns A/B/C exist
stopifnot(all(c("A", "B", "C") %in% names(bc_wide)))

min_total_nodes <- 10

bc_wide[, total_nodes := A + B + C]
bc_wide[total_nodes < min_total_nodes, c("A", "B", "C") := .(NA_real_, NA_real_, NA_real_)]

bc_wide[, `:=`(
  logA = log1p(A),
  logB = log1p(B),
  logC = log1p(C)
)]

bc_wide[, `:=`(
  d_AB = abs(logA - logB),
  d_AC = abs(logA - logC),
  d_BC = abs(logB - logC)
)]

bc_wide[, `:=`(
  BC_A = (d_AB + d_AC - d_BC) / 2,
  BC_B = (d_AB + d_BC - d_AC) / 2,
  BC_C = (d_AC + d_BC - d_AB) / 2
)]

for (v in c("BC_A", "BC_B", "BC_C")) {
  bc_wide[, (v) := pmax(get(v), 0)]
}

setorder(bc_wide, chrom_num, win)

roll_n <- 5   # 5 windows; increase if you want broader smoothing

bc_wide[, BC_A_roll := slide_dbl(BC_A, mean, .before = roll_n - 1, .complete = TRUE), by = chrom]
bc_wide[, BC_B_roll := slide_dbl(BC_B, mean, .before = roll_n - 1, .complete = TRUE), by = chrom]
bc_wide[, BC_C_roll := slide_dbl(BC_C, mean, .before = roll_n - 1, .complete = TRUE), by = chrom]

bc_long <- melt(
  bc_wide,
  id.vars = c("chrom_num", "chr_len_Mb", "chrom_start_Mb", "chrom_end_Mb",
              "chrom", "win", "win_start", "win_mid_Mb", "cum_Mb", "total_nodes"),
  measure.vars = c("BC_A_roll", "BC_B_roll", "BC_C_roll"),
  variable.name = "metric",
  value.name = "branch_contrast"
)

bc_long[, species_label := fifelse(metric == "BC_A_roll", "C. caryae",
                                   fifelse(metric == "BC_B_roll", "C. nanulus",
                                           "C. glandium"))]


bc_long <- bc_long[is.finite(branch_contrast)]

sp_colors <- c(
  "C. glandium" = "#E69F00",
  "C. nanulus"  = "#CC79A7",
  "C. caryae"   = "#7F3C00"
)

plot_p_bc <- bc_long[species_label %in% c("C. caryae", "C. nanulus")]

p_bc <- ggplot(
  plot_p_bc,
  aes(x = cum_Mb, y = branch_contrast, colour = species_label, group = species_label)
) +
  geom_line(linewidth = 1.0, alpha = 0.95) +
  scale_colour_manual(values = sp_colors) +
  labs(
    x = "Chromosomes (concatenated)",
    y = sprintf("Branch contrast (rolling mean, %d windows)", roll_n),
    colour = "Species",
    title = sprintf("Concatenated branch-contrast trace (%d-window rolling mean)", roll_n)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = c(0.84, 0.84),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.85), colour = "black", linewidth = 0.3),
    legend.key = element_rect(fill = "white", colour = NA),
    aspect.ratio = 1
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1)))

#plot_p_bc <- bc_long[species_label %in% c("C. caryae", "C. nanulus")]
print(p_bc)

chr_len <- unique(plot_p_bc[, .(
  chrom_num,
  chrom,
  chr_len_Mb,
  chrom_start_Mb,
  chrom_end_Mb
)])

setorder(chr_len, chrom_num)

# x-axis ticks at chromosome centers
chr_ticks <- chr_len[, .(
  x = (chrom_start_Mb + chrom_end_Mb) / 2,
  lab = as.character(chrom_num)
)]

chr_len <- unique(plot_p_bc[, .(chrom_num, chrom, chr_len_Mb)])
setorder(chr_len, chrom_num)
chr_len[, chrom_start_Mb := c(0, cumsum(head(chr_len_Mb, -1)))]
chr_len[, chrom_end_Mb := chrom_start_Mb + chr_len_Mb]

chr_ticks <- chr_len[, .(
  x = (chrom_start_Mb + chrom_end_Mb) / 2,
  lab = as.character(chrom_num)
)]

percentiles <- c(0.90, 0.95, 0.99)

vals <- plot_p_bc$branch_contrast[is.finite(plot_p_bc$branch_contrast)]
qs <- quantile(vals, probs = percentiles, na.rm = TRUE)

bands <- data.table(
  ymin = c(qs[1], qs[2], qs[3]),
  ymax = c(qs[2], qs[3], max(vals, na.rm = TRUE)),
  alpha = c(0.08, 0.12, 0.18)
)

sp_colors <- c(
  "C. caryae"  = "#7F3C00",
  "C. nanulus" = "#CC79A7"
)

p <- ggplot() +
  geom_rect(
    data = bands,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey60",
    alpha = bands$alpha
  ) +
  # traces
  geom_line(
    data = plot_p_bc,
    aes(x = cum_Mb, y = branch_contrast, colour = species_label, group = species_label),
    linewidth = 0.9,
    alpha = 0.95
  ) +
  
  # per-chromosome means (optional)
  geom_segment(
    data = plot_p_bc[, .(
      mean_chr = mean(branch_contrast, na.rm = TRUE),
      chrom_start_Mb = first(chrom_start_Mb),
      chrom_end_Mb = first(chrom_end_Mb)
    ), by = .(species_label, chrom_num)],
    aes(x = chrom_start_Mb, xend = chrom_end_Mb,
        y = mean_chr, yend = mean_chr, colour = species_label),
    linewidth = 0.45,
    alpha = 0.9
  ) +
  
  # genome-wide means (optional)
  geom_hline(
    data = plot_p_bc[, .(mean_genome = mean(branch_contrast, na.rm = TRUE)), by = species_label],
    aes(yintercept = mean_genome, colour = species_label),
    linetype = 2,
    linewidth = 0.45,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  
  # chromosome boundaries
  geom_vline(
    data = chr_len[chrom_start_Mb > 0],
    aes(xintercept = chrom_start_Mb),
    colour = "grey70",
    linewidth = 0.3
  ) +
  
  scale_colour_manual(values = sp_colors) +
  scale_x_continuous(
    breaks = chr_ticks$x,
    labels = chr_ticks$lab,
    expand = c(0, 0.002),
    name = "Chromosomes (concatenated; labels = chromosome number)"
  ) +
  labs(
    y = "Branch contrast (rolling mean, 5 windows)",
    colour = "Species",
    title = "Concatenated branch-contrast trace (chromosomes ordered)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = c(0.88, 0.84),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.85), colour = "black", linewidth = 0.3),
    legend.key = element_rect(fill = "white", colour = NA),
    aspect.ratio = 1
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1)))

print(p)


### find overlaps:

# example thresholds
thr_rss_A <- quantile(rss_wide$rss_A, 0.95, na.rm = TRUE)
thr_rss_B <- quantile(rss_wide$rss_B, 0.95, na.rm = TRUE)

thr_bc_A  <- quantile(bc_wide$BC_A_roll, 0.95, na.rm = TRUE)
thr_bc_B  <- quantile(bc_wide$BC_B_roll, 0.95, na.rm = TRUE)

# define outliers
rss_wide[, rss_outlier_AB := (rss_A >= thr_rss_A) & (rss_B >= thr_rss_B)]
bc_wide[,  bc_outlier_AB  := (BC_A_roll >= thr_bc_A) & (BC_B_roll >= thr_bc_B)]

# merge by chrom + win
cmp <- merge(
  rss_wide[, .(chrom, win, rss_outlier_AB)],
  bc_wide[,  .(chrom, win, bc_outlier_AB)],
  by = c("chrom", "win"),
  all = TRUE
)

cmp[is.na(rss_outlier_AB), rss_outlier_AB := FALSE]
cmp[is.na(bc_outlier_AB),  bc_outlier_AB  := FALSE]

# overlap table
table(cmp$rss_outlier_AB, cmp$bc_outlier_AB)
