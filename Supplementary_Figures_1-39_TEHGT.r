### Generate TE & HGT plots 
library(tidyverse)
library(GenomicRanges)

hgt_tbl <- read_tsv("/Users/zacharycohen/Desktop/FAU_ZPC/EarlGrey_output/Ccaryae_ALL_hgt_curated_Rinput.tsv", col_names = FALSE, show_col_types = FALSE) %>%
  setNames(c("chrom","start","end","n_hits","best_id","best_evalue","best_bitscore","size")) %>%
  filter(best_evalue <= 1e-10, best_bitscore >= 100)  # tune for your use case

HGT_gr <- GRanges(
  seqnames = hgt_tbl$chrom,
  ranges   = IRanges(start = hgt_tbl$start, end = hgt_tbl$end)
)

# Optional: merge HGT loci that are very close (improves visibility)
#HGT_gr <- reduce(HGT_gr, min.gapwidth = 50001)  # merge within 50 kb

HGT_gr <- GenomicRanges::reduce(HGT_gr, min.gapwidth =2000001)  # merge within 500 kb
HGT_gr <- resize(HGT_gr, width = width(HGT_gr) + 2*200000, fix="center")  # pad ±200 kb
head(HGT_gr)

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(readr)

# Load repeats GFF and convert to GRanges (intervals only)
gff_path <- "/Users/zacharycohen/Desktop/FAU_ZPC/EarlGrey_output/Ccaryae_UMTE_3426divergence.gff"
gff <- rtracklayer::import(gff_path)
gff_gr <- GRanges(seqnames = seqnames(gff), ranges = ranges(gff))
gff_gr <- gff_gr[width(gff_gr) > 0]

# Load fai and build genome windows
fai <- read_tsv("/Users/zacharycohen/Desktop/FAU_ZPC/EarlGrey_output/Ccaryae_unmasked.reformatted.fa.fai",
                col_names = c("chrom","seqlength","offset","linebases","linewidth"),
                show_col_types = FALSE)

seqlens <- setNames(fai$seqlength, fai$chrom)

# Keep only chromosomes in fai (and drop weird seqlevels in repeats not in assembly)
gff_gr <- keepSeqlevels(gff_gr, names(seqlens), pruning.mode="coarse")

win_size <- 2e6
gr_windows <- tileGenome(seqlens, tilewidth = win_size, cut.last.tile.in.chrom = TRUE)

# TE bp overlap per window
hits <- findOverlaps(gr_windows, gff_gr)
bp_overlap <- width(pintersect(gr_windows[queryHits(hits)], gff_gr[subjectHits(hits)]))

window_te_bp <- tibble(
  win_id = queryHits(hits),
  bp     = bp_overlap) %>%
  group_by(win_id) %>%
  summarise(te_bp_total = sum(bp), .groups = "drop") %>%
  right_join(
    tibble(
      win_id    = seq_along(gr_windows),
      chrom     = as.character(seqnames(gr_windows)),
      win_start = start(gr_windows),
      win_end   = end(gr_windows)
    ),
    by = "win_id"
  ) %>%
  mutate(te_bp_total = replace_na(te_bp_total, 0)) %>%
  dplyr::select(chrom, win_start, win_end, te_bp_total) %>%
  mutate(mid = (win_start + win_end) / 2)

search()
find("select")

# Genome-wide z-score (as you had)
window_TE <- window_te_bp %>%
  mutate(
    genome_mean = mean(te_bp_total, na.rm = TRUE),
    genome_sd   = sd(te_bp_total, na.rm = TRUE),
    genome_z    = ifelse(genome_sd == 0, NA_real_, (te_bp_total - genome_mean) / genome_sd)
  )

library(ggplot2)

plot_chr_te <- function(chr, window_TE, HGT_gr, tail_mb = 1, title_prefix = "TE landscape around HGT") {
  df <- window_TE %>% filter(chrom == chr)
  
  if (nrow(df) == 0) return(NULL)
  
  # optional trim the last tail_mb
  x_max_mb <- max(df$mid, na.rm = TRUE) / 1e6 - tail_mb
  df <- df %>% filter(mid/1e6 <= x_max_mb)
  
  # quantiles per chromosome (matches your figure style)
  q <- quantile(df$genome_z, probs = c(0.025, 0.975), na.rm = TRUE)
  lo <- unname(q[1]); hi <- unname(q[2])
  
  # HGT ranges for this chr
  hgt_chr <- HGT_gr[seqnames(HGT_gr) == chr]
  hgt_df  <- as.data.frame(hgt_chr)
  # as.data.frame gives start/end columns, plus seqnames
  
  ggplot(df, aes(x = mid/1e6, y = genome_z)) +
    geom_line(linewidth = 0.6) +
    {if (nrow(hgt_df) > 0)
      geom_rect(
        data = hgt_df,
        fill = "#D55E00",
        aes(xmin = start/1e6, xmax = end/1e6, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        alpha = 0.18
      )
    } +
    geom_hline(yintercept = c(lo, hi), linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = Inf, y = hi, label = "97.5%", hjust = 1.1, vjust = -0.4) +
    annotate("text", x = Inf, y = lo, label = "2.5%",  hjust = 1.1, vjust =  1.4) +
    theme_bw() +
    labs(
      x = paste0("Position on ", chr, " (Mb)"),
      y = "TE density z-score",
      title = paste(title_prefix, "on", chr)
    )
}


plot_chr_te_Cc <- function(chr, window_TE, HGT_gr, tail_mb = 1,
                        title_prefix = "TE landscape around HGT") {
  
  df <- window_TE %>% filter(chrom == chr)
  if (nrow(df) == 0) return(NULL)
  
  x_max_mb <- max(df$mid, na.rm = TRUE) / 1e6 - tail_mb
  df <- df %>% filter(mid/1e6 <= x_max_mb)
  
  q <- quantile(df$genome_z, probs = c(0.025, 0.975), na.rm = TRUE)
  lo <- unname(q[1]); hi <- unname(q[2])
  
  hgt_chr <- HGT_gr[seqnames(HGT_gr) == chr]
  hgt_df  <- as.data.frame(hgt_chr)
  
  yr <- range(df$genome_z, finite = TRUE)
  track_h <- 0.08 * diff(yr)          # thickness of track
  track_ymax <- yr[2]                 # top of plot
  track_ymin <- yr[2] - track_h
  
  ggplot(df, aes(x = mid/1e6, y = genome_z)) +
    geom_line(linewidth = 0.6) +
    # HGT track (thin band)
    {if (nrow(hgt_df) > 0)
      geom_rect(
        data = hgt_df,
        aes(xmin = start/1e6, xmax = end/1e6, ymin = track_ymin, ymax = track_ymax),
        inherit.aes = FALSE,
        fill = "#D55E00", alpha = 0.85
      )
    } +
    geom_hline(yintercept = c(lo, hi), linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = Inf, y = hi, label = "97.5%", hjust = 1.1, vjust = -0.4) +
    annotate("text", x = Inf, y = lo, label = "2.5%",  hjust = 1.1, vjust =  1.4) +
    theme_bw() +
    labs(
      x = paste0("Position on ", chr, " (Mb)"),
      y = "TE density z-score",
      title = paste(title_prefix, "on", chr)
    )
}

view(plot_chr_te_Cc)

outdir <- "/Users/zacharycohen/Desktop/FAU_ZPC/EarlGrey_output/Ccaryae_TE_HGT_zoomMay7/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

chroms <- sort(unique(window_TE$chrom))

plot_chr_te <- function(chr, window_TE2, HGT_gr, tail_mb = 1,
                        title_prefix = "TE landscape around HGT") {
  
  df <- window_TE2 %>% filter(chrom == chr)
  if (nrow(df) == 0) return(NULL)
  
  x_max_mb <- max(df$mid, na.rm = TRUE) / 1e6 - tail_mb
  df <- df %>% filter(mid/1e6 <= x_max_mb)
  
  # chromosome-specific quantiles (your style)
  q <- quantile(df$chrom_z, probs = c(0.025, 0.975), na.rm = TRUE)
  lo <- unname(q[1]); hi <- unname(q[2])
  
  # HGT for this chr
  hgt_chr <- HGT_gr[seqnames(HGT_gr) == chr]
  hgt_df  <- as.data.frame(hgt_chr)
  #hgt_df <- hgt_df %>%
  #  mutate(
   #   hgt_mid_mb = (start + end) / 2 / 1e6,
   #   hgt_width_mb = pmax((end - start) / 1e6, 1),  # minimum 1 Mb visible width
    #  plot_start_mb = hgt_mid_mb - hgt_width_mb / 2,
    #  plot_end_mb   = hgt_mid_mb + hgt_width_mb / 2
   # )
  yr <- range(df$chrom_z, finite = TRUE)
  track_h <- 0.08 * diff(yr)
  track_ymax <- yr[2]
  track_ymin <- yr[2] - track_h
  #ggplot(df, aes(x = mid/1e6, y = genome_z)) +
    #{if (nrow(hgt_df) > 0)
      #geom_rect(
        #data = hgt_df,
        #aes(xmin = start/1e6, xmax = end/1e6,
       #     ymin = -Inf, ymax = Inf, fill = hgt_class),
      #  inherit.aes = FALSE,
     #   alpha = 0.18
    #  )
    #} +
    #geom_line(linewidth = 0.6) +
  ggplot(df, aes(x = mid/1e6, y = chrom_z)) +
    geom_line(linewidth = 0.6) +
    
    # HGT track (colored by high/low/mid association)
    {if (nrow(hgt_df) > 0)
      geom_rect(
        data = hgt_df,
        aes(xmin = start/1e6, xmax = end/1e6,
            ymin = track_ymin, ymax = track_ymax,
            fill = hgt_class),
        inherit.aes = FALSE,
        alpha = 0.9
      )
    } +
    scale_fill_manual(
      values = c(high_TE = "#D55E00", low_TE = "#0072B2", mid_TE = "grey70"),
      name = "HGT overlaps"
    ) +
    # remove these for Fig4:
    geom_hline(yintercept = c(lo, hi), linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = Inf, y = hi, label = "97.5%", hjust = 1.1, vjust = -0.4) +
    annotate("text", x = Inf, y = lo, label = "2.5%",  hjust = 1.1, vjust =  1.4) +
    
    theme_bw() +
    theme(
      panel.grid = element_blank()
    ) +
    labs(
      x = paste0("Position on ", chr, " (Mb)"),
      y = "TE density z-score",
      title = paste(title_prefix, "on", chr)
    )
}

for (chr in chroms) {
  p <- plot_chr_te(chr, window_TE2, HGT_gr, tail_mb = 1,
                   title_prefix = "TE landscape for C. caryae around HGT events")
  if (is.null(p)) next
  p <- p + theme(aspect.ratio = 1)
  
  ggsave(
    filename = file.path(outdir, paste0("TE_HGT_", chr, ".png")),
    plot = p, width = 5, height = 5, dpi = 300
  )
}
