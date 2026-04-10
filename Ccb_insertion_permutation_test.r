library(dplyr)
library(data.table)

set.seed(1)

chr_lengths <- c(
  CM117440.1 = 299321432,   # C1
  CM117442.1 = 261115320    # C3
)

obs <- data.table(
  insert    = c("C1", "C3"),
  chrom     = c("CM117440.1", "CM117442.1"),
  win_start = c(134001568, 136515532),   # real coordinates here
  win_end   = c(134533379, 137385348)    # real coordinates here
)

flank_bp <- 2000000   # example; replace with your actual flank size
nperm <- 10000

sample_interval <- function(chr_len, insert_len, flank_bp) {
  min_start <- flank_bp + 1
  max_start <- chr_len - flank_bp - insert_len + 1
  start <- sample.int(max_start - min_start + 1, 1) + min_start - 1
  end <- start + insert_len - 1
  c(start = start, end = end)
}

get_region_stat <- function(windows, chr, region_start, region_end, flank_bp) {
  core <- windows %>%
    filter(
      .data$chrom == chr,
      .data$win_end >= region_start,
      .data$win_start <= region_end
    )
  
  left_flank <- windows %>%
    filter(
      .data$chrom == chr,
      .data$win_end >= (region_start - flank_bp),
      .data$win_start <= (region_start - 1)
    )
  
  right_flank <- windows %>%
    filter(
      .data$chrom == chr,
      .data$win_end >= (region_end + 1),
      .data$win_start <= (region_end + flank_bp)
    )
  
  tibble(
    core_mean  = mean(core$chrom_z, na.rm = TRUE),
    flank_mean = mean(c(left_flank$chrom_z, right_flank$chrom_z), na.rm = TRUE)
  )
}

obs_stats <- lapply(seq_len(nrow(obs)), function(i) {
  x <- obs[i]
  out <- get_region_stat(window_TE2, x$chrom, x$win_start, x$win_end, flank_bp)
  cbind(x, out)
}) %>% bind_rows()

obs_stats[, contrast := flank_mean - core_mean]

obs_score <- mean(obs_stats$flank_mean - obs_stats$core_mean, na.rm = TRUE)
obs_score

#permuation test:
set.seed(1)
nperm <- 10000
null_scores <- numeric(nperm)
names(chr_lengths)
setdiff(obs$chrom, names(chr_lengths))

obs[, len := win_end - win_start + 1]
obs$len

for (i in seq_len(nrow(obs))) {
  x <- obs[i]
  print(x$chrom)
  print(chr_lengths[[x$chrom]])
  print(x$len)
}
obs$chrom

for (p in seq_len(nperm)) {
  perm_contrasts <- numeric(nrow(obs))
  
  for (i in seq_len(nrow(obs))) {
    x <- obs[i]
    
    chr_i <- as.character(x$chrom)
    chr_len_i <- unname(chr_lengths[chr_i])
    len_i <- as.numeric(x$len)
    
    pos <- sample_interval(
      chr_len = chr_len_i,
      insert_len = len_i,
      flank_bp = flank_bp
    )
    
    perm_stat <- get_region_stat(
      window_TE2,
      chr = chr_i,
      region_start = pos["start"],
      region_end = pos["end"],
      flank_bp = flank_bp
    )
    
    perm_contrasts[i] <- perm_stat$flank_mean - perm_stat$core_mean
  }
  
  null_scores[p] <- mean(perm_contrasts, na.rm = TRUE)
}

p_emp <- (1 + sum(null_scores >= obs_score, na.rm = TRUE)) / (nperm + 1)
p_emp
