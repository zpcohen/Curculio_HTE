library(data.table)
library(GenomicRanges)
library(IRanges)


nodes_df <- fread("/Users/zacharycohen/Desktop/FAU_ZPC/90pc_nodes_Mar7/Nobreaks_sorted_Ccaryae_All_90pcMar9.bed", header = FALSE)

setnames(nodes_df, c("chrom","start","end"))
getwd()
# BED is 0-based, half-open. Convert to 1-based inclusive for GRanges.
nodes_gr <- GRanges(
  seqnames = nodes_df$chrom,
  ranges   = IRanges(start = nodes_df$start + 1L, end = nodes_df$end),
  strand   = "*"
)
nodes_gr

union_len <- function(starts, ends) {
  o <- order(starts, ends)
  starts <- starts[o]; ends <- ends[o]
  cur_s <- starts[1]; cur_e <- ends[1]
  total <- 0L
  if (length(starts) == 1) return(ends[1] - starts[1] + 1L)
  for (i in 2:length(starts)) {
    s <- starts[i]; e <- ends[i]
    if (s <= cur_e + 1L) {
      if (e > cur_e) cur_e <- e
    } else {
      total <- total + (cur_e - cur_s + 1L)
      cur_s <- s; cur_e <- e
    }
  }
  total + (cur_e - cur_s + 1L)
}

# te_gr @ line ~759
get_te_age_tbl <- function(regions_gr, label, te_gr) {
  hits <- findOverlaps(regions_gr, te_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(data.table())
  
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  seg <- pintersect(regions_gr[q], te_gr[s])
  seg_start <- start(seg)
  seg_end   <- end(seg)
  seg_bp    <- seg_end - seg_start + 1L
  
  family <- as.character(mcols(te_gr)$type)[s]
  kimura <- as.numeric(mcols(te_gr)$KIMURA80)[s]
  
  keep <- !is.na(kimura) & seg_bp > 0
  q <- q[keep]; family <- family[keep]; kimura <- kimura[keep]
  seg_start <- seg_start[keep]; seg_end <- seg_end[keep]; seg_bp <- seg_bp[keep]
  
  if (!length(q)) return(data.table())
  
  key <- paste(q, family, sep="||")
  
  st_list <- base::split(seg_start, key)
  en_list <- base::split(seg_end,   key)
  bp_list <- base::split(seg_bp,    key)
  km_list <- base::split(kimura,    key)
  
  keys <- names(st_list)
  
  te_bp_union <- vapply(seq_along(keys), function(i) union_len(st_list[[i]], en_list[[i]]), integer(1))
  mean_kim    <- vapply(seq_along(keys), function(i) weighted.mean(km_list[[i]], w = bp_list[[i]]), numeric(1))
  
  frac_young  <- vapply(seq_along(keys), function(i) {
    bp <- bp_list[[i]]; km <- km_list[[i]]
    sum(bp[km <= 0.02]) / sum(bp)
  }, numeric(1))
  
  frac_old    <- vapply(seq_along(keys), function(i) {
    bp <- bp_list[[i]]; km <- km_list[[i]]
    sum(bp[km > 0.02]) / sum(bp)
  }, numeric(1))
  
  region_id <- as.integer(sub("\\|\\|.*$", "", keys))
  fam_out   <- sub("^.*\\|\\|", "", keys)
  
  data.table(
    group      = label,
    region_id  = region_id,
    chrom      = as.character(seqnames(regions_gr))[region_id],
    start      = start(regions_gr)[region_id],
    end        = end(regions_gr)[region_id],
    region_len = width(regions_gr)[region_id],
    family     = fam_out,
    te_bp      = te_bp_union,
    te_frac    = te_bp_union / width(regions_gr)[region_id],
    mean_kimura = mean_kim,
    frac_young = frac_young,
    frac_old   = frac_old
  )
}

te_te_age <- data.table(
  group       = "TE",
  region_id   = seq_along(te_gr),
  chrom       = as.character(seqnames(te_gr)),
  start       = start(te_gr),
  end         = end(te_gr),
  region_len  = width(te_gr),
  family      = as.character(mcols(te_gr)$type),
  te_bp       = width(te_gr),
  te_frac     = 1,
  kimura = as.numeric(mcols(te_gr)$KIMURA80),
  young  = as.numeric(as.numeric(mcols(te_gr)$KIMURA80) <= 0.05),
  old    = as.numeric(as.numeric(mcols(te_gr)$KIMURA80) > 0.05)
)

                        family_counts <- te_te_age %>%
  filter(young == 1 | old == 1) %>%
  mutate(age_class = case_when(
    young == 1 ~ "young",
    old   == 1 ~ "old"
  )) %>%
  group_by(family, age_class) %>%
  summarise(
    n_nodes = n(),
    total_te_bp = sum(te_bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = age_class,
    values_from = c(n_nodes, total_te_bp),
    values_fill = list(n_nodes = 0, total_te_bp = 0)
  ) %>%
  mutate(
    young_frac_nodes = n_nodes_young / (n_nodes_young + n_nodes_old),
    young_frac_bp    = total_te_bp_young / (total_te_bp_young + total_te_bp_old)
  ) %>%
  arrange(desc(total_te_bp_young))
family_counts <- node_age_calls %>%
  filter(age_class %in% c("young","old")) %>%
  group_by(family, age_class) %>%
  summarise(
    n_nodes = n(),
    total_te_bp = sum(te_bp),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = age_class,
    values_from = c(n_nodes, total_te_bp),
    values_fill = 0
  ) %>%
  mutate(
    young_frac_nodes = n_nodes_young / (n_nodes_young + n_nodes_old),
    young_frac_bp    = total_te_bp_young / (total_te_bp_young + total_te_bp_old)
  ) %>%
  arrange(desc(total_te_bp_young))
write.table(family_counts, "/Users/zacharycohen/Desktop/FAU_ZPC/EarlGrey_output/K0.05_young_CgTEs.tsv", sep = "\t")
write.table(family_counts, "/Users/zacharycohen/Desktop/FAU_ZPC/EarlGrey_output/K0.05_young_CnTEs.tsv", sep = "\t")
write.table(family_counts, "/Users/zacharycohen/Desktop/FAU_ZPC/EarlGrey_output/K0.05_young_CcTEs.tsv", sep = "\t")
