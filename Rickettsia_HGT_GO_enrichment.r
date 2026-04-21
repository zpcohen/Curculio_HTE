full_c1 <- fread("/Users/zacharycohen/Desktop/FAU_ZPC/C1_hgt.terminology", sep="\t", header=FALSE)
head(full_c1$V2)


fg_long <- full_c1 %>%
  transmute(gene = V1, go_list = V2) %>%
  mutate(go_list = na_if(go_list, "")) %>%
  filter(!is.na(go_list)) %>%
  separate_rows(go_list, sep = ";") %>%
  filter(grepl("^GO:[0-9]+$", go_list)) %>%
  rename(GO = go_list)


go2gene <- split(full_cc_unique$gene, full_cc_unique$GO)

enrich_tbl_c1 <- lapply(names(go2gene), function(go) {
  
  genes_go <- unique(go2gene[[go]])
  
  a <- sum(fg_long %in% genes_go)                 # FG ∩ GO
  b <- length(fg_long) - a                        # FG not GO
  c <- sum(bg_genes %in% genes_go) - a             # BG ∩ GO minus FG
  d <- length(bg_genes) - a - b - c                # BG not GO
  
  if (a == 0) return(NULL)  # no signal
  
  ft <- fisher.test(matrix(c(a,b,c,d), nrow=2))
  
  data.frame(
    GO = go,
    fg_hits = a,
    bg_hits = a + c,
    fg_total = length(fg_long),
    bg_total = length(bg_genes),
    odds_ratio = ft$estimate,
    p_value = ft$p.value
  )
}) %>% bind_rows()

enrich_tbl_c1 <- enrich_tbl_c1 %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    enrichment = (fg_hits / fg_total) / (bg_hits / bg_total)
  ) %>%
  arrange(p_adj, desc(enrichment))
