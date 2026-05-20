# =============================================================================
# Mouth-Tumor-Gut (MTG) Axis: Transmission Index Computation (Zhang Cohort)
# Purpose : Quantify ASV-level co-occurrence across four matched sample types
#           (Mouth → Tumor → Paratumor → Feces) for each individual in the
#           Zhang cohort, and decompose each transmission pathway into
#           oral-derived vs. non-oral-derived fractions.
#
# Transmission indices computed per individual (GC and CRC separately):
#   MTT  — Mouth-to-Tumor           : shared ASV fecal count in tumor
#   MTP  — Mouth-to-Paratumor       : shared ASV count in paratumor
#   TTG  — Tumor-to-Gut (feces)     : shared ASV count in feces
#   PTG  — Paratumor-to-Gut (feces) : shared ASV count in feces
#   MTG  — Mouth-to-Gut (feces)     : shared ASV count in feces
#
# Each index is further decomposed:
#   *_not_mouth / *_from_mouth  — tissue-to-feces bacteria absent/present in mouth
#   *_not_feces / *_to_feces    — mouth-to-tissue bacteria absent/present in feces
#   MTG_not_TP / MTG_via_tissue — direct vs. tissue-mediated oral transmission to gut
#
# Inclusion criterion: individuals with all four sample types available (Sum == 4)
#
# INPUT  : data/zhang.15000.RDS — rarefied Zhang phyloseq object (GC + CRC)
# OUTPUT : data/gc.trans.index.csv  — per-subject transmission indices (GC)
#          data/crc.trans.index.csv — per-subject transmission indices (CRC)
#          fig4_GC_transmission.svg / fig4_CRC_transmission.svg
# =============================================================================

library(phyloseq)
library(tidyverse)
library(rlang)      # !!sym() for dynamic column name reference
library(ggalluvial)
library(ggpubr)

zhang <- readRDS("data/zhang.15000.RDS")

# Color scheme for transmission pathway fractions
trans_fill_colors <- c(
  "Absent in feces"     = "gray80",
  "Tissues-to-feces"    = "#CC6600",
  "Absent in mouth"     = "gray80",
  "Oral-derived"        = "#C1E3FF",
  "Absent in tissues"   = "gray80",
  "Via tissues"         = "#B266B2"
)

# =============================================================================
# Helper function: compute all transmission indices for one cancer group
# =============================================================================
compute_transmission_index <- function(ps_cancer, cancer_label,
                                       n_shared_samples) {
  # ── Step 1. Subset and prune each sample type ──────────────────────────────
  ps_mouth <- subset_samples(ps_cancer, Type == "Mouth") %>%
    prune_taxa(rowSums(otu_table(.)) > 0, .)
  ps_tumor <- subset_samples(ps_cancer, Type == "Tumor") %>%
    prune_taxa(rowSums(otu_table(.)) > 0, .)
  ps_para  <- subset_samples(ps_cancer, Type == "Para")  %>%
    prune_taxa(rowSums(otu_table(.)) > 0, .)
  ps_feces <- subset_samples(ps_cancer, Type == "Feces") %>%
    prune_taxa(rowSums(otu_table(.)) > 0, .)

  cat(cancer_label, "sample counts — Mouth:", nsamples(ps_mouth),
      "Tumor:", nsamples(ps_tumor),
      "Para:",  nsamples(ps_para),
      "Feces:", nsamples(ps_feces), "\n")

  # ── Step 2. Extract ASV tables and rename columns to SampleID2 ─────────────
  asv_mouth <- otu_table(ps_mouth); colnames(asv_mouth) <- sample_data(ps_mouth)$SampleID2
  asv_tumor <- otu_table(ps_tumor); colnames(asv_tumor) <- sample_data(ps_tumor)$SampleID2
  asv_para  <- otu_table(ps_para);  colnames(asv_para)  <- sample_data(ps_para)$SampleID2
  asv_feces <- otu_table(ps_feces); colnames(asv_feces) <- sample_data(ps_feces)$SampleID2

  # ── Step 3. Identify individuals with all four sample types available ──────
  sample_summary <- sample_data(ps_cancer) %>%
    data.frame() %>%
    rownames_to_column("SampleID") %>%
    mutate(Value = TRUE) %>%
    pivot_wider(id_cols     = SampleID2,
                names_from  = Type,
                values_from = Value,
                values_fill = FALSE) %>%
    mutate(Sum = Para + Mouth + Feces + Tumor)

  shared_ids <- sample_summary %>%
    filter(Sum == 4) %>%
    pull(SampleID2)

  cat("  Individuals with all 4 sample types (n =", length(shared_ids), ")\n")

  # Restrict ASV tables to fully matched individuals
  asv_mouth <- asv_mouth[, shared_ids]; asv_mouth <- asv_mouth[rowSums(asv_mouth) > 0, ]
  asv_tumor <- asv_tumor[, shared_ids]; asv_tumor <- asv_tumor[rowSums(asv_tumor) > 0, ]
  asv_para  <- asv_para[,  shared_ids]; asv_para  <- asv_para[rowSums(asv_para)   > 0, ]
  asv_feces <- asv_feces[, shared_ids]; asv_feces <- asv_feces[rowSums(asv_feces) > 0, ]

  # ── Step 4. Per-individual transmission index computation ──────────────────
  # Column index reference in tmp_table:
  #   2 = mouth,  3 = tumor,  4 = para,  5 = feces
  trans_df <- data.frame(
    Type           = cancer_label,
    ID             = shared_ids,
    MTT            = 0, MTP            = 0, TTG = 0, PTG = 0, MTG = 0,
    TTG_not_mouth  = 0, PTG_not_mouth  = 0,
    MTG_not_TP     = 0,
    MTT_not_feces  = 0, MTP_not_feces  = 0,
    MTT_to_feces   = 0, MTP_to_feces   = 0,
    TTG_from_mouth = 0, PTG_from_mouth = 0,
    MTG_via_tissue = 0
  )

  for (i in seq_along(shared_ids)) {
    tmp_id <- shared_ids[i]

    # Wide table: rows = ASVs, columns = mouth / tumor / para / feces
    tmp_table <- rbind(
      asv_mouth[, tmp_id] %>% data.frame() %>% rownames_to_column("ASV") %>% mutate(Type = "mouth"),
      asv_tumor[, tmp_id] %>% data.frame() %>% rownames_to_column("ASV") %>% mutate(Type = "tumor"),
      asv_para[,  tmp_id] %>% data.frame() %>% rownames_to_column("ASV") %>% mutate(Type = "para"),
      asv_feces[, tmp_id] %>% data.frame() %>% rownames_to_column("ASV") %>% mutate(Type = "feces")
    ) %>%
      pivot_wider(id_cols = ASV, names_from = Type, values_from = !!sym(tmp_id))
    tmp_table[is.na(tmp_table)] <- 0
    # Column indices: mouth=2, tumor=3, para=4, feces=5

    # 1. MTT: ASVs shared in mouth AND tumor → tumor count
    trans_df$MTT[i] <- tmp_table[rowSums(tmp_table[, 2:3] > 0) == 2, 3] %>% sum(na.rm = TRUE)

    # 2. MTP: ASVs shared in mouth AND paratumor → paratumor count
    trans_df$MTP[i] <- tmp_table[rowSums(tmp_table[, c(2,4)] > 0) == 2, 4] %>% sum(na.rm = TRUE)

    # 3. TTG: ASVs shared in tumor AND feces → fecal count
    trans_df$TTG[i] <- tmp_table[rowSums(tmp_table[, c(3,5)] > 0) == 2, 5] %>% sum(na.rm = TRUE)

    # 4. PTG: ASVs shared in paratumor AND feces → fecal count
    trans_df$PTG[i] <- tmp_table[rowSums(tmp_table[, c(4,5)] > 0) == 2, 5] %>% sum(na.rm = TRUE)

    # 5. MTG: ASVs shared in mouth AND feces → fecal count (direct oral-fecal)
    trans_df$MTG[i] <- tmp_table[rowSums(tmp_table[, c(2,5)] > 0) == 2, 5] %>% sum(na.rm = TRUE)

    # Subset: ASVs absent in the mouth (oral count = 0)
    mouth0 <- tmp_table[tmp_table[, 2] == 0, ]

    # 6. TTG_not_mouth: fecal bacteria from tumor that were absent in mouth
    trans_df$TTG_not_mouth[i] <- mouth0[rowSums(mouth0[, c(3,5)] > 0) == 2, 5] %>% sum(na.rm = TRUE)

    # 7. PTG_not_mouth: fecal bacteria from paratumor absent in mouth
    trans_df$PTG_not_mouth[i] <- mouth0[rowSums(mouth0[, c(4,5)] > 0) == 2, 5] %>% sum(na.rm = TRUE)

    # Subset: ASVs absent in both tumor AND paratumor
    tp0 <- tmp_table[tmp_table[, 3] + tmp_table[, 4] == 0, ]

    # 8. MTG_not_TP: oral bacteria in feces NOT detected in any tissue (direct route)
    trans_df$MTG_not_TP[i] <- tp0[rowSums(tp0[, c(2,5)] > 0) == 2, 5] %>% sum(na.rm = TRUE)

    # Subset: ASVs absent in feces
    feces0 <- tmp_table[tmp_table[, 5] == 0, ]

    # 9. MTT_not_feces: oral bacteria in tumor that did NOT reach feces
    trans_df$MTT_not_feces[i] <- feces0[rowSums(feces0[, c(2,3)] > 0) == 2, 3] %>% sum(na.rm = TRUE)

    # 10. MTP_not_feces: oral bacteria in paratumor that did NOT reach feces
    trans_df$MTP_not_feces[i] <- feces0[rowSums(feces0[, c(2,4)] > 0) == 2, 4] %>% sum(na.rm = TRUE)

    # 11. MTT_to_feces: oral bacteria in tumor that ALSO appeared in feces
    trans_df$MTT_to_feces[i] <- tmp_table[rowSums(tmp_table[, c(2,3,5)] > 0) == 3, 3] %>% sum(na.rm = TRUE)

    # 12. MTP_to_feces: oral bacteria in paratumor that ALSO appeared in feces
    trans_df$MTP_to_feces[i] <- tmp_table[rowSums(tmp_table[, c(2,4,5)] > 0) == 3, 4] %>% sum(na.rm = TRUE)

    # 13. TTG_from_mouth: fecal bacteria from tumor that ALSO appeared in mouth
    trans_df$TTG_from_mouth[i] <- tmp_table[rowSums(tmp_table[, c(2,3,5)] > 0) == 3, 5] %>% sum(na.rm = TRUE)

    # 14. PTG_from_mouth: fecal bacteria from paratumor that ALSO appeared in mouth
    trans_df$PTG_from_mouth[i] <- tmp_table[rowSums(tmp_table[, c(2,4,5)] > 0) == 3, 5] %>% sum(na.rm = TRUE)

    # 15. MTG_via_tissue: oral-fecal bacteria that were ALSO detected in at least one tissue
    tp1 <- tmp_table[tmp_table[, 3] + tmp_table[, 4] > 0, ]
    trans_df$MTG_via_tissue[i] <- tp1[rowSums(tp1[, c(2,5)] > 0) == 2, 5] %>% sum(na.rm = TRUE)
  }

  # ── Step 5. Reshape to long format for visualization ───────────────────────
  trans_long <- trans_df %>%
    pivot_longer(
      cols      = c("MTT_not_feces", "MTT_to_feces",
                    "MTP_not_feces", "MTP_to_feces",
                    "TTG_not_mouth", "TTG_from_mouth",
                    "PTG_not_mouth", "PTG_from_mouth",
                    "MTG_not_TP",    "MTG_via_tissue"),
      names_to  = "Trans",
      values_to = "Count"
    ) %>%
    mutate(
      Group = case_when(
        grepl("TTG", Trans) ~ "TTG",
        grepl("PTG", Trans) ~ "PTG",
        grepl("MTG", Trans) ~ "MTG",
        grepl("MTT", Trans) ~ "MTT",
        grepl("MTP", Trans) ~ "MTP"
      ) %>% factor(levels = c("MTT", "MTP", "TTG", "PTG", "MTG")),

      Trans2 = case_when(
        grepl("from",    Trans) ~ "Oral-derived",
        grepl("not_mouth", Trans) ~ "Absent in mouth",
        grepl("not_feces", Trans) ~ "Absent in feces",
        grepl("to_feces",  Trans) ~ "Tissues-to-feces",
        grepl("via",       Trans) ~ "Via tissues",
        grepl("not_TP",    Trans) ~ "Absent in tissues"
      ) %>% factor(levels = c("Absent in feces", "Absent in mouth",
                               "Absent in tissues", "Oral-derived",
                               "Tissues-to-feces",  "Via tissues"))
    )

  # ── Step 6. Summarize for stacked bar plot (relative abundance %) ──────────
  trans_summary <- trans_long %>%
    group_by(Trans) %>%
    summarise(Rel.abund = sum(Count) / (n_shared_samples * 15000) * 100,
              .groups = "drop") %>%
    mutate(
      Group = case_when(
        grepl("TTG", Trans) ~ "TTG",
        grepl("PTG", Trans) ~ "PTG",
        grepl("MTG", Trans) ~ "MTG",
        grepl("MTT", Trans) ~ "MTT",
        grepl("MTP", Trans) ~ "MTP"
      ) %>% factor(levels = c("MTT", "MTP", "TTG", "PTG", "MTG")),

      Trans2 = case_when(
        grepl("from",      Trans) ~ "Oral-derived",
        grepl("not_mouth", Trans) ~ "Absent in mouth",
        grepl("not_feces", Trans) ~ "Absent in feces",
        grepl("to_feces",  Trans) ~ "Tissues-to-feces",
        grepl("via",       Trans) ~ "Via tissues",
        grepl("not_TP",    Trans) ~ "Absent in tissues"
      ) %>% factor(levels = c("Absent in feces", "Absent in mouth",
                               "Absent in tissues", "Oral-derived",
                               "Tissues-to-feces",  "Via tissues"))
    )

  # ── Step 7. Report pathway fractions (values cited in paper) ───────────────
  s <- trans_long %>% group_by(Trans) %>% summarise(Total = sum(Count), .groups = "drop")
  get_frac <- function(num_key, denom_keys) {
    num   <- s$Total[s$Trans == num_key]
    denom <- sum(s$Total[s$Trans %in% denom_keys])
    sprintf("  %s / (%s): %.1f%%", num_key,
            paste(denom_keys, collapse = " + "),
            num / denom * 100)
  }
  cat("\n--- Pathway fraction summary:", cancer_label, "---\n")
  cat(get_frac("MTT_to_feces",  c("MTT_not_feces",  "MTT_to_feces")),  "\n")
  cat(get_frac("MTP_to_feces",  c("MTP_not_feces",  "MTP_to_feces")),  "\n")
  cat(get_frac("TTG_from_mouth",c("TTG_from_mouth", "TTG_not_mouth")), "\n")
  cat(get_frac("PTG_from_mouth",c("PTG_from_mouth", "PTG_not_mouth")), "\n")
  cat(get_frac("MTG_via_tissue",c("MTG_not_TP",     "MTG_via_tissue")),"\n")

  # ── Step 8. Stacked bar plot ───────────────────────────────────────────────
  p <- trans_summary %>%
    ggplot(aes(x = Group, y = Rel.abund)) +
    geom_bar(aes(fill = Trans2), color = "white", linewidth = 1.5,
             stat = "identity") +
    scale_fill_manual(values = trans_fill_colors) +
    theme_pubr() +
    theme(axis.title.x  = element_blank(),
          legend.position = "right",
          legend.title    = element_blank()) +
    labs(y = "Rel. abund. of transmissible bacteria (%)")

  return(list(trans_df = trans_df, trans_long = trans_long,
              trans_summary = trans_summary, plot = p))
}

# =============================================================================
# Run for GC (Figure 4A) and CRC (Figure 4B)
# n_shared_samples: number of fully matched individuals (all 4 types present)
#   GC: n = 33 (verified from sample_data inspection)
#   CRC: n = 23
# =============================================================================
gc_ps  <- subset_samples(zhang, Group == "GC")
crc_ps <- subset_samples(zhang, Group == "CRC")

res_gc  <- compute_transmission_index(gc_ps,  "GC",  n_shared_samples = 33)
res_crc <- compute_transmission_index(crc_ps, "CRC", n_shared_samples = 23)

# Save per-subject transmission index tables
write.csv(res_gc$trans_df,  "data/gc.trans.index.csv",  row.names = FALSE)
write.csv(res_crc$trans_df, "data/crc.trans.index.csv", row.names = FALSE)
cat("Saved: data/gc.trans.index.csv, data/crc.trans.index.csv\n")

# Save figures (SVG for vector editing in Illustrator)
ggsave("fig4A_GC_transmission.svg",  res_gc$plot,  device = "svg", width = 6, height = 5)
ggsave("fig4B_CRC_transmission.svg", res_crc$plot, device = "svg", width = 6, height = 5)
cat("Saved: fig4A_GC_transmission.svg, fig4B_CRC_transmission.svg\n")
