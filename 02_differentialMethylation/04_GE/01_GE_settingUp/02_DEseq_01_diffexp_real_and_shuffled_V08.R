#!/usr/bin/env Rscript
# DATE:       2025-10-11
# AUTHOR:     MZF (revised structure with contrast order control)
# SCRIPT:     02_DEseq_01_diffexp_real_and_shuffled_V08.R
# VERSION:    08
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Stage A — Use methyl.tsv to define (condition, comparison) pairs and group
#   memberships, map methyl IDs → RNA IDs, subset by Tissue/Treatment, run DESeq2
#   (design ~ Group), and export results.
#   
#   V08 — Adds support for explicit contrast order control via orderOfcontrasts.tsv
#         This allows matching the ref/alt direction of existing results by
#         specifying which group should be reference (first) vs alternate (second).
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

# BiocManager::install("apeglm")

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
})


# ------------------------- 0) FIXED PARAMETERS -------------------------------
ts     <- "Muscle"  # Tissue filter
trmnt  <- "veh"     # Treatment filter
coveragethrd <- 0.0 # Min avg reads per sample (0 disables additional raw-count filtering)
Wald_FDR     <- 0.05

# ------------------------- 1) INPUT PATHS ------------------------------------
counts_path      <- "Combined_Counts.tsv"              # genes x samples
sample_info_path <- "sample_info_with_comparisons.tsv" # has Sample, Tissue, Treatment, Morph, ...
methyl_path      <- "methyl.tsv"                       # condition, comparison, group, Rep1, Rep2, Rep3
map_path         <- "rnaseq_to_methyl_map.tsv"         # RNASample, methylationSample
gene_names       <- "gene_names.tsv"                   # not used yet

# NEW: Optional contrast order specification
contrast_order_path <- "orderOfcontrasts.tsv"          # ref, new, first, second

# ------------------------- 2) OUTPUT PATHS -----------------------------------
long_out_path    <- "DESeq2_allComparisons.tsv"
wide_out_path    <- paste0(ts, "_real-and-shuffled_combined_DEGs_CPMfilter.tsv")
objects_root     <- "deseq_objects"

# ------------------------- 3) READ INPUTS ------------------------------------
df <- read.table(counts_path, header = TRUE, row.names = 1, sep = "\t",
                 stringsAsFactors = FALSE, strip.white = TRUE, check.names = FALSE)

dkey_base <- read.table(sample_info_path, header = TRUE, sep = "\t",
                        stringsAsFactors = TRUE, check.names = FALSE)

methyl <- read_tsv(methyl_path, show_col_types = FALSE)

rna2meth <- read_tsv(map_path, show_col_types = FALSE,
                     col_names = c("RNASample","methylationSample"))

# NEW: Load contrast order table (optional)
contrast_order <- NULL
if (file.exists(contrast_order_path)) {
  contrast_order <- read_tsv(contrast_order_path, show_col_types = FALSE) %>%
    mutate(across(everything(), ~str_trim(as.character(.))))
  
  message("✓ Loaded contrast order specifications from: ", contrast_order_path)
  message("  Contrast order will be controlled for comparisons listed in 'new' column")
} else {
  message("ℹ No contrast order file found (", contrast_order_path, ")")
  message("  Using default alphabetical group ordering")
}

# ------------------------- 4) BUILD ASSIGNMENTS FROM METHYL -------------------
methyl_long <- methyl %>%
  pivot_longer(cols = starts_with("Rep"),
               names_to = "rep_label",
               values_to = "methylationSample") %>%
  filter(!is.na(methylationSample), methylationSample != "") %>%
  select(condition, comparison, group, methylationSample)

assignments_rna <- methyl_long %>%
  left_join(rna2meth, by = "methylationSample") %>%
  filter(!is.na(RNASample)) %>%
  transmute(condition, comparison, group, Sample = RNASample)

# ------------------------- 5) CPM FILTER SETUP (GLOBAL) ----------------------
## Tunables
cpm_cut              <- 1      # CPM threshold
min_samples_any      <- 3      # global rule: in >= 3 samples
use_pair_group_rule  <- FALSE  # also require CPM per group within each pair
min_per_group        <- 2      # when use_pair_group_rule=TRUE, require >=2 samples per group

## Build a *global* whitelist for this tissue/treatment universe
samples_tt <- dplyr::filter(dkey_base, Tissue == ts, Treatment == trmnt) %>%
  dplyr::pull(Sample) %>%
  unique()

mat_tt <- as.matrix(df[, samples_tt, drop = FALSE])
keep_global_any <- rowSums(edgeR::cpm(mat_tt) >= cpm_cut) >= min_samples_any

# ------------------------- 6) HELPER FUNCTIONS -------------------------------
collect_pairs <- function(assignments_rna) {
  assignments_rna %>%
    distinct(condition, comparison) %>%
    arrange(condition, comparison)
}

prepare_dkey_for_pair <- function(assignments_rna, dkey_base, condition_id, comparison_id, ts, trmnt) {
  members <- assignments_rna %>%
    filter(condition == condition_id, comparison == comparison_id) %>%
    distinct(Sample, group)
  
  dk <- dkey_base %>%
    inner_join(members, by = "Sample") %>%
    filter(Tissue == ts, Treatment == trmnt) %>%
    mutate(Group = group)
  
  groups <- sort(unique(dk$Group))
  list(dkey_subset = dk, groups = groups)
}

# NEW: Function to determine contrast order
determine_contrast_order <- function(comparison_id, groups, contrast_order_table = NULL) {
  # Default: alphabetical order
  ordered_groups <- sort(groups)
  ref_group <- ordered_groups[1]
  alt_group <- ordered_groups[2]
  method <- "alphabetical (default)"
  
  # Check if we have explicit ordering for this comparison
  if (!is.null(contrast_order_table)) {
    match_row <- contrast_order_table %>%
      filter(new == comparison_id)
    
    if (nrow(match_row) == 1) {
      specified_first <- match_row$first[1]
      specified_second <- match_row$second[1]
      
      # Validate that specified groups exist in our data
      if (specified_first %in% groups && specified_second %in% groups) {
        ref_group <- specified_first
        alt_group <- specified_second
        method <- paste0("explicit (from orderOfcontrasts.tsv: ", 
                         match_row$ref[1], " → ", comparison_id, ")")
        
        message("    ✓ Using explicit contrast order: ", ref_group, " (ref) vs ", alt_group, " (alt)")
      } else {
        warning("    ⚠ Specified groups (", specified_first, ", ", specified_second, 
                ") don't match detected groups (", paste(groups, collapse = ", "), 
                "). Using alphabetical order.")
      }
    } else if (nrow(match_row) > 1) {
      warning("    ⚠ Multiple entries for comparison '", comparison_id, 
              "' in contrast order table. Using alphabetical order.")
    }
  }
  
  list(
    ref = ref_group,
    alt = alt_group,
    method = method
  )
}

# Filter count matrix (separate from DESeq2)
filter_count_matrix_for_pair <- function(df_counts, dkey_subset, keep_global_any, 
                                         use_pair_group_rule = FALSE, 
                                         min_per_group = 2, 
                                         cpm_cut = 1,
                                         coverage_thr = 0.0) {
  # STEP 1: Get samples in the CORRECT ORDER (matching dkey_subset row order)
  samples_ordered <- dkey_subset$Sample
  
  # Verify all samples exist in count matrix
  missing <- setdiff(samples_ordered, colnames(df_counts))
  if (length(missing) > 0) {
    stop("Samples missing from count matrix: ", paste(missing, collapse = ", "))
  }
  
  # STEP 2: Subset and REORDER columns to match dkey_subset
  count_mat <- df_counts[keep_global_any, samples_ordered, drop = FALSE]
  
  # CRITICAL VALIDATION: Ensure column order matches
  if (!identical(colnames(count_mat), samples_ordered)) {
    stop("Column order mismatch after subsetting!")
  }
  
  genes_after_global <- nrow(count_mat)
  
  # STEP 3: Optional per-pair per-group CPM rule
  genes_after_pair <- genes_after_global
  if (use_pair_group_rule) {
    # Group vector is now correctly aligned because count_mat columns match dkey_subset rows
    grp_vec <- dkey_subset$Group
    
    cpm_ok <- edgeR::cpm(count_mat) >= cpm_cut
    per_group_ok <- sapply(unique(grp_vec), function(g) {
      rowSums(cpm_ok[, grp_vec == g, drop = FALSE]) >= min_per_group
    })
    keep_pair <- apply(per_group_ok, 1, all)
    count_mat <- count_mat[keep_pair, , drop = FALSE]
    genes_after_pair <- nrow(count_mat)
  }
  
  # STEP 4: Legacy raw-count mean filter
  genes_after_coverage <- genes_after_pair
  if (coverage_thr > 0) {
    min_sum  <- ncol(count_mat) * coverage_thr
    keep_cov <- rowSums(count_mat) >= min_sum
    count_mat <- count_mat[keep_cov, , drop = FALSE]
    genes_after_coverage <- nrow(count_mat)
  }
  
  # Return filtered matrix + diagnostic info
  list(
    count_mat = count_mat,
    genes_after_global = genes_after_global,
    genes_after_pair = genes_after_pair,
    genes_after_coverage = genes_after_coverage
  )
}

# MODIFIED: Run DESeq2 with explicit ref/alt specification
run_deseq_on_filtered_matrix <- function(count_mat, dkey_subset, ref_group, alt_group, fdr = 0.05) {
  # CRITICAL VALIDATION: Verify alignment before DESeq2
  if (!identical(colnames(count_mat), dkey_subset$Sample)) {
    stop("FATAL: Count matrix columns do NOT match dkey_subset$Sample order!")
  }
  
  # Validate that specified groups exist
  if (!ref_group %in% dkey_subset$Group || !alt_group %in% dkey_subset$Group) {
    stop("FATAL: Specified groups (", ref_group, ", ", alt_group, 
         ") not found in data! Available: ", paste(unique(dkey_subset$Group), collapse = ", "))
  }
  
  # Build DESeqDataSet - NO reordering needed since alignment is guaranteed
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = dkey_subset,  # Already in correct order
    design    = ~ Group
  )
  
  # Set reference level EXPLICITLY as specified
  dds$Group <- factor(dds$Group, levels = c(ref_group, alt_group))
  
  # Run DESeq2
  dds <- DESeq2::DESeq(dds)
  
  # Extract results
  res <- DESeq2::results(dds, alpha = fdr)
  res <- res[order(res$padj), ]
  
  # Shrink log2FC
  coef_name  <- paste0("Group_", alt_group, "_vs_", ref_group)
  res_shrunk <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
  res_shrunk <- res_shrunk[order(res_shrunk$padj), ]
  
  list(dds = dds, res = res, res_shrunk = res_shrunk,
       ref = ref_group, alt = alt_group, coef_name = coef_name)
}

result_to_tbl <- function(res, res_shrunk, condition, comparison, ref, alt, ts, trmnt) {
  tibble(
    gene_id         = rownames(as.data.frame(res)),
    condition       = condition,
    comparison      = comparison,
    group_ref       = ref,
    group_alt       = alt,
    tissue          = ts,
    treatment       = trmnt
  ) %>%
    bind_cols(
      as.data.frame(res) %>%
        select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
        as_tibble()
    ) %>%
    mutate(log2FoldChange_shrunk = res_shrunk$log2FoldChange[match(gene_id, rownames(res_shrunk))])
}

save_dds_and_normcounts <- function(dds, condition, comparison, objects_root = "deseq_objects") {
  out_dir <- file.path(objects_root, condition, comparison)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(dds, file = file.path(out_dir, "dds.rds"))
  
  norm_mat <- counts(dds, normalized = TRUE)
  norm_df  <- as.data.frame(norm_mat) %>%
    rownames_to_column(var = "gene_id")
  write_tsv(norm_df, file.path(out_dir, "normalized_counts.tsv"))
  
  invisible(out_dir)
}

# ------------------------- 7) MAIN LOOP (WITH CONTRAST ORDER CONTROL) --------
pairs_tbl <- collect_pairs(assignments_rna)

all_results <- vector("list", length = nrow(pairs_tbl))
names(all_results) <- paste(pairs_tbl$condition, pairs_tbl$comparison, sep = "::")

for (i in seq_len(nrow(pairs_tbl))) {
  cond_id <- pairs_tbl$condition[i]
  comp_id <- pairs_tbl$comparison[i]
  
  # Prepare metadata for this pair
  prep <- prepare_dkey_for_pair(assignments_rna, dkey_base, cond_id, comp_id, ts, trmnt)
  
  if (length(prep$groups) != 2) {
    message("Skipping ", cond_id, " / ", comp_id, " (requires exactly 2 groups after filters).")
    next
  }
  
  # NEW: Determine contrast order (alphabetical or explicit)
  contrast_spec <- determine_contrast_order(comp_id, prep$groups, contrast_order)
  ref <- contrast_spec$ref
  alt <- contrast_spec$alt
  
  title_hint <- paste(ts, ref, "vs", alt)
  
  cat("====================================================================\n")
  cat("Running DESeq2 for:", cond_id, "/", comp_id, "—", title_hint, "\n")
  cat("Contrast order method:", contrast_spec$method, "\n")
  cat("Samples in pair:", paste(prep$dkey_subset$Sample, collapse = ", "), "\n")
  cat("Sample order (CRITICAL):", paste(prep$dkey_subset$Sample, collapse = " -> "), "\n")
  
  # Group membership summary
  group_summary <- prep$dkey_subset %>%
    group_by(Group) %>%
    summarise(samples = paste(Sample, collapse = ", "), .groups = "drop")
  cat("Group assignments:\n")
  for (j in seq_len(nrow(group_summary))) {
    cat("  ", group_summary$Group[j], ": ", group_summary$samples[j], "\n", sep = "")
  }
  
  # ===== STEP 1: FILTER COUNT MATRIX (separate from DESeq2) =====
  filtered <- filter_count_matrix_for_pair(
    df_counts = df,
    dkey_subset = prep$dkey_subset,
    keep_global_any = keep_global_any,
    use_pair_group_rule = use_pair_group_rule,
    min_per_group = min_per_group,
    cpm_cut = cpm_cut,
    coverage_thr = coveragethrd
  )
  
  # Print filtering summary
  total_genes_universe <- nrow(df)
  cat("---- Gene filtering summary ----\n")
  cat("Universe (all genes):                           ", total_genes_universe, "\n")
  cat("After GLOBAL CPM filter (>=", cpm_cut, "in >=", min_samples_any, "samples): ",
      sum(keep_global_any), "\n", sep = "")
  if (use_pair_group_rule) {
    cat("After PER-PAIR group CPM rule (>=", min_per_group, "/group):  ",
        filtered$genes_after_pair, "\n", sep = "")
  }
  if (coveragethrd > 0) {
    cat("After legacy raw-count coverage filter:         ", filtered$genes_after_coverage, "\n", sep = "")
  }
  cat("Genes entering DESeq2:                          ", nrow(filtered$count_mat), "\n")
  
  # Verify column alignment
  cat("Column alignment check: ", 
      ifelse(identical(colnames(filtered$count_mat), prep$dkey_subset$Sample), 
             "✓ PASS", "✗ FAIL"), "\n\n")
  
  # ===== STEP 2: RUN DESEQ2 WITH EXPLICIT CONTRAST ORDER =====
  fit <- run_deseq_on_filtered_matrix(
    count_mat = filtered$count_mat,
    dkey_subset = prep$dkey_subset,
    ref_group = ref,
    alt_group = alt,
    fdr = Wald_FDR
  )
  
  # Summary stats
  n_significant <- sum(fit$res$padj <= Wald_FDR, na.rm = TRUE)
  cat("Significant DEGs (padj <=", Wald_FDR, "): ", n_significant, "\n", sep = "")
  cat("Contrast direction: ", ref, " (reference) vs ", alt, " (alternate)\n", sep = "")
  cat("  → Positive log2FC = higher in ", alt, "\n", sep = "")
  cat("  → Negative log2FC = higher in ", ref, "\n\n", sep = "")
  
  # Store results
  all_results[[i]] <- result_to_tbl(fit$res, fit$res_shrunk,
                                    condition = cond_id,
                                    comparison = comp_id,
                                    ref = fit$ref, alt = fit$alt,
                                    ts = ts, trmnt = trmnt)
  
  # Save dds + normalized counts
  save_dds_and_normcounts(fit$dds, cond_id, comp_id, objects_root = objects_root)
}

joint_df <- bind_rows(all_results[!vapply(all_results, is.null, logical(1))])

# ------------------------- 8) WRITE LONG OUTPUT ------------------------------
write_tsv(joint_df, long_out_path)

# ------------------------- 9) BUILD & WRITE WIDE OUTPUT ----------------------
wide_df <- joint_df %>%
  select(gene_id, comparison, baseMean, log2FoldChange, stat, pvalue, padj) %>%
  pivot_longer(
    cols = c(baseMean, log2FoldChange, stat, pvalue, padj),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = dplyr::recode(metric, log2FoldChange = "l2fc"),
    comp_metric = paste0(comparison, "_", metric)
  ) %>%
  select(gene_id, comp_metric, value) %>%
  distinct() %>%
  pivot_wider(names_from = comp_metric, values_from = value) %>%
  arrange(gene_id) %>%
  mutate(gene_name = gene_id, .before = 2)

write_tsv(wide_df, wide_out_path)

cat("====================================================================\n")
cat("Pipeline complete!\n")
cat("Long output: ", long_out_path, "\n")
cat("Wide output: ", wide_out_path, "\n")
cat("DESeq objects saved to: ", objects_root, "/\n")
if (!is.null(contrast_order)) {
  cat("Contrast order controlled via: ", contrast_order_path, "\n")
}

