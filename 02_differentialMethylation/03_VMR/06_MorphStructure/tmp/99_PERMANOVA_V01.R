#!/usr/bin/env Rscript

# DATE:   20260803
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
# 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

library(data.table)
library(vegan)


# Select one condition
    # condition <- "real_morph"
    # condition <- "shuffled_morph"
    
    
    # condition <- "real_sex"
    condition <- "shuffled_sex"

out_dir <- paste0("./dataset_Structure/",condition,"_vmr_methylation")
dir.create(out_dir, recursive = TRUE)


# Read VMR table
dt <- fread(paste0("./dataset_Structure/",condition,"_vmr_methylation/vmr.",condition,".morph_and_sample_observed.tsv"))

# Individual-sample methylation columns
meth_cols <- grep("^meth_ppar", names(dt), value = TRUE)

# Optionally retain only VMRs measured in all 12 samples
dt <- dt[n_meth_nonNA_samples == 12]

# Convert to numeric matrix:
# current orientation = VMR rows × sample columns
meth_vmr_by_sample <- as.matrix(dt[, ..meth_cols])
storage.mode(meth_vmr_by_sample) <- "numeric"

# Transpose:
# required orientation = sample rows × VMR columns
meth_sample_by_vmr <- t(meth_vmr_by_sample)

# Clean sample names
rownames(meth_sample_by_vmr) <- sub(
  "^meth_",
  "",
  rownames(meth_sample_by_vmr)
)

# Sample metadata
samples <- rownames(meth_sample_by_vmr)


metadata <- data.frame(
  sample = samples,
  morph = factor(
    c(
      rep("female",      3),
      rep("immaculata",  3),
      rep("parae",       3),
      rep("yellow",      3)
    ),
    levels = c("female", "immaculata", "parae", "yellow")
  )
)

# metadata <- data.frame(
#   sample = samples,
#   stringsAsFactors = FALSE
# )
# 
# metadata$sex <- ifelse(
#   grepl("^pparfmem", metadata$sample),
#   "female",
#   "male"
# )
# 
# metadata$phenotype <- ifelse(
#   grepl("^pparfmem", metadata$sample), "female",
#   ifelse(
#     grepl("^pparimem", metadata$sample), "immaculata",
#     ifelse(
#       grepl("^pparpmem", metadata$sample), "parae",
#       ifelse(
#         grepl("^pparymem", metadata$sample), "yellow",
#         NA_character_
#       )
#     )
#   )
# )
# 
# metadata$male_morph <- ifelse(
#   metadata$sex == "male",
#   metadata$phenotype,
#   NA_character_
# )
# 
# metadata$sex <- factor(
#   metadata$sex,
#   levels = c("female", "male")
# )
# 
# metadata$phenotype <- factor(
#   metadata$phenotype,
#   levels = c("female", "immaculata", "parae", "yellow")
# )
# 
# metadata$male_morph <- factor(
#   metadata$male_morph,
#   levels = c("immaculata", "parae", "yellow")
# )
# 
# stopifnot(!anyNA(metadata$sex))
# stopifnot(!anyNA(metadata$phenotype))

# Euclidean distances between samples
d <- dist(meth_sample_by_vmr, method = "euclidean")

# Standard PERMANOVA
set.seed(123)
permanova <- adonis2(
  d ~ morph,
  data = metadata,
  permutations = 9999
)

permanova

# Convert to data frame
res_df <- as.data.frame(permanova)


# Save results

write.csv(res_df, file = file.path(out_dir,paste0("PERMANOVA.results.",condition,".csv")), row.names = TRUE)
write.csv(metadata, file = file.path(out_dir,paste0("PERMANOVA.metadata.",condition,".csv")), row.names = TRUE)
