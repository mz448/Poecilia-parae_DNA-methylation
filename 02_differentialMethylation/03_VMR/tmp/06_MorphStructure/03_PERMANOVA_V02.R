#!/usr/bin/env Rscript

# DATE:   20260803
# AUTHOR: MZF
# VERSION:V02
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:   Test whether the effect of sample identity (in a VMR dataset) is 
#         relevant for the methylation structure
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# QUESTIONS:
#  1) Phenotype PERMANOVA: Is there any multivariate difference among females, 
#                          immaculata, parae and yellow?
#
#  2) Sex PERMANOVA: Do females and males have different VMR methylation profiles?
#     IMPORTANT TO CONSIDER: 
#        An unbalanced PERMANOVA is vulnerable to unequal within-group dispersion:
#        a significant result can reflect differences in spread as well as 
#        differences between group centroids. PERMANOVA is generally more robust 
#        to dispersion heterogeneity in balanced designs than in unbalanced ones.
#        REF: 10.1002/9781118445112.stat07841Digital (DOI)
#     INTERPRETATION:
#         A) significant sex PERMANOVA, nonsignificant PERMDISP: stronger evidence
#            for female–male centroid separation;
#         B) significant sex PERMANOVA and significant PERMDISP: the result may 
#            partly reflect different within-sex variability;
#         C) nonsignificant sex PERMANOVA: insufficient evidence for female–male 
#            separation in that VMR set.
#  3) Male-morph PERMANOVA: Do immaculata, parae and yellow males differ in 
#                           their multivariate VMR methylation profiles?
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

library(data.table)
library(vegan)


# Select one condition ---------------------------------------------------------
    condition <- "real_morph"
    # condition <- "shuffled_morph"
    # condition <- "real_sex"
    # condition <- "shuffled_sex"
    # condition <- "shuffled_all"

# Prepare data -----------------------------------------------------------------
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

# Populate Metadata ------------------------------------------------------------

# Sample metadata
samples <- rownames(meth_sample_by_vmr)

metadata <- data.frame(
  sample = samples,
  stringsAsFactors = FALSE
)

metadata$sex <- ifelse(
  grepl("^pparfmem", metadata$sample),
  "female",
  "male"
)

metadata$phenotype <- ifelse(
  grepl("^pparfmem", metadata$sample), "female",
  ifelse(
    grepl("^pparimem", metadata$sample), "immaculata",
    ifelse(
      grepl("^pparpmem", metadata$sample), "parae",
      ifelse(
        grepl("^pparymem", metadata$sample), "yellow",
        NA_character_
      )
    )
  )
)

metadata$male_morph <- ifelse(
  metadata$sex == "male",
  metadata$phenotype,
  NA_character_
)

metadata$sex <- factor(
  metadata$sex,
  levels = c("female", "male")
)

metadata$phenotype <- factor(
  metadata$phenotype,
  levels = c("female", "immaculata", "parae", "yellow")
)

metadata$male_morph <- factor(
  metadata$male_morph,
  levels = c("immaculata", "parae", "yellow")
)

stopifnot(!anyNA(metadata$sex))
stopifnot(!anyNA(metadata$phenotype))

# Save metadata
write.csv(metadata, file = file.path(out_dir,paste0("PERMANOVA.metadata.",condition,".csv")), row.names = TRUE)

# Euclidean distances between samples
d <- dist(meth_sample_by_vmr, method = "euclidean")

# 1) PERMANOVA for phenotype ---------------------------------------------------

set.seed(123)
permanova_phenotype <- adonis2(
  d ~ phenotype,
  data = metadata,
  permutations = 9999
)
permanova_phenotype
p_p <- as.data.frame(permanova_phenotype)

# Save results
write.csv(p_p, file = file.path(out_dir,paste0("PERMANOVA.phen.results.",condition,".csv")), row.names = TRUE)

# 2.1) PERMANOVA for sex -------------------------------------------------------

set.seed(123)

permanova_sex <- adonis2(
  d ~ sex,
  data = metadata,
  permutations = 9999
)
permanova_sex

# save results
p_s <- as.data.frame(permanova_sex)
write.csv(p_s, file = file.path(out_dir,paste0("PERMANOVA.sex.results.",condition,".csv")), row.names = TRUE)

# 2.2) PERMDISP for sex --------------------------
dispersion_sex <- betadisper(
  d,
  metadata$sex,
  type = "median",
  bias.adjust = TRUE # because dispersion estimates are downward biased when group centroids are estimated from small, unequal sample sizes. Recomended by Vegan
)
dispersion_sex

set.seed(123)

permdisp_sex <- permutest(
  dispersion_sex,
  permutations = 9999
)

dispersion_sex$group.distances
permdisp_sex

# Save the group distances
sex_group_distances <- data.frame(
  sex = names(dispersion_sex$group.distances),
  mean_distance_to_median = as.numeric(
    dispersion_sex$group.distances
  )
)

write.csv(sex_group_distances,
  file = file.path(out_dir,paste0("PERMDISP.sex.group_distances.", condition, ".csv")),row.names = FALSE)

#  Save permutation-test table
pd_s <- permdisp_sex$tab
write.csv(pd_s, file = file.path(out_dir,paste0("PERMDISP.sex.results.",condition,".csv")), row.names = TRUE)



# 3.1) PERMANOVA for male-morph ------------------------------------------------

# Subset metadata and methylaiton matrix:
male_index <- metadata$sex == "male"

meth_males <- meth_sample_by_vmr[
  male_index,
  ,
  drop = FALSE
]

metadata_males <- droplevels(
  metadata[male_index, ]
)

d_males <- dist(
  meth_males,
  method = "euclidean"
)

# Run PERMANOVA
set.seed(123)
permanova_male_morph <- adonis2(
  d_males ~ male_morph,
  data = metadata_males,
  permutations = 9999
)
permanova_male_morph 

# save results
p_m <- as.data.frame(permanova_male_morph)
write.csv(p_m, file = file.path(out_dir,paste0("PERMANOVA.morph.results.",condition,".csv")), row.names = TRUE)


# 3.2) PERMDISP for male-morph -------------------
dispersion_male_morph <- betadisper(
  d_males,
  metadata_males$male_morph,
  type = "median"
)

set.seed(123)

permdisp_male_morph <- permutest(
  dispersion_male_morph,
  permutations = 9999
)

dispersion_male_morph$group.distances
permdisp_male_morph


# save results
pd_m <- permdisp_male_morph$tab
write.csv(pd_m, file = file.path(out_dir,paste0("PERMDISP.morph.results.",condition,".csv")), row.names = TRUE)

# Save mean distance to the group median for each male morph
male_morph_group_distances <- data.frame(
  male_morph = names( dispersion_male_morph$group.distances),
  mean_distance_to_median = as.numeric(dispersion_male_morph$group.distances)
  )

write.csv(male_morph_group_distances,file = file.path(out_dir,
    paste0("PERMDISP.morph.group_distances.",condition,".csv")),row.names = FALSE
  )
