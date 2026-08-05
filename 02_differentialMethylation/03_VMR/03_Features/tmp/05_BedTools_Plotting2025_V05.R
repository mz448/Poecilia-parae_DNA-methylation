# DATE:   20250113
# AUTHOR: MZF
# GOAL #########################################################################
# Plot % of count of genes overlapping with ≥ 1 DMR by feature
# Using BedTools output
#
# LIBRARIES: ####################################################################
library(tidyverse)
library(readr)

# USAGE: #######################################################################
# Import the summary Dataframe
input_dir <-  getwd()
DMR.overlaping <- read.table(paste0(input_dir,"/summary.normalized.perChr.countPerGene.DMRsOvelapping.allGroups.allFeatures.PparFemVer2024.transcriptome.6col.bed"), header = FALSE)

# Label columns
DMR.overlaping <-  rename(DMR.overlaping,
              chromosome = V1, 
              count = V2,
              reference = V3,
              precentage = V4, 
              group = V5,
              feature = V6) 
# Clean the chromosome names
DMR.overlaping$chromosome <- sub("Parae_", "", DMR.overlaping$chromosome) # Remove "Parae_" prefix
# Include only the chromosomes in range [1:23] (this excludes any scaffolds)
DMR.overlaping <- DMR.overlaping[DMR.overlaping$chromosome %in% sprintf("%02d", 1:23), ]


# Add flags in new columns 
DMR.overlaping <- DMR.overlaping %>% mutate(
    chromosome_type = ifelse(chromosome == "12", "Sex_Ch", "Autosome"),
    diversity_axis = ifelse(group %in% c("f_vs_i", "f_vs_p", "f_vs_y"), "sex", "morph")
  )

#Save Summary
write_tsv(DMR.overlaping, "summary.normalized.perChr.countPerGene.DMRsOvelapping.allGroups.allFeatures.PparFemVer2024.transcriptome.8col.bed")

# PLOTING ######################################################################
# Individual plot :WORKING
# DMR.overlaping %>%
#   filter(feature == "CDS") %>%
#   ggplot(aes(x = chromosome, y = precentage, fill = group)) + # Use `fill` instead of `color` for bars
#   geom_col(position = "dodge", alpha = 0.7) +                # Use `geom_col` for bar plot
#   labs(
#     title = "% of genes overlaping on the CDS with ≥ 1 DMR",
#     x = "Chromosome",
#     y = "Percentage",
#     fill = "Group"
#   ) +
#   theme_minimal() + # Minimal theme for cleaner visuals
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
#     plot.title = element_text(hjust = 0.5)            # Center align title
#   )

# # Relevant features in a Group plot :WORKING
# DMR.overlaping %>%
#   filter(feature %in% c("five_prime_UTR", "gene", "three_prime_UTR")) %>% # Filter desired features
#   ggplot(aes(x = chromosome, y = precentage, fill = group)) +
#   geom_col(position = "dodge", alpha = 0.7) +                # Use `geom_col` for bar plot
#   labs(
#     title = "% of Genes Overlapping with ≥ 1 DMR by Feature",
#     x = "Chromosome",
#     y = "% of genes per chromosome",
#     fill = "Group"
#   ) +
#   facet_wrap(~feature, scales = "fixed") +                   # Facet by feature with fixed axes
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
#     plot.title = element_text(hjust = 0.5)            # Center align title
#   )
# ggsave("PercentOfGenes.Overlaping.byFeature.pdf", width = 8, height = 6)
rm (mean_lines)

# Plot relevant features in a sinlge plot. Adding specific colors :WORKING
# set colors for groups
group_colors <- c(
  "f_vs_i" = "blue", 
  "f_vs_p" = "blue", 
  "f_vs_y" = "blue", 
  "i_vs_y" = "red",
  "p_vs_i" = "red",
  "y_vs_p" = "red")

# Calculate means for each feature
mean_lines <- DMR.overlaping %>%
  group_by(feature) %>%
  summarise(mean_percentage = mean(precentage, na.rm = TRUE))

# Calculate medians for each feature
median_lines <- DMR.overlaping %>%
  group_by(feature) %>%
  summarise(median_percentage = median(precentage, na.rm = TRUE))

# Define the desired order for features and chromosomes
desired_features <- c("gene",
                "five_prime_UTR",
                "geneBody",
                "geneBody_exon",
                "geneBody_intron",
                "three_prime_UTR",
                "u2000",
                "d2000",
                "intron",
                "exon",
                "five_prime_UTR_exon",
                "three_prime_UTR_exon")

# filter mean Lines
mean_lines_filtered <- mean_lines %>%
  filter(feature %in% desired_features)

# filter median Lines
# median_lines_filtered <- median_lines %>%
#   filter(feature %in% desired_features)

# 1) Plot bar plot + Means ##########################################################
DMR.overlaping %>%
  filter(feature %in% desired_features) %>% # Filter desired features
  ggplot(aes(x = chromosome, y = precentage, fill = group)) +
  geom_col(position = "dodge", alpha = 0.7) +                # Use `geom_col` for bar plot
  scale_fill_manual(values = group_colors) +                # Manually set colors for groups
  geom_hline(data = mean_lines_filtered, aes(yintercept = mean_percentage), 
             linetype = "dashed", color = "black") +         # Add mean line for each panel
  labs(
    title = "% of Genes Overlapping with ≥ 1 DMR by Feature",
    x = "Chromosome",
    y = "% of genes per chromosome",
    fill = "Group"
  ) +
  #facet_wrap(~feature, scales = "fixed",ncol=3) +       # Facet by feature with fixed axes
  facet_grid(~feature, scales = "fixed") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5 ), # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5)            # Center align title
  )
# Save the plot
ggsave("./PercentOfGenes.Overlaping.byFeature.pdf", width = 30, height = 6)


# 1) Plot bar plot + Medians #####################################################
DMR.overlaping %>%
  filter(feature %in% desired_features) %>% # Filter desired features
  ggplot(aes(x = chromosome, y = precentage, fill = group)) +
  geom_col(position = "dodge", alpha = 0.7) +                # Use `geom_col` for bar plot
  scale_fill_manual(values = group_colors) +                # Manually set colors for groups
  geom_hline(data = median_lines_filtered, aes(yintercept = median_percentage), 
             linetype = "dashed", color = "black") +         # Add mean line for each panel
  labs(
    title = "% of Genes Overlapping with ≥ 1 DMR by Feature",
    x = "Chromosome",
    y = "% of genes per chromosome",
    fill = "Group"
  ) +
  #facet_wrap(~feature, scales = "fixed",ncol=3) +       # Facet by feature with fixed axes
  facet_grid(~feature, scales = "fixed") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5 ), # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5)            # Center align title
  )
# Save the plot
ggsave("./PercentOfGenes.Overlaping.byFeature.medians.pdf", width = 30, height = 6)

