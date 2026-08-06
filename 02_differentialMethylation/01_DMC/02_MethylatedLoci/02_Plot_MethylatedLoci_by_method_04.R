# DATE:   20260806
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:   Plotting %mCpG per chromosomes or chromosome type form 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# INPUT:  1) =Normalized CpG Methylation percentage files for each sample
# OUTPUT: 2 column files with normalized (%) values of CpG methylation
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Load necessary libraries
library(tidyverse)
library(ggplot2)

plots <- "../plots/methylatedLoci/"
dir.create(plots, recursive = TRUE)

# Load normalized counts per chromosome
work_dir <- ("./methylatedLoci/")
f1  <- read.table(paste0(work_dir,"norm.postProb_meth.fmem001.sig_0.9.tsv"), header = FALSE)
f2  <- read.table(paste0(work_dir,"norm.postProb_meth.fmem002.sig_0.9.tsv"), header = FALSE)
f3  <- read.table(paste0(work_dir,"norm.postProb_meth.fmem003.sig_0.9.tsv"), header = FALSE)
p4  <- read.table(paste0(work_dir,"norm.postProb_meth.pmem004.sig_0.9.tsv"), header = FALSE)
p5  <- read.table(paste0(work_dir,"norm.postProb_meth.pmem005.sig_0.9.tsv"), header = FALSE)
p6  <- read.table(paste0(work_dir,"norm.postProb_meth.pmem006.sig_0.9.tsv"), header = FALSE)
y7  <- read.table(paste0(work_dir,"norm.postProb_meth.ymem007.sig_0.9.tsv"), header = FALSE)
y8  <- read.table(paste0(work_dir,"norm.postProb_meth.ymem008.sig_0.9.tsv"), header = FALSE)
y9  <- read.table(paste0(work_dir,"norm.postProb_meth.ymem009.sig_0.9.tsv"), header = FALSE)
i10 <- read.table(paste0(work_dir,"norm.postProb_meth.imem010.sig_0.9.tsv"), header = FALSE)
i11 <- read.table(paste0(work_dir,"norm.postProb_meth.imem011.sig_0.9.tsv"), header = FALSE)
i12 <- read.table(paste0(work_dir,"norm.postProb_meth.imem012.sig_0.9.tsv"), header = FALSE)

df<- bind_rows (f1  %>% mutate(sample="fmem001",sex = "female",phenotype = "Female",    male_morph= "NA"),
                f2  %>% mutate(sample="fmem002",sex = "female",phenotype = "Female",    male_morph= "NA"),
                f3  %>% mutate(sample="fmem003",sex = "female",phenotype = "Female",    male_morph= "NA"),
                p4  %>% mutate(sample="pmem004",sex = "male"  ,phenotype = "Parae",     male_morph= "parae"),
                p5  %>% mutate(sample="pmem005",sex = "male"  ,phenotype = "Parae",     male_morph= "parae"),
                p6  %>% mutate(sample="pmem006",sex = "male"  ,phenotype = "Parae",     male_morph= "parae"),
                y7  %>% mutate(sample="ymem007",sex = "male"  ,phenotype = "Yellow",    male_morph= "yellow"),
                y8  %>% mutate(sample="ymem008",sex = "male"  ,phenotype = "Yellow",    male_morph= "yellow"),
                y9  %>% mutate(sample="ymem009",sex = "male"  ,phenotype = "Yellow",    male_morph= "yellow"),
                i10 %>% mutate(sample="imem010",sex = "male"  ,phenotype = "Immaculata",male_morph= "immaculata"),
                i11 %>% mutate(sample="imem011",sex = "male"  ,phenotype = "Immaculata",male_morph= "immaculata"),
                i12 %>% mutate(sample="imem012",sex = "male"  ,phenotype = "Immaculata",male_morph= "immaculata")
                )

df <-  rename(df,
    chromosome = V1, 
    CpG_proportion = V2)


df <- df %>%
  mutate(
    CpG_pct = CpG_proportion * 100,
    chromosome_type = case_when(
      str_starts(chromosome, "Parae") ~ "Genomic",
      chromosome == "LambdaNEB" ~ "(- control)",
      chromosome == "P_parae_Mitochondria" ~ "Mitochondria",
      chromosome == "pUC19" ~ "(+ control)",
      TRUE ~ NA_character_
    )
  )

# df$chromosome <- sub("Parae_", "", df$chromosome) # Remove "Parae_" prefix

colors <- c(
            Female = "#000000",
            Yellow = "#ffd800",
            Immaculata = "#888888",
            Parae = "#804ca5"
            )

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Plot the bar plot Per chromosome ---------------------------------------------
filename="normalized_counts_DML_PerChromosome_barplot.pdf"

bar_plot <- ggplot(df, aes(x = chromosome, y = CpG_pct, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Normalized Counts of DML per Chromosome ~DNMTools~",
    x = "Chromosome",
    y = "Normalized methylation (mCpG/CpG)") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(bar_plot)
# Save the plot
  ggsave(paste0(plots, "/", filename), plot = bar_plot, width = 8, height = 6)

# Print success message
cat("Bar plot saved to", plots, "\n")


# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Plot mean percentage by chromosome type and group ----------------------------
# Genomic = mean across Parae chromosomes with SEM, separately for each group

colors <- c(
  Female = "#399957",
  Yellow = "#ffd800",
  Immaculata = "#888888",
  Parae = "#804ca5"
)

chromosome_type_summary <- df %>%
  mutate(
    chromosome_type = factor(
      chromosome_type,
      # levels = c("(- control)", "Mitochondria", "Genomic", "(+ control)")
      levels = c("Genomic", "Mitochondria", "(- control)", "(+ control)")
    ),
    group = factor(
      phenotype,
      levels = c("Female", "Parae", "Yellow", "Immaculata")
    )
  ) %>%
  group_by(chromosome_type, phenotype) %>%
  summarise(
    mean_percentage = mean(CpG_pct, na.rm = TRUE),
    sem_percentage = sd(CpG_pct, na.rm = TRUE) / sqrt(sum(!is.na(CpG_pct))),
    n = sum(!is.na(CpG_pct)),
    .groups = "drop"
  ) 

# plot  
filename <- "normalized_counts_DML_PerChromosomeType_byGroup_barplot.pdf"


bar_plot_type <- ggplot(
  chromosome_type_summary,
  aes(x = chromosome_type, y = mean_percentage, fill = phenotype)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_percentage - sem_percentage,
      ymax = mean_percentage + sem_percentage
    ),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f ± %.3f",
        mean_percentage,
        sem_percentage
      ),
      y = mean_percentage + sem_percentage
    ),
    position = position_dodge(width = 0.8),
    angle = 45,
    hjust = 0,
    vjust = 0.5,
    size = 3
  ) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 100), clip = "off")+
  labs(
    title = "Normalized Counts of DML by Chromosome Type ~DNMTools~",
    x = "Chromosome type",
    y = "Normalized methylation (mCpG/CpG)",
    fill = "Phenotype"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot(bar_plot_type)


# Save plot
ggsave(
  paste0(plots, "/", filename),
  plot = bar_plot_type,
  width = 7,
  height = 4
)
