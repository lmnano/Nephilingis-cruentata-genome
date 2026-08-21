######################## SPIDER GENOMICS ############################

# load libraries
library(tidyverse)

##### SEX CHROMOSOME DETERMINATION NEPHILINGIS #####
setwd("chr_det/00_data/")
## load mosdepth summary files
files <- list.files(pattern = "mosdepth.summary.txt$", full.names = TRUE)
## read all files
cov_list <- lapply(files, function(f){
  df <- read.delim(f)
  # keep only chromosome-scale scaffolds
  df <- df %>%
    filter(grepl("^Chr", chrom))
  # extract sample name
  sample <- basename(f)
  sample <- sub(".mosdepth.summary.txt", "", sample)
  # assign sex manually
  # edit these IDs to match your samples
  sex <- case_when(
    sample %in% c("F13_0051","F13_0080","F13_0157","F13_0116","F13_0109") ~ "male",
    sample %in% c("F13_0018","F13_0008","F12_0279","F11_0576","F11_0574") ~ "female",
    TRUE ~ NA_character_)
  df %>%
    mutate(sample = sample, sex = sex) %>%
    select(sample, sex, chrom, mean)})
## combine all samples
cov <- bind_rows(cov_list)
## normalize coverage within each sample
# removes library size / sequencing depth effects
cov <- cov %>%
  group_by(sample) %>%
  mutate(norm_cov = mean / median(mean)) %>%
  ungroup()
## calculate mean normalized coverage by sex
sex_cov <- cov %>%
  group_by(sex, chrom) %>%
  summarise(mean_cov = mean(norm_cov),sd_cov = sd(norm_cov), n = n(), se_cov = sd_cov / sqrt(n), .groups = "drop")
## convert to wide format
sex_cov_wide <- sex_cov %>%
  select(sex, chrom, mean_cov) %>%
  pivot_wider(names_from = sex, values_from = mean_cov)
## calculate male:female coverage ratio
sex_cov_wide <- sex_cov_wide %>%
  mutate(male_female_ratio = male / female)
## inspect table
print(sex_cov_wide)

### FIGURE 1A: mean normalized coverage per chromosome by sex (barplot)
ggplot(sex_cov, aes(x = chrom, y = mean_cov, fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_cov - se_cov, ymax = mean_cov + se_cov), position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  theme_classic() +
  ylab("mean normalized coverage (± SE)") +
  xlab("chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("female" = "#E64B35", "male" = "#4DBBD5"))

### FIGURE 1B: male vs female difference per chromosome
sex_diff <- sex_cov_wide %>%
  mutate(delta = male - female, log2_ratio = log2(male / female))
ggplot(sex_diff, aes(x = chrom, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_col(fill = "black", alpha = 0.7) +
  theme_classic() +
  ylab("male - female normalized coverage") +
  xlab("chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

