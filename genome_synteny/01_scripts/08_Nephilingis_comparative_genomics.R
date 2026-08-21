######################## SPIDER GENOMICS ############################

# load libraries
library(tidyverse)
library(igraph)
library(lme4)
library(lmerTest)
library(ape)
library(patchwork)
library(phytools)
library(emmeans)
library(corHMM)
library(phylolm)
library(cocor)

##### COMPARATIVE GENOMICS #####
## organise species into groups with sex and without sex chromosomes
# Species with explicitly annotated sex chromosomes
species_with_known_sex <- c(
  "Amaurobius_ferox",
  "Argiope_bruennichi",
  "Gibbaranea_gibbosa",
  "Larinioides_cornutus",
  "Meta_bourneti",
  "Nephilingis_cruentata",
  "Uloborus_diversus",
  "Oedothorax_gibbosus",
  "Parasteatoda_lunata",
  "Tetragnatha_montana")

# Species currently lacking confirmed sex chromosome annotation treated as biologically unresolved
# absence of annotation ≠ absence of sex chromosomes
species_with_unknown_sex <- c(
  "Ectatosticta_davidi",
  "Pardosa_pseudoannulata",
  "Latrodectus_elegans",
  "Araneus_marmoreus",
  "Hylyphantes_graminicola",
  "Trichonephila_antipodiana",
  "Dysdera_silvatica")

# combined lookup table for explicit status tracking
species_sex_status <- tibble(species = c(species_with_known_sex, species_with_unknown_sex),
  sex_annotation_status = c(rep("known", length(species_with_known_sex)), rep("unknown", length(species_with_unknown_sex))))

### create edge list from all pairwise ref and target species tables (synteny results from cactus)
setwd("genome_synteny/00_data/02_cactus_data")
path <- getwd()
files <- list.files(path, pattern="_chr_match_table\\.tsv$", full.names=TRUE)
parse_chr_table <- function(file){
  df <- read.delim(file, check.names=FALSE, header=FALSE)
  ref_species <- as.character(df[1,1])
  ref_chr_labels <- as.character(df[1,])
  cat("\nProcessing:", ref_species, "| rows:", nrow(df), "| cols:", ncol(df), "\n")
  data_rows <- df[-1,]
  out <- list()
  row_idx <- seq(1,nrow(data_rows),by=2)
  for(r in row_idx){
    chr_row <- data_rows[r,]
    perc_row <- data_rows[r+1,]
    target_species <- as.character(chr_row[[1]])
    if(grepl("percCountMax", target_species)) next
    for(i in 2:ncol(data_rows)){
      out[[length(out)+1]] <- tibble(
        ref_species=ref_species,
        target_species=target_species,
        ref_chr=ref_chr_labels[i],
        target_chr=as.character(chr_row[[i]]),
        weight=as.numeric(perc_row[[i]]))}}
  bind_rows(out)}
all_edges <- map_dfr(files, parse_chr_table) %>%
  filter(!is.na(weight))
# exclude mitochondria
all_edges <- all_edges %>%
  filter(!grepl("^mt$", ref_chr, ignore.case = TRUE), !grepl("^mt$", target_chr, ignore.case = TRUE), !grepl("mito", ref_chr, ignore.case = TRUE), !grepl("mito", target_chr, ignore.case = TRUE))
# re-label confirmed sex chromosomes applied globally before all downstream analyses
# Nephilingis_cruentata: chr 2 = x1, chr 7 = x2
# Uloborus_diversus: chr 3 = x1, chr 10 = x2
all_edges <- map_dfr(files, parse_chr_table) %>%
  filter(!is.na(weight)) %>%
  mutate(ref_chr = case_when(
      ref_species == "Nephilingis_cruentata" & ref_chr == "2" ~ "x1",
      ref_species == "Nephilingis_cruentata" & ref_chr == "7" ~ "x2",
      ref_species == "Uloborus_diversus" & ref_chr == "3" ~ "x1",
      ref_species == "Uloborus_diversus" & ref_chr == "10" ~ "x2",TRUE ~ ref_chr),
         target_chr = case_when(
      target_species == "Nephilingis_cruentata" & target_chr == "2" ~ "x1",
      target_species == "Nephilingis_cruentata" & target_chr == "7" ~ "x2",
      target_species == "Uloborus_diversus" & target_chr == "3" ~ "x1",
      target_species == "Uloborus_diversus" & target_chr == "10" ~ "x2",
      TRUE ~ target_chr),
    ref_node = paste(ref_species, ref_chr, sep="_"),
    target_node = paste(target_species, target_chr, sep="_"))

## filter strong edges
strong_edges <- all_edges %>% filter(weight>=0.5)
#write.csv(strong_edges,"sex_chr_edge_list.csv",row.names=FALSE)
# each edge = synteny signal between chromosomes
# weight = conserved gene order fraction

### create chromosome homology groups ###
g <- graph_from_data_frame(strong_edges %>% select(ref_node,target_node,weight), directed=FALSE)
E(g)$weight <- strong_edges$weight
cl <- cluster_louvain(g, weights=E(g)$weight)
membership_df <- tibble(chromosome=names(membership(cl)), homolog_group=membership(cl))
membership_df$has_X <- grepl("x[0-9]+$", membership_df$chromosome, ignore.case=TRUE)
#write.csv(membership_df,"chromosome_homolog_groups_strong_edges.csv",row.names=FALSE)

### create homology groups including all edges ###
## "un"filter strong edges
### USE ALL EDGES (overwrite strong_edges for compatibility) ###
strong_edges <- all_edges
### create chromosome homology groups ###
g <- graph_from_data_frame(strong_edges %>% select(ref_node,target_node,weight), directed=FALSE)
E(g)$weight <- strong_edges$weight
cl <- cluster_louvain(g, weights=E(g)$weight)
membership_df <- tibble(chromosome=names(membership(cl)), homolog_group=membership(cl))
membership_df$has_X <- grepl("x[0-9]+$", membership_df$chromosome, ignore.case=TRUE)
#write.csv(membership_df,"chromosome_homolog_groups_all_edges.csv",row.names=FALSE)

### sex chromosome homology groups ###
sex_groups_strict <- membership_df %>%
  filter(has_X) %>%
  pull(homolog_group) %>%
  unique()
sex_groups_prop <- sex_groups_strict

### annotate edges ###
edge_annot <- all_edges %>%
  left_join(membership_df,by=c("ref_node"="chromosome")) %>% rename(ref_group=homolog_group) %>%
  left_join(membership_df,by=c("target_node"="chromosome")) %>% rename(target_group=homolog_group)

### TABLE SX: HOMOLOGY GROUPS AND CORRESPONDING CHROMOOMES FOR EACH SPECIES ###
Table_SX_homology_groups_species_chromosome_correspondence <- membership_df %>%
  separate(chromosome, into = c("species", "chromosome"), sep = "_(?=[^_]+$)") %>%
  group_by(species, homolog_group) %>%
  summarise(chromosomes = paste(sort(unique(chromosome)), collapse = ";"), .groups = "drop") %>%
  pivot_wider(names_from = homolog_group, values_from = chromosomes)
Table_SX_homology_groups_species_chromosome_correspondence
#write.csv(Table_SX_homology_groups_species_chromosome_correspondence, "Table_SX_homology_groups_species_chromosome_correspondence.csv", row.names = FALSE)

### FIGURE X: CONSERVATION OF X CHROMOSOME HOMOLOGY GROUPS
## prepare chromosome homology table using propagated sex chromosome assignments
homology_plot <- membership_df %>%
  mutate(species=sub("_[^_]+$","",chromosome),
         chromosome=sub("^.*_","",chromosome),
         chromosome_type=ifelse(homolog_group %in% sex_groups_prop,"sex chromosome","autosome")) %>%
  group_by(species,homolog_group,chromosome_type) %>%
  summarise(n_chromosomes=n(),.groups="drop")
## plot number of chromosomes per homology group
homology_plot %>%
  ggplot(aes(x=homolog_group,y=species,size=n_chromosomes,fill=chromosome_type))+
  geom_point(shape=21,color="black")+
  scale_fill_manual(values=c("autosome"="#BDBDBD","sex chromosome"="#E64B35"))+
  scale_size_continuous(range=c(2,8))+
  theme_classic()+
  xlab("Chromosome homology group")+
  ylab("")+
  theme(axis.text.y=element_text(size=8),
        axis.text.x=element_text(angle=90,hjust=1))


#### ANALYSIS WITH ONLY SPECIES THAT HAVE RELIABLE SEX CHROMOSOME IDENTIFICATION ###
### species-level annotation (uses only known set)
edge_annot_strict <- edge_annot %>%
  mutate(ref_class = ifelse(ref_chr %in% c("x1","x2","x3"), "sex", "auto"),
    target_class = ifelse(target_chr %in% c("x1","x2","x3"), "sex", "auto"),
    pair_class = paste(pmin(ref_class, target_class), pmax(ref_class, target_class), sep="-"),
    ref_has_known_sex = ref_species %in% species_with_known_sex,
    target_has_known_sex = target_species %in% species_with_known_sex,
    ref_unknown_sex = ref_species %in% species_with_unknown_sex,
    target_unknown_sex = target_species %in% species_with_unknown_sex)
### filter here only species that have sex chromosomes identified reliably ###
edge_annot_strict_filtered <- edge_annot_strict %>%
  filter(ref_has_known_sex & target_has_known_sex)
# BOTH species must have known sex chromosomes

### transition candidates (expected empty if high degree of conservation)
sex_transition_candidates <- edge_annot_strict_filtered %>%
  filter(pair_class == "auto-sex") %>%
  group_by(ref_species, target_species, ref_chr, target_chr) %>%
  summarise(mean_weight = mean(weight), n = n(), .groups="drop") %>%
  arrange(desc(mean_weight))
#write.csv(sex_transition_candidates, "sex_autosome_transition_candidates_filtered.csv", row.names=FALSE)
# interpretation:
# this is a *synteny asymmetry table*, not a transition inference table
# RESULT: Table empty, so no switches between sex chromosomes and autosomes 

### PROPAGATED ANALYSIS ###
edge_annot_prop <- edge_annot_strict_filtered %>%   # <<< FIX IS HERE
  mutate(ref_class=ifelse(ref_group %in% sex_groups_prop,"sex","auto"),
    target_class=ifelse(target_group %in% sex_groups_prop,"sex","auto")) %>%
  mutate(pair_class=paste(pmin(ref_class,target_class),
                          pmax(ref_class,target_class),sep="-"))

### STATISTICS 1:
model_final <- lmer(weight ~ pair_class + (1|ref_species) + (1|target_species), data=edge_annot_prop)
summary(model_final)
anova(model_final)

# FIGURE 1:
# statistical comparison of synteny conservation across:
# auto-auto | auto-sex | sex-sex
# hypothesis: sex chromosomes differ in evolutionary constraint regime
emmeans(model_final,pairwise~pair_class)

### FIGURE 3A: DISTRIBUTION OF SYNTENY ###
edge_annot_prop %>%
  ggplot(aes(pair_class,weight))+
  geom_boxplot()+
  geom_jitter(width=0.15,alpha=0.15)+
  theme_classic()+
  ylab("Synteny / homology weight")+
  xlab("Chromosome comparison class")
# interpretation:
# direct visualization of homology conservation distributions

### HOMOLOGY STABILITY ###
edge_annot_prop %>%
  filter(!is.na(ref_group),!is.na(target_group)) %>%
  mutate(same_group=ref_group==target_group) %>%
  group_by(pair_class) %>%
  summarise(stability=mean(same_group),n=n(),.groups="drop")
# interpretation:
# stability = conservation of inferred chromosomal identity

### PHYLOGENY ###
tree_raw <- read.nexus(file.path(path,"iqtree_full_run_ML_consensus.tree"))
tree <- if(inherits(tree_raw,"multiPhylo")) tree_raw[[1]] else tree_raw
tree_mid <- midpoint.root(tree)
species_trait <- edge_annot_prop %>%
  group_by(ref_species) %>%
  summarise(trait=mean(weight),.groups="drop")
trait_vec <- species_trait$trait
names(trait_vec) <- species_trait$ref_species
trait_vec <- trait_vec[tree_mid$tip.label]

phylosig(tree_mid,trait_vec,method="K",test=TRUE)
# interpretation:
# phylogenetic structuring of synteny conservation

## Plot the tree
sex_state <- edge_annot_prop %>%
  group_by(ref_species) %>%
  summarise(n_sex=sum(ref_class=="sex"),.groups="drop")

plot(tree_mid,cex=0.7)
tiplabels(sex_state$n_sex[match(tree_mid$tip.label,sex_state$ref_species)],
          frame="none",adj=-0.2)
# interpretation:
# exploratory mapping of inferred sex-linked chromosomal content

### SIZE EVOLUTION ###
# rename chromosome of Nephilingis and Uloborus
chr_sizes <- read.csv("chromosome_lengths.csv") %>%
  mutate(chr = case_when(
    species == "Nephilingis_cruentata" & chr == "2" ~ "x1",
    species == "Nephilingis_cruentata" & chr == "7" ~ "x2",
    species == "Uloborus_diversus" & chr == "3" ~ "x1",
    species == "Uloborus_diversus" & chr == "10" ~ "x2",
    TRUE ~ as.character(chr)))
edge_sizes <- edge_annot_prop %>%
  left_join(chr_sizes, by=c("ref_species"="species","ref_chr"="chr")) %>%
  rename(ref_length=length) %>%
  left_join(chr_sizes, by=c("target_species"="species","target_chr"="chr")) %>%
  rename(target_length=length) %>%
  filter(!is.na(ref_length),!is.na(target_length))
edge_sizes <- edge_sizes %>%
  mutate(size_diff=abs(ref_length-target_length),
         log_size_diff=abs(log(ref_length+1)-log(target_length+1)))

### STATISTICS 2: CHROMOSOME LENGTH DIVERGENCE ###
model_size <- lmer(log_size_diff ~ pair_class +
                     (1|ref_species) +
                     (1|target_species),
                   data=edge_sizes)
summary(model_size)
anova(model_size)
emmeans(model_size, pairwise ~ pair_class)
# interpretation: tests whether homologous chromosome pairs differ in structural size divergence between:
# auto-auto = homologous autosome pairs
# sex-sex   = homologous sex chromosome pairs
# lower values = stronger chromosome size conservation
# if sex-sex > auto-auto: sex chromosomes show elevated structural divergence
# if sex-sex < auto-auto: sex chromosomes are more structurally conserved

### FIGURE 3B: Size differences between auto and sex chromosomes ###
edge_sizes %>%
  ggplot(aes(pair_class,log_size_diff))+
  geom_boxplot()+
  geom_jitter(width=0.2,alpha=0.2)+
  theme_classic()
# interpretation: structural chromosome divergence across classes

## keep only homologous chromosome comparisons
edge_sizes <- edge_annot_prop %>%
  filter(ref_group == target_group) %>%
  left_join(chr_sizes, by=c("ref_species"="species","ref_chr"="chr")) %>%
  rename(ref_length=length) %>%
  left_join(chr_sizes, by=c("target_species"="species","target_chr"="chr")) %>%
  rename(target_length=length) %>%
  filter(!is.na(ref_length), !is.na(target_length))
## chromosome size divergence
edge_sizes <- edge_sizes %>%
  mutate(size_diff = abs(ref_length-target_length), log_size_diff = abs(log(ref_length+1)-log(target_length+1)))

### STATISTICS 2b: CHROMOSOME LENGTH DIVERGENCE (ONLY HOMOLGOUS)###
model_size <- lmer(log_size_diff ~ pair_class + (1|ref_species) + (1|target_species), data=edge_sizes)
summary(model_size)
anova(model_size)
emmeans(model_size, pairwise ~ pair_class)

### FIGURE 3 ###
edge_sizes %>%
  ggplot(aes(pair_class, log_size_diff))+
  geom_boxplot()+
  geom_jitter(width=0.2, alpha=0.2)+
  theme_classic()
# interpretation: structural chromosome divergence across homologous chromosome classes

### CONTROL ANALYSIS: ONLY SPECIES WITH TWO SEX CHROMOSOMES ###
## Derive species-level counts of sex-linked chromosomes from homology groups (based on synteny-inferred sex chromosome clusters)
sex_group_counts <- membership_df %>%
  filter(homolog_group %in% sex_groups_prop) %>%
  separate(chromosome, into = c("species", "chr"), sep = "_(?=[^_]+$)") %>%
  group_by(species, homolog_group) %>%
  summarise(n_chr = n(), .groups = "drop")
species_sex_counts <- sex_group_counts %>%
  group_by(species) %>%
  summarise(total_sex_chr = sum(n_chr),
            n_sex_groups = n_distinct(homolog_group),
            .groups = "drop")
## identify species with exactly two inferred sex chromosomes
twoX_species <- species_sex_counts %>%
  filter(total_sex_chr == 2) %>%
  pull(species)
## restrict dataset
edge_sizes_2X <- edge_sizes %>%
  filter(ref_species %in% twoX_species,
         target_species %in% twoX_species)
### STATISTICS 2c: LENGTH DIVERGENCE IN TWO-X SPECIES ONLY ###
model_size_2X <- lmer(log_size_diff ~ pair_class +
                        (1|ref_species) +
                        (1|target_species), data=edge_sizes_2X)
summary(model_size_2X)
anova(model_size_2X)
emmeans(model_size_2X, pairwise ~ pair_class)
# interpretation: control analysis removing potential bias from species with 1 or 3 sex chromosomes tests whether elevated sex chromosome divergence remains within uniform 2X systems

### FIGURE 3b ###
edge_sizes_2X %>%
  ggplot(aes(pair_class, log_size_diff))+
  geom_boxplot()+
  geom_jitter(width=0.2, alpha=0.2)+
  theme_classic()
# interpretation: chromosome size divergence restricted to species with two sex chromosomes


#### Phylogenetic test for sex chromosome evolution
### SPECIES-LEVEL SEX STATE FROM HOMOLOGY (PREFERRED) ###
species_sex_state <- edge_annot_prop %>%
  mutate(ref_is_sex = ref_group %in% sex_groups_prop,
         target_is_sex = target_group %in% sex_groups_prop) %>%
  select(ref_species, ref_is_sex, target_species, target_is_sex) %>%
  pivot_longer(cols = c(ref_species, target_species), names_to = "role", values_to = "species") %>%
  mutate(sex_state = c(ref_is_sex, target_is_sex)[seq_len(n())]) %>%
  group_by(species) %>%
  summarise(has_sex_chr = any(sex_state, na.rm = TRUE), n_sex_edges = sum(sex_state, na.rm = TRUE), .groups = "drop")


#### SEX CHROMOSOME COPY NUMBER EVOLUTION ####
## count sex-linked homologous groups per species (derived from synteny-based homolog clustering)
sex_group_counts <- membership_df %>%
  filter(homolog_group %in% sex_groups_prop) %>%
  separate(chromosome, into=c("species","chr"), sep="_(?=[^_]+$)") %>%
  group_by(species, homolog_group) %>%
  summarise(n_chr=n(), .groups="drop")
# write it out
#write.csv(sex_group_counts, "sex_homolog_group_copy_numbers.csv", row.names=FALSE)
## collapse to species-level sex chromosome copy number
species_sex_counts <- sex_group_counts %>%
  group_by(species) %>%
  summarise(total_sex_chr=sum(n_chr),
            n_sex_groups=n_distinct(homolog_group),
            .groups="drop")
count_vec <- species_sex_counts$total_sex_chr
names(count_vec) <- species_sex_counts$species
## align trait to tree
count_vec_aligned <- count_vec[tree_mid$tip.label]
## make tree ultrametric
tree_mid <- chronos(tree_mid)

### phylogenetic discrete models (chromosome count evolution)
## corHMM input
dat_corhmm <- data.frame(species = names(count_vec_aligned), state = count_vec_aligned)
## match tree and data
tree_cor <- drop.tip(tree_mid, setdiff(tree_mid$tip.label, dat_corhmm$species))
dat_corhmm <- dat_corhmm[match(tree_cor$tip.label, dat_corhmm$species), ]
dat_corhmm$state <- as.factor(dat_corhmm$state)
## model 1: equal rates (er)
fit_ER <- corHMM(phy = tree_cor, data = dat_corhmm, rate.cat = 1, model = "ER")
## model 2: asymmetric rates (ard)
fit_AR <- corHMM(phy = tree_cor, data = dat_corhmm, rate.cat = 1, model = "ARD")
## model 3: use corHMM
# prepare data frame for corHMM
## ordered state space
states <- sort(unique(dat_corhmm$state))
k <- length(states)
## build adjacency-constrained rate matrix
Q <- matrix(0, k, k)
rownames(Q) <- colnames(Q) <- states
for (i in 1:(k - 1)) {
  Q[i, i + 1] <- 1
  Q[i + 1, i] <- 1}
diag(Q) <- 0
fit_STEP <- corHMM(phy = tree_cor, data = dat_corhmm, rate.cat = 1, model = "custom", rate.mat = Q)
## run corHMM with fixed topology constraint
fit_STEP <- corHMM(phy = tree_cor, data = dat_corhmm, rate.cat = 1, model = "custom", rate.mat = Q)
### model comparison 
model_table <- data.frame(model = c("ER", "ARD", "STEP"), logLik = c(fit_ER$loglik, fit_AR$loglik, fit_STEP$loglik), AIC = c(fit_ER$AIC, fit_AR$AIC, fit_STEP$AIC))
print(model_table)
## choose best model (lowest AIC)
best_fit <- if (fit_ER$AIC < fit_AR$AIC) fit_ER else fit_AR

### ancestral state reconstruction (best model)
## extract marginal ancestral states from best model
anc_states <- best_fit$states
### FIGURE X: ancestral reconstruction of chromosome number
## colour palette for states (more distinct, publication-style)
state_cols <- c("1" = "#edf6f9", "2" = "#83c5be", "3" = "#006d77")
plot(tree_cor, cex = 0.7, label.offset = 0.02)
## coloured tip symbols (circles only, slightly larger than nodes)
tip_state <- as.character(count_vec_aligned[tree_cor$tip.label])
tiplabels(pch = 21, bg = state_cols[tip_state], col = "black", cex = 2.5)
## larger ancestral pies (slightly smaller than tips for contrast)
nodelabels(pie = anc_states, piecol = state_cols, cex = 0.6)
# interpretation: circles at tips = observed sex chromosome copy number per species node pies = ancestral state probabilities under corHMM model

### CONTINUOUS TRAIT: SEX CHROMOSOME LENGTH MAPPED ON PHYLOGENY
### define sex chromosomes based on HOMOLOGY (not annotation)
sex_chr_lengths <- membership_df %>%
  filter(homolog_group %in% sex_groups_prop) %>%
  separate(chromosome, into = c("species", "chr"), sep = "_(?=[^_]+$)") %>%
  left_join(chr_sizes, by = c("species", "chr")) %>%
  filter(!is.na(length)) %>%
  group_by(species) %>%
  summarise(total_sex_chr_length = sum(length), .groups = "drop")
## build named vector
sex_vec <- sex_chr_lengths$total_sex_chr_length
names(sex_vec) <- sex_chr_lengths$species
## make tree ultrametric
tree_ultra <- chronos(tree_mid)
## confirm ultrametric conversion worked
stopifnot(is.ultrametric(tree_ultra))
## align to tree
common <- intersect(names(sex_vec), tree_ultra$tip.label)
tree_len <- drop.tip(tree_ultra, setdiff(tree_ultra$tip.label, common))
sex_vec <- sex_vec[tree_len$tip.label]
stopifnot(all(names(sex_vec) == tree_len$tip.label))
## log trasnform chromosome size
sex_log <- log(sex_vec + 1)
### FIGURE X: ancestral reconstruction of chromosome size
anc <- phytools::fastAnc(tree = tree_len, x = sex_log, vars = TRUE, CI = TRUE)
## build contMap
obj <- phytools::contMap(tree_len, sex_log, plot = FALSE, res = 200)
## use stable colour gradient
my_cols <- colorRampPalette(c("#edf6f9", "#83c5be", "#006d77"))(200)
obj <- phytools::setMap(obj, colors = my_cols)
## plot mapped tree
plot(obj, fsize = 0.7, lwd = 4, legend = 0.8 * max(nodeHeights(tree_ultra)))

### check if sex chromosome length correlates with genome size and repeat content across the genome
## load table with genome size and repeat content
path <- "genome_synteny/00_data/02_cactus_data/Table_S5_spider_chromosomes.csv"
dat <- read.csv(path, stringsAsFactors = FALSE)
dat
## clean species 
clean_species <- function(x) {
  x %>%
    trimws() %>%
    gsub("_+$", "", .) %>%   # remove trailing underscores
    gsub("\\s+", "_", .)}
dat <- dat %>%
  mutate(species = clean_species(species))
sex_chr_lengths <- sex_chr_lengths %>%
  mutate(species = clean_species(species))
## clean species names and create repeat fraction
dat <- dat %>%
  mutate(species = gsub(" ", "_", species), genome_bp = as.numeric(genome_size) * 1e9,
    repeat_fraction = as.numeric(repeat_content) / 100)
## merge datasets with sex chromosome length and genome size 
### join sex chromosome lengths
dat2 <- dat %>%
  left_join(sex_chr_lengths, by = "species") %>%
  mutate(sex_chr_bp = total_sex_chr_length, sex_chr_fraction = sex_chr_bp / genome_bp) %>%
  filter(!is.na(genome_bp), !is.na(sex_chr_bp))
nrow(dat2)
setdiff(dat$species, sex_chr_lengths$species)
## check summary
summary(dat2$genome_bp)
summary(dat2$sex_chr_bp)
summary(dat2$repeat_fraction)

### SEX CHROMOSOME SIZE vs GENOME SIZE
m_genome <- lm(sex_chr_bp ~ genome_bp, data = dat2)
summary(m_genome)
anova(m_genome)
# visualization
p0 <- ggplot(dat2, aes(genome_bp / 1e6, sex_chr_bp / 1e6)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  theme_classic() +
  xlab("Genome size (bp)") +
  ylab("Sex chromosome size (bp)")

## SEX CHROMOSOME SIZE vs REPEAT CONTENT
dat_rep <- dat2 %>%
  filter(!is.na(repeat_fraction))
m_repeat <- lm(sex_chr_bp ~ repeat_fraction, data = dat_rep)
summary(m_repeat)
anova(m_repeat)
# visualization
p1 <- ggplot(dat_rep, aes(repeat_fraction, sex_chr_bp)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  theme_classic() +
  xlab("Repeat fraction") +
  ylab("Sex chromosome size (bp)")

## GEONOME SIZE vs REPEAT CONTENT
m_repeat_genome <- lm(repeat_fraction ~ genome_bp, data = dat2)
summary(m_repeat_genome)
anova(m_repeat_genome)

# visualization
p2 <- ggplot(dat2, aes(genome_bp, repeat_fraction)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  theme_classic() +
  xlab("Genome size (bp)") +
  ylab("Repeat fraction")

(p0 | p1 | p2)

### DO THE SAME WITH OWN REPEAT LIBRARY
## SEX CHROMOSOME SIZE vs REPEAT CONTENT
dat_rep <- dat2 %>%
  filter(!is.na(masked_perc))
m_repeat <- lm(masked_perc ~ sex_chr_bp, data = dat_rep)
summary(m_repeat)
anova(m_repeat)
# visualization
p1 <- ggplot(dat_rep, aes(sex_chr_bp / 1e6, masked_perc)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  theme_classic() +
  xlab("Sex chromosome size (Mbp)") +
  ylab("Repeat fraction")
## GENOME SIZE vs REPEAT CONTENT
m_repeat_genome <- lm(masked_perc ~ genome_bp, data = dat2)
summary(m_repeat_genome)
anova(m_repeat_genome)
# visualization
p2 <- ggplot(dat2, aes(genome_bp / 1e6, masked_perc)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  theme_classic() +
  xlab("Genome size (Mbp)") +
  ylab("Repeat fraction")

(p0 | p1 | p2)

### CHECK IF AUTOSOMES CORRELATE LESS WITH REPEAT CONTENT THAN SEX CHROMOSOMES
## Calculate the total autosomal genome size by subtracting the total sex chromosome size from the total genome size
dat2 <- dat2 %>%
  mutate(autosome_bp = genome_bp - sex_chr_bp)
## Check the resulting autosomal genome sizes
summary(dat2$autosome_bp)
## Check for biologically impossible negative values
dat2 %>%
  filter(autosome_bp < 0)
## Prepare a common dataset containing only species with complete data
## for all variables used in the comparison
dat_test <- dat2 %>%
  filter(!is.na(masked_perc), !is.na(sex_chr_bp), !is.na(autosome_bp))
## SEX CHROMOSOME SIZE vs REPEAT CONTENT
m_sex <- lm(masked_perc ~ sex_chr_bp, data = dat_test)
summary(m_sex)
anova(m_sex)
## AUTOSOMAL GENOME SIZE vs REPEAT CONTENT
m_auto <- lm(masked_perc ~ autosome_bp, data = dat_test)
summary(m_auto)
anova(m_auto)
## COMPARE THE TWO CORRELATIONS
## Correlation between repeat content and sex chromosome size
r_y_sex <- cor(dat_test$masked_perc, dat_test$sex_chr_bp)
## Correlation between repeat content and autosomal genome size
r_y_auto <- cor(dat_test$masked_perc, dat_test$autosome_bp)
## Correlation between the two predictors
## This is required for the test of two overlapping correlations
r_sex_auto <- cor(dat_test$sex_chr_bp, dat_test$autosome_bp)
r_y_sex
r_y_auto
r_sex_auto
## Test whether the association between repeat content and sex chromosome
## size differs significantly from the association between repeat content
## and autosomal genome size
cocor.dep.groups.overlap(r.jk = r_y_sex, r.jh = r_y_auto, r.kh = r_sex_auto, n = nrow(dat_test), test = "steiger1980")

### VISUALIZATION
## Repeat content vs sex chromosome size
p1 <- ggplot(dat_test, aes(sex_chr_bp / 1e6, masked_perc)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  theme_classic() +
  xlab("Sex chromosome size (Mbp)") +
  ylab("Repeat fraction")
## Repeat content vs autosomal genome size
p2 <- ggplot(dat_test, aes(autosome_bp / 1e6, masked_perc)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  theme_classic() +
  xlab("Autosomal genome size (Mbp)") +
  ylab("Repeat fraction")

## Display the two plots side by side
p1 | p2