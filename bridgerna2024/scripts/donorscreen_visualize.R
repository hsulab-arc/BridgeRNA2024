#setwd("~/ANALYSIS/28.DONOR_SCREEN/05.TRULY_FINAL/")
require(tidyverse)
require(ggseqlogo)
require(patchwork)
#require(combinat)


args = commandArgs(trailingOnly=TRUE)

OLIGO_COV <- read_tsv(args[1])
#OLIGO_COV <- read_tsv("~/ANALYSIS/28.DONOR_SCREEN/05.TRULY_FINAL/WORKDIR/08.screen_oligo_counts/screen_oligo_counts.tsv")
OLIGO_INFO <- read_tsv(args[2])
#OLIGO_INFO <- read_tsv("~/ANALYSIS/28.DONOR_SCREEN/05.TRULY_FINAL/WORKDIR/oligo_info.tsv")
SINGLE_MISMATCH_SETS <- read_tsv(args[3])
#SINGLE_MISMATCH_SETS <- read_tsv("~/ANALYSIS/28.DONOR_SCREEN/05.TRULY_FINAL/WORKDIR/single_mismatch_sets.tsv")
VISUALIZE_DIR <- args[4]
#VISUALIZE_DIR <- "~/ANALYSIS/28.DONOR_SCREEN/05.TRULY_FINAL/WORKDIR/09.visualize/"
#dir.create(VISUALIZE_DIR)

#####################
# Functions
#####################

nummismatches <- function(seq1, seq2){
  seq1_list <- strsplit(seq1, '')[[1]]
  seq2_list <- strsplit(seq2, '')[[1]]
  mm <- 0
  for (i in 1:length(seq1_list)){
    if (seq1_list[i] != seq2_list[i]){
      mm <- mm + 1
    }
  }
  return(mm)
}

make_motif <- function(seqs){

  counts <- matrix(0, 4, 11)
  for (s in seqs){
    ssplit <- strsplit(s, "")[[1]]
    for (i in 1:length(ssplit)){
      nuc <- ssplit[[i]]
      if (nuc == 'A'){
        counts[1, i] <- counts[1, i] + 1
      }else if (nuc == 'C'){
        counts[2, i] <- counts[2, i] + 1
      }else if (nuc == 'G'){
        counts[3, i] <- counts[3, i] + 1
      }else if (nuc == 'T'){
        counts[4, i] <- counts[4, i] + 1
      }
    }
  }

  motif <- counts / length(seqs)
  rownames(motif) <- c('A', 'C', 'G', 'T')
  return(motif)
}

reverse_complement <- function(seq){
  out <- c()
  seqlist <- strsplit(seq, "")[[1]]
  for (nuc in rev(seqlist)){
    if (nuc == 'A'){
      out <- c(out, 'T')
    }else if (nuc == 'T'){
      out <- c(out, 'A')
    }else if (nuc == 'G'){
      out <- c(out, 'C')
    }else if (nuc == 'C'){
      out <- c(out, "G")
    }else{
      out <- c(out, nuc)
    }
  }
  return(paste(out, collapse=''))
}

get_single_mutants <- function(seq, p){
  nucs <- c('A', 'C', 'G', 'T')
  seqlist <- strsplit(seq, '')[[1]]
  all_mutants <- c()
  for (nuc in nucs){
    out <- seqlist
    out[[p]] <- nuc
    all_mutants <- c(all_mutants, paste(out, collapse=''))
  }
  return(unique(all_mutants))
}

######################
# General CPM analysis
######################

# Convert barcode counts to counts per million (CPM)
umis.cpm <- OLIGO_COV %>%
  mutate(oligo_id=factor(oligo_id, levels=OLIGO_INFO$oligo_id)) %>%
  mutate(sample_id=factor(sample_id, levels=unique(sample_id))) %>%
  complete(sample_id, oligo_id, fill=list(oligo_count=0)) %>%
  separate(sample_id, into=c("screen", "condition", "biorep"), sep='-') %>%
  dplyr::select(-screen) %>%
  group_by(condition, biorep) %>%
  mutate(`cpm`=`oligo_count`/sum(`oligo_count`)*1e6) %>%
  dplyr::select(-oligo_count) %>% group_by()


# Plot CPM comparison across bioreps
ggplot(umis.cpm %>% spread(biorep, cpm),
       aes(x=log10(rep1+1), y=log10(rep2+1))) +
  geom_point(alpha=0.1) + theme_classic() +
  facet_grid(~condition) +
  xlab("rep1 - log10(CPM+1)") +
  ylab("rep2 - log10(CPM+1)") +
  geom_abline(linetype='dashed', color='blue') +
  coord_equal()
ggsave(paste0(VISUALIZE_DIR, '/cpm_biorep_comparison.pdf'), width=200*(4/3), height=100*(4/3), units='mm')

correl <- umis.cpm %>% mutate(cpm=log10(cpm+1)) %>%
  spread(biorep, cpm) %>%
  group_by(condition) %>%
  summarize(r=cor(rep1, rep2))
write_tsv(correl, paste0(VISUALIZE_DIR, '/bioreplicate_correlations.tsv'))

# Average CPM values across bioreps for each condition
umis.cpm.avg <- umis.cpm %>%
  spread(biorep, cpm) %>%
  mutate(avg_cpm=(rep1+rep2)/2) %>%
  dplyr::select(-rep1, -rep2) %>% group_by()

#######################################################
# Analysis of CPM Distribution in recombinant condition
#######################################################

# Correct recombinant CPM by control CPM
umis.cpm.avg.corrected <- umis.cpm.avg %>%
  spread(condition, avg_cpm) %>%
  mutate(control_expected=(1/n())*1e6) %>%
  mutate(correction_factor=control_expected/control) %>%
  mutate(recombinant_corrected=recombinant*correction_factor) %>%
  dplyr::select(oligo_id, cpm_corrected=recombinant_corrected)

# Calculate number of target/target loop mismatches for each oligo
oligo_mismatches.cpm <- OLIGO_INFO %>%
  mutate(donor_match=donor==donor_loop) %>%
  dplyr::select(oligo_category, oligo_id, donor_match) %>%
  filter(grepl("single_mismatch", oligo_category) |
           grepl("double_mismatch", oligo_category) |
           grepl('perfect_match', oligo_category) |
           oligo_category=='lib1_negative_control'
           ) %>%
    mutate(`Donor Loop/Donor Mismatches`=ifelse(
    donor_match, 0, ifelse(
      grepl("single_mismatch", oligo_category), 1, ifelse(
        grepl("double_mismatch", oligo_category), 2, 9
      )))) %>% dplyr::select(-donor_match) %>%
  left_join(umis.cpm.avg.corrected) %>%
  mutate(`Donor Loop/Donor Mismatches`=factor(`Donor Loop/Donor Mismatches`))

# Get number of oligos per mismatch category
oligo_mismatches.cpm.N <- oligo_mismatches.cpm %>%
  dplyr::count(`Donor Loop/Donor Mismatches`) %>%
  mutate(topy=c(3.8, 3.4, 2.5, 1.45))

# Plot the CPM distribution per mismatch category
ggplot(oligo_mismatches.cpm, aes(x=factor(`Donor Loop/Donor Mismatches`),
                                 y=log10(cpm_corrected+1))) +
  geom_violin(scale='width') +
  geom_text(data=oligo_mismatches.cpm.N,
            aes(x=`Donor Loop/Donor Mismatches`, label=n, y=topy)) +
  theme_classic() +
  xlab("Number of\nDonor Loop/Donor\nMismatches") +
  ylab("Selection - log10(CPM+1)")
ggsave(paste0(VISUALIZE_DIR, '/cpm_per_mismatch_category.pdf'), width=100*(4/3), height=100*(4/3), units='mm')

# Repeat but without any GO or WT sequences
oligo_mismatches.cpm <- OLIGO_INFO %>%
  mutate(donor_match=donor==donor_loop) %>%
  dplyr::select(oligo_category, oligo_id, donor_match) %>%
  filter(oligo_category %in% c('select100_single_mismatch', 'select5000_perfect_match', 'lib1_negative_control')) %>%
  mutate(`Donor Loop/Donor Mismatches`=ifelse(
    donor_match, 0, ifelse(
      grepl("single_mismatch", oligo_category), 1, ifelse(
        grepl("double_mismatch", oligo_category), 2, 9
      )))) %>% dplyr::select(-donor_match) %>%
  left_join(umis.cpm.avg.corrected) %>%
  mutate(`Donor Loop/Donor Mismatches`=factor(`Donor Loop/Donor Mismatches`))

# Get number of oligos per mismatch category
oligo_mismatches.cpm.N <- oligo_mismatches.cpm %>%
  dplyr::count(`Donor Loop/Donor Mismatches`) %>%
  mutate(topy=c(3.8, 3.4, 2.5))

ggplot(oligo_mismatches.cpm, aes(x=factor(`Donor Loop/Donor Mismatches`),
                                 y=log10(cpm_corrected+1))) +
  geom_violin(scale='width') +
  geom_boxplot(outlier.shape=NA, fill='light grey', width=0.2) +
  geom_text(data=oligo_mismatches.cpm.N,
            aes(x=`Donor Loop/Donor Mismatches`, label=n, y=topy)) +
  theme_classic() +
  xlab("Number of\nDonor Loop/Donor\nMismatches") +
  ylab("Selection - log10(CPM+1)")


########################################
# Single Mismatch Sets Analysis
########################################

# Extract target/target loop nucleotides and combine with CPM
single_mismatch_sets.cpm <- SINGLE_MISMATCH_SETS %>%
  inner_join(umis.cpm.avg.corrected) %>%
  filter(!is.na(cpm_corrected)) %>% group_by(baseseq) %>%
  filter(n()==144) %>%
  dplyr::rename(cpm=cpm_corrected)

# Calculate the percentage CPM for each mismatch set/pair
single_mismatch_sets.cpm_perc <- single_mismatch_sets.cpm %>%
  group_by(baseseq, mmpos, donor_loop_nuc) %>%
  mutate(cpm_perc=cpm/sum(cpm)) %>%
  dplyr::select(baseseq, mmpos, donor_loop_nuc, donor_nuc, cpm, cpm_perc) %>%
  group_by()

# Determine the top CPM quintile of single mismatch sets
top_quintile_oligo_set <- single_mismatch_sets.cpm_perc %>%
  group_by(baseseq) %>% filter(sum(is.nan(cpm_perc)) == 0) %>% group_by() %>%
  filter(donor_loop_nuc==donor_nuc) %>%
  group_by(baseseq) %>% filter(cpm==max(cpm)) %>% group_by() %>%
  arrange(desc(cpm)) %>% filter(row_number() < 0.2*n()) %>%
  dplyr::select(baseseq)

# Averages CPMs across each position in the top quintile
avg_cpm_perc_top_quintile <- single_mismatch_sets.cpm_perc %>%
  inner_join(top_quintile_oligo_set) %>%
  group_by(mmpos, donor_loop_nuc, donor_nuc) %>%
  summarize(cpm_perc=mean(cpm_perc)*100) %>% group_by() %>%
  mutate(mmpos=paste0("Position: ", mmpos)) %>%
  mutate(mmpos=factor(mmpos, levels=unique(mmpos))) %>%
  mutate(donor_loop_nuc=gsub("T", "U", donor_loop_nuc))

# Plot
ggplot(avg_cpm_perc_top_quintile %>%
         dplyr::rename(`% Recovered\nrecombinants`=cpm_perc),
       aes(x=donor_loop_nuc, y=donor_nuc, fill=`% Recovered\nrecombinants`)) +
  geom_tile() +
  facet_grid(~mmpos) +
  theme_classic() +
  scale_fill_gradient(low='white', high='red', limits=c(0, 100)) +
  xlab("Donor Guide Nucleotide") +
  ylab("Donor Nucleotide") +
  coord_equal()
ggsave(paste0(VISUALIZE_DIR, '/single_mismatch_efficiency_grid.pdf'), width=200*(4/3), height=50*(4/3), units='mm')


###########################################
# Top Quintile Enrichment - Perfect Matches
###########################################

# Get the perfect match oligos from 4 different categories
perfect_matches <- OLIGO_INFO %>% filter(oligo_category %in% c(
  "select100_single_mismatch", "select5000_perfect_match"
)) %>% filter(donor==donor_loop) %>% dplyr::select(oligo_id, donor) %>%
  left_join(umis.cpm.avg.corrected) %>%
  mutate(cpm_corrected=ifelse(is.na(cpm_corrected), 0, cpm_corrected)) %>%
  arrange(desc(cpm_corrected)) %>% filter(!is.infinite(cpm_corrected))

# Filter to the top CPM
perfect_matches.top_quintile <- perfect_matches %>% filter(row_number() < 0.2*n())

# Make a motif of top quintile targets
top_quintile_motif <- make_motif(perfect_matches.top_quintile$donor)
# Make a motif of all targets
background_motif <- make_motif(perfect_matches$donor)

# Create a PWM corrected by background frequencies
motif_pwm <- log2(top_quintile_motif/background_motif)

# Plot
ggseqlogo(motif_pwm, method='custom', seq_type='dna') +
  geom_hline(yintercept=0) +
  ggtitle("Enrichment of Nucleotides in Top CPM Quintile - WT Donor: ACAGTATCTTG") +
  ylab("log2(Enrichment)")
ggsave(paste0(VISUALIZE_DIR, '/enrichment_top_quintile.pdf'), width=120*(4/3), height=50*(4/3), units='mm')


#######################################################
# Distance from WT
#######################################################

# Calculate number of donor/donor loop mismatches for each oligo
wt_mismatches.cpm <- OLIGO_INFO %>%
  filter(donor==donor_loop) %>%
  filter(oligo_category != 'promoter_box_mutants') %>% rowwise() %>%
  mutate(wt_mismatches=nummismatches(donor, "ACAGTATCTTG")) %>%
  group_by() %>%
  dplyr::select(oligo_category, oligo_id, wt_mismatches) %>%
  left_join(umis.cpm.avg.corrected) %>%
  filter(!is.na(cpm_corrected)) %>% filter(!is.infinite(cpm_corrected))

# Get number of oligos per mismatch category
wt_mismatches.cpm.N <- wt_mismatches.cpm %>%
  dplyr::count(wt_mismatches)

# Plot the CPM distribution per mismatch category
ggplot(wt_mismatches.cpm %>% filter(wt_mismatches!=0), aes(x=factor(wt_mismatches),
                                 y=log10(cpm_corrected+1))) +
  geom_violin(scale='width', color='black', fill='#ade6ae') +
  geom_boxplot(outlier.shape=NA, fill='#e2f6e2', width=0.1) +
  geom_hline(yintercept=wt_mismatches.cpm %>% filter(wt_mismatches==0) %>% mutate(cpm_corrected=log10(cpm_corrected+1)) %>% .$cpm_corrected,
             linetype='dashed', color='red') +
  geom_text(data=wt_mismatches.cpm.N %>% filter(wt_mismatches!=0),
            aes(x=wt_mismatches, label=n, y=3.9)) +
  theme_classic() +
  xlab("Number of Mismatches with WT") +
  ylab("Read Abundance - log10(CPM+1)")
ggsave(paste0(VISUALIZE_DIR, '/cpm_distance_from_wildtype.pdf'), width=100*(4/3), height=100*(4/3), units='mm')

percentile <- wt_mismatches.cpm %>%
  filter(wt_mismatches==0 | oligo_category=='select5000_perfect_match') %>%
  arrange(desc(cpm_corrected)) %>%
  mutate(perc=(n()-row_number())/n()) %>%
  filter(wt_mismatches==0)
write_tsv(percentile, paste0(VISUALIZE_DIR, '/wildtype_single_mismatch_percentile_5k.tsv'))