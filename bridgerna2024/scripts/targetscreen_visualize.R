require(tidyverse)
require(ggseqlogo)
require(patchwork)

args = commandArgs(trailingOnly=TRUE)

OLIGO_COV <- read_tsv(args[1])
#OLIGO_COV <- read_tsv("~/ANALYSIS/21.POOLED_SPECIFICITY_ASSAY/06.PAPER_FINAL/WORKDIR/04.oligo_coverage/all_oligo_coverage.tsv")
OLIGO_INFO <- read_tsv(args[2])
#OLIGO_INFO <- read_tsv("~/ANALYSIS/21.POOLED_SPECIFICITY_ASSAY/06.PAPER_FINAL/WORKDIR/oligo_info.tsv")
VISUALIZE_DIR <- args[3]
#VISUALIZE_DIR <- "~/ANALYSIS/21.POOLED_SPECIFICITY_ASSAY/06.PAPER_FINAL/WORKDIR/05.visualize/"
#dir.create(VISUALIZE_DIR)


#####################
# Functions
#####################

mismatchpos <- function(seq){
  pos <- unique(as.numeric(str_locate_all(seq, "N")[[1]]))
  if (length(pos) == 1){
    return(4)
  }else{
    return(pos[pos!=4])
  }
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


######################
# General CPM analysis
######################

OLIGO_COV %>% 
  mutate(barcode=factor(barcode, levels=OLIGO_INFO$barcode)) %>% 
  mutate(library=factor(library, levels=unique(library))) %>% 
  complete(library, barcode, fill=list(barcode_count=0)) %>% 
  filter(grepl("Recomb", library)) %>% group_by(library) %>% 
  summarize(dropout_perc=sum(barcode_count==0)/n()*100)

# Convert barcode counts to counts per million (CPM)
barcodes.cpm <- OLIGO_COV %>% 
  mutate(barcode=factor(barcode, levels=OLIGO_INFO$barcode)) %>% 
  mutate(library=factor(library, levels=unique(library))) %>% 
  complete(library, barcode, fill=list(barcode_count=0)) %>% 
  separate(library, into=c("screen", "condition", "biorep"), sep='-') %>%
  dplyr::select(-screen) %>%
  group_by(condition, biorep) %>% 
  mutate(`cpm`=`barcode_count`/sum(`barcode_count`)*1e6) %>% 
  dplyr::select(-barcode_count) %>% group_by()

# Plot CPM comparison across bioreps
ggplot(barcodes.cpm %>% spread(biorep, cpm), 
       aes(x=log10(rep1+1), y=log10(rep2+1))) +
  geom_point(alpha=0.1) + theme_classic() +
  facet_grid(~condition) + 
  xlab("rep1 - log10(CPM+1)") +
  ylab("rep2 - log10(CPM+1)") +
  geom_abline(linetype='dashed', color='blue')
ggsave(paste0(VISUALIZE_DIR, '/cpm_biorep_comparison.pdf'), width=200*(4/3), height=100*(4/3), units='mm')

# Average CPM values across bioreps for each condition
barcodes.cpm.avg <- barcodes.cpm %>% spread(biorep, cpm) %>% 
  mutate(avg_cpm=(rep1+rep2)/2) %>% 
  dplyr::select(-rep1, -rep2) %>% group_by()

#######################################################
# Analysis of CPM Distribution in recombinant Condition
#######################################################

# Correct recombinant CPM by control CPM
barcodes.cpm.avg.corrected <- barcodes.cpm.avg %>% 
  spread(condition, avg_cpm) %>% 
  mutate(control_expected=(1/n())*1e6) %>%
  mutate(correction_factor=control_expected/control) %>% 
  mutate(recombinant_corrected=recombinant*correction_factor) %>% 
  dplyr::select(barcode, cpm_corrected=recombinant_corrected)

# Calculate number of target/target loop mismatches for each oligo
oligo_mismatches.cpm <- OLIGO_INFO %>% mutate(target_match=target==target_loop) %>% 
  dplyr::select(oligo_category, barcode, target_match) %>% 
  filter(oligo_category %in% c(
    "single_mismatch_pairs_other", "single_mismatch_pairs_saturated", "total_mismatch_pairs",
    "double_mismatch_pairs"
  )) %>% mutate(`Target Loop/Target Mismatches`=ifelse(
    target_match, 0, ifelse(
      grepl("single_mismatch", oligo_category), 1, ifelse(
        oligo_category=='total_mismatch_pairs', 9, 2
      )))) %>% dplyr::select(-target_match) %>% 
  left_join(barcodes.cpm.avg.corrected) %>% 
  mutate(`Target Loop/Target Mismatches`=factor(`Target Loop/Target Mismatches`))

# Get number of oligos per mismatch category
oligo_mismatches.cpm.N <- oligo_mismatches.cpm %>% 
  dplyr::count(`Target Loop/Target Mismatches`) %>% 
  mutate(topy=c(3.8, 3.4, 2.5, 1.45))

# Plot the CPM distribution per mismatch category
ggplot(oligo_mismatches.cpm, aes(x=factor(`Target Loop/Target Mismatches`), 
                                 y=log10(cpm_corrected+1))) +
  geom_violin(scale='width') +
  geom_boxplot(outlier.shape=NA, fill='light grey', width=0.2) + 
  geom_text(data=oligo_mismatches.cpm.N, 
            aes(x=`Target Loop/Target Mismatches`, label=n, y=topy)) +
  theme_classic() +
  xlab("Number of\nTarget Loop/Target\nMismatches") + 
  ylab("Selection - log10(CPM+1)")
ggsave(paste0(VISUALIZE_DIR, '/cpm_per_mismatch_category.pdf'), width=100*(4/3), height=100*(4/3), units='mm')

# T test of 0 mismatch category vs. 9 mismatch category
t.test(oligo_mismatches.cpm %>% filter(`Target Loop/Target Mismatches`==0) %>% 
         mutate(cpm_corrected=cpm_corrected+1) %>% 
         .$cpm_corrected %>% log10,
       oligo_mismatches.cpm %>% filter(`Target Loop/Target Mismatches`==9) %>% 
         mutate(cpm_corrected=cpm_corrected+1) %>% 
         .$cpm_corrected %>% log10, alternative='greater')

# T test of 1 mismatch category vs. 9 mismatch category
t.test(oligo_mismatches.cpm %>% filter(`Target Loop/Target Mismatches`==1) %>% 
         mutate(cpm_corrected=cpm_corrected+1) %>% 
         .$cpm_corrected %>% log10,
       oligo_mismatches.cpm %>% filter(`Target Loop/Target Mismatches`==9) %>% 
         mutate(cpm_corrected=cpm_corrected+1) %>% 
         .$cpm_corrected %>% log10, alternative='greater')

# T test of 2 mismatch category vs. 9 mismatch category
t.test(oligo_mismatches.cpm %>% filter(`Target Loop/Target Mismatches`==2) %>% 
         mutate(cpm_corrected=cpm_corrected+1) %>% 
         .$cpm_corrected %>% log10,
       oligo_mismatches.cpm %>% filter(`Target Loop/Target Mismatches`==9) %>% 
         mutate(cpm_corrected=cpm_corrected+1) %>% 
         .$cpm_corrected %>% log10, alternative='greater')


########################################
# Single Mismatch Sets Analysis
########################################

# Extract target/target loop nucleotides and combine with CPM
single_mismatch_sets.cpm <- OLIGO_INFO %>% 
  filter(oligo_category=='single_mismatch_pairs_saturated') %>% 
  dplyr::select(barcode, oligo_name, target, target_loop) %>% 
  mutate(oligo_name=gsub("_.+", "", oligo_name)) %>% rowwise() %>% 
  mutate(mmpos=mismatchpos(oligo_name)) %>% group_by() %>% 
  mutate(lp1=gsub("(.)..........", "\\1", target_loop)) %>% 
  mutate(lp2=gsub(".(.).........", "\\1", target_loop)) %>% 
  mutate(lp3=gsub("..(.)........", "\\1", target_loop)) %>% 
  mutate(lp4=gsub("...(.).......", "\\1", target_loop)) %>% 
  mutate(lp5=gsub("....(.)......", "\\1", target_loop)) %>% 
  mutate(lp6=gsub(".....(.).....", "\\1", target_loop)) %>% 
  mutate(lp7=gsub("......(.)....", "\\1", target_loop)) %>% 
  mutate(lp10=gsub(".........(.).", "\\1", target_loop)) %>%
  mutate(lp11=gsub("..........(.)", "\\1", target_loop)) %>% 
  mutate(tp1=gsub("(.)..........", "\\1", target)) %>% 
  mutate(tp2=gsub(".(.).........", "\\1", target)) %>% 
  mutate(tp3=gsub("..(.)........", "\\1", target)) %>% 
  mutate(tp4=gsub("...(.).......", "\\1", target)) %>% 
  mutate(tp5=gsub("....(.)......", "\\1", target)) %>% 
  mutate(tp6=gsub(".....(.).....", "\\1", target)) %>% 
  mutate(tp7=gsub("......(.)....", "\\1", target)) %>% 
  mutate(tp10=gsub(".........(.).", "\\1", target)) %>%
  mutate(tp11=gsub("..........(.)", "\\1", target)) %>% 
  dplyr::select(-target, -target_loop) %>% 
  left_join(barcodes.cpm.avg.corrected) %>% 
  mutate(cpm_corrected=ifelse(is.na(cpm_corrected), 0, cpm_corrected))

# Calculate the percentage CPM for each mismatch set/pair
single_mismatch_sets.cpm_perc <- single_mismatch_sets.cpm %>% 
  gather(pos, nuc, -oligo_name, -barcode, -mmpos, -cpm_corrected) %>% 
  mutate(this_pos=as.integer(gsub("lp", "", pos) %>% gsub("tp", "", .))) %>% 
  filter(mmpos==this_pos) %>% dplyr::select(-this_pos) %>% 
  mutate(pos=ifelse(grepl("lp", pos), "target_loop_nuc", "target_nuc")) %>% 
  spread(pos, nuc) %>% arrange(oligo_name, target_loop_nuc, target_nuc) %>% 
  dplyr::select(-barcode) %>% dplyr::rename(cpm=cpm_corrected) %>% 
  group_by(oligo_name, target_loop_nuc) %>% 
  mutate(cpm_perc=cpm/sum(cpm)) %>% 
  dplyr::select(oligo_name, mmpos, target_loop_nuc, target_nuc, cpm, cpm_perc) %>% 
  group_by()

# Determine the top CPM quintile of single mismatch sets
# Filters cases where target loop matches target
# For each single mismatch set, it selects the nucleotides with the highest CPM
# Identifies the top CPM quintile from list
top_quintile_oligo_set <- single_mismatch_sets.cpm_perc %>% 
  filter(target_loop_nuc==target_nuc) %>% 
  group_by(oligo_name) %>% filter(cpm==max(cpm)) %>% group_by() %>% 
  arrange(desc(cpm)) %>% filter(row_number() < 0.2*n()) %>% 
  dplyr::select(oligo_name)

# Averages CPMs across each position in the top quintile
avg_cpm_perc_top_quintile <- single_mismatch_sets.cpm_perc %>% 
  inner_join(top_quintile_oligo_set) %>% 
  group_by(mmpos, target_loop_nuc, target_nuc) %>% 
  summarize(cpm_perc=mean(cpm_perc)*100) %>% group_by() %>% 
  mutate(mmpos=paste0("Position: ", mmpos)) %>% 
  mutate(mmpos=factor(mmpos, levels=unique(mmpos))) %>% 
  mutate(target_loop_nuc=gsub("T", "U", target_loop_nuc))

# Plot
ggplot(avg_cpm_perc_top_quintile %>% 
         dplyr::rename(`% Recovered\nrecombinants`=cpm_perc), 
       aes(x=target_loop_nuc, y=target_nuc, fill=`% Recovered\nrecombinants`)) +
  geom_tile() +
  facet_grid(~mmpos) +
  theme_classic() +
  scale_fill_gradient(low='white', high='red', limits=c(0, 100)) +
  xlab("Target Guide Nucleotide") +
  ylab("Target Nucleotide") +
  coord_equal()
ggsave(paste0(VISUALIZE_DIR, '/single_mismatch_efficiency_grid.pdf'), width=200*(4/3), height=50*(4/3), units='mm')


########################################
# Double Mismatches
########################################

# Get mean CPM for all double mismatch pair positions
double_mismatch_sets.cpm <- OLIGO_INFO %>% 
  filter(oligo_category=='double_mismatch_pairs') %>% 
  filter(target_loop!=target) %>% 
  dplyr::select(oligo_name, barcode) %>% 
  mutate(mm1=as.integer(gsub("^[^_]+_[^_]+_([^_]+)_([^_]+).+", "\\1", oligo_name))) %>% 
  mutate(mm2=as.integer(gsub("^[^_]+_[^_]+_([^_]+)_([^_]+).+", "\\2", oligo_name))) %>% 
  mutate(oligo_name=gsub("^([^_]+_[^_]+).+", "\\1", oligo_name)) %>% 
  left_join(barcodes.cpm.avg.corrected) %>% 
  mutate(cpm_corrected=ifelse(is.na(cpm_corrected), 0, cpm_corrected)) %>% 
  group_by(mm1, mm2) %>% 
  summarize(mean_cpm=mean(cpm_corrected)) %>% group_by()

ggplot() +
  geom_tile(data=double_mismatch_sets.cpm %>% dplyr::mutate(`log10(CPM+1)`=log10(mean_cpm+1)), 
            aes(x=mm1, y=mm2, fill=`log10(CPM+1)`)) + 
  scale_fill_gradient(low='white', high='red') +
  scale_y_continuous(breaks=2:11, expand=c(0, 0)) +
  scale_x_continuous(breaks=1:10, expand=c(0, 0)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = 'light grey')) +
  xlab('Mismatch Position 1') + ylab("Mismatch Position 2") +
  coord_equal()
ggsave(paste0(VISUALIZE_DIR, '/double_mismatch_cpm_grid.pdf'), width=130*(4/3), height=100*(4/3), units='mm')


###########################################
# Top Quintile Enrichment - Perfect Matches
###########################################

# Get the perfect match oligos from 4 different categories
perfect_matches <- OLIGO_INFO %>% filter(oligo_category %in% c(
  "single_mismatch_pairs_other", "single_mismatch_pairs_saturated", "total_mismatch_pairs",
  "double_mismatch_pairs"
)) %>% filter(target==target_loop) %>% dplyr::select(barcode, target) %>% 
  left_join(barcodes.cpm.avg.corrected) %>% 
  mutate(cpm_corrected=ifelse(is.na(cpm_corrected), 0, cpm_corrected)) %>% 
  arrange(desc(cpm_corrected))

# Filter to the top CPM
perfect_matches.top_quintile <- perfect_matches %>% filter(row_number() < 0.2*n())

# Make a motif of top quintile targets
top_quintile_motif <- make_motif(perfect_matches.top_quintile$target)
# Make a motif of all targets
background_motif <- make_motif(perfect_matches$target)

# Create a PWM corrected by background frequencies
motif_pwm <- log2(top_quintile_motif/background_motif)

# Plot
ggseqlogo(motif_pwm, method='custom', seq_type='dna') +
  geom_hline(yintercept=0) +
  ggtitle("Enrichment of Nucleotides in Top CPM Quintile - WT Target: ATCGGGCCTAC") +
  ylab("log2(Enrichment)")
ggsave(paste0(VISUALIZE_DIR, '/enrichment_top_quintile.pdf'), width=120*(4/3), height=50*(4/3), units='mm')

###########################################
# Core Mismatches
###########################################

# Make dataframe of wild-type nucleotides at boundary positions
core_mismatches.cpm <- OLIGO_INFO %>% 
  filter(oligo_category=='core_mismatches') %>% 
  dplyr::select(barcode, oligo_name, target, target_loop) %>% 
  mutate(ltg_core1=ifelse(!grepl("\\*", target_loop), 'C', gsub("^.......(.).+", "\\1", target_loop)),
         ltg_core2=ifelse(!grepl("\\*", target_loop), 'T', gsub("^........(.).+", "\\1", target_loop)),
         rtg_core1=ifelse(!grepl("\\*", target_loop), 'C', gsub("^..........(.).+", "\\1", target_loop)),
         rtg_core2=ifelse(!grepl("\\*", target_loop), 'T', gsub("^...........(.).+", "\\1", target_loop))) %>% 
  left_join(barcodes.cpm.avg.corrected) %>% 
  mutate(cpm_corrected=ifelse(is.na(cpm_corrected), 0, cpm_corrected)) %>% 
  mutate(core=paste0(ltg_core1, ltg_core2, "|", rtg_core1, rtg_core2)) %>% 
  dplyr::select(-ltg_core1, -ltg_core2, -rtg_core1, -rtg_core2) %>% 
  rbind(
    OLIGO_INFO %>% 
      filter(oligo_category=='total_mismatch_pairs') %>% 
      filter(target != target_loop) %>% 
      dplyr::select(barcode, oligo_name, target, target_loop) %>% 
      left_join(barcodes.cpm.avg.corrected) %>% 
      mutate(cpm_corrected=ifelse(is.na(cpm_corrected), 0, cpm_corrected)) %>%
      mutate(core='Neg. control')
  )

core_mismatches.cpm.avg <- core_mismatches.cpm %>% group_by(core) %>% 
  summarize(cpm_corrected=mean(cpm_corrected)) %>% 
  arrange(cpm_corrected)

core_order <- core_mismatches.cpm.avg$core
core_order <- core_order[core_order!='Neg. control']
core_order <- c('Neg. control', core_order)
ggplot(core_mismatches.cpm %>% mutate(core=factor(core, levels=core_order)), 
       aes(x=log10(cpm_corrected+1), y=core)) +
  geom_boxplot(outlier.size=0.1) +
  xlab("log10(CPM+1)") +
  theme_classic() +
  ylab("LTG-Core | RTG-Core") +
  theme(axis.text.y=element_text(size=9, family='mono'))
ggsave(paste0(VISUALIZE_DIR, '/guide_core_mismatch_cpm_distribution.pdf'), width=120, height=175, units='mm')

core_mismatches.cpm <- core_mismatches.cpm %>%
  mutate(core2=ifelse(core=='CT|CT', 'CT|CT', ifelse(
  grepl("CT\\|[^C]T", core), 'CT|DT', ifelse(
    grepl('[^C]T\\|CT', core), 'DT|CT', ifelse(
      grepl('C[^V]\\|CT', core), 'CV|CT', ifelse(
        grepl('CT\\|C[^T]', core), 'CT|CV', ifelse(
          grepl('^CT', core), 'CT|NN', ifelse(
            grepl('CT$', core), 'NN|CT', ifelse(
            grepl('Neg', core), 'Neg. control', 'NN|NN'))))))))) %>%
  mutate(core2=factor(core2, levels=rev(c(
    'CT|CT', 'CT|DT', 'DT|CT', 'CV|CT', 'CT|CV', 'CT|NN', 'NN|CT', 'NN|NN', 
    'Neg. control'))))

ggplot(core_mismatches.cpm %>% mutate(core=factor(core, levels=core_order)), 
       aes(x=log10(cpm_corrected+1), y=core2)) +
  geom_boxplot(outlier.size=0.1) +
  xlab("log10(CPM+1)") +
  theme_classic() +
  ylab("LTG-Core | RTG-Core") +
  theme(axis.text.y=element_text(size=9, family='mono'))
ggsave(paste0(VISUALIZE_DIR, '/guide_core_mismatch_cpm_distribution_grouped.pdf'), width=120, height=100, units='mm')