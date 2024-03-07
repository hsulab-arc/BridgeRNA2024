require(tidyverse)
require(stringdist)

args <- commandArgs(trailingOnly = TRUE)
insertions_sites_dir <- args[1]
sample_info_path <- args[2]
outdir <- args[3]

sample_info <- read_tsv(sample_info_path) %>% dplyr::select(
  sample=sample_id, exp_target_11mer, exp_donor_11mer, exp_target_14mer, exp_donor_14mer,
)

# Functions
get_kmers <- function(seq, k){
  kmers <- sapply(1:(nchar(seq)-k+1), function(i) substr(seq, i, i+k-1))
  return(kmers)
}

shares_kmer <- function(seq1, seq2, k){
  seq1_kmers <- get_kmers(seq1, k)
  matches <- 0
  for (kmer in seq1_kmers){
    if (grepl(kmer, seq2)){
      matches <- matches +1
    }
  }
  if (matches > 0){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# Load the insertion sites
all_insertions <- NULL
for (f in Sys.glob(paste0(insertions_sites_dir, "/*merged.tsv"))){

  sample_id <- basename(f) %>% gsub(".insertion_sites.merged.tsv", "", .)
  df <- read_tsv(f)
  if (nrow(df) == 0){
    next
  }
  df <- data.frame(sample_id=sample_id) %>%
    mutate(sample=gsub(".biorep.*", "", sample_id)) %>%
    mutate(biorep=gsub(".+(biorep.).*", "\\1", sample_id)) %>%
    cbind(df) %>% group_by()
  all_insertions <- rbind(all_insertions, df)

}

# Get representative core sequences for each insertion window
repr_cores <- all_insertions %>%
  dplyr::count(sample, biorep, contig_id, window_start, window_end, core_start, core_end, strand) %>%
  group_by(sample, biorep, contig_id, window_start, window_end, strand) %>%
  filter(n==max(n)) %>% group_by() %>% dplyr::select(-n)

# Get the insertion counts for each window
window_counts <- all_insertions %>%
  group_by(sample, biorep, contig_id, window_start, window_end, strand) %>%
  summarize(
    lf_count=length(unique(read_id[flankside=='LF'])),
    rf_count=length(unique(read_id[flankside=='RF'])),
    total_count=length(unique(read_id)),
    lf_clipped_count=length(unique(read_id[flankside=='LF' & softclipped_length > 0])),
    rf_clipped_count=length(unique(read_id[flankside=='RF' & softclipped_length > 0])),
    total_clipped_count=length(unique(read_id[softclipped_length > 0]))
    ) %>% group_by() %>%
  mutate(perc_lf_clipped=lf_clipped_count/lf_count*100,
         perc_rf_clipped=rf_clipped_count/rf_count*100,
         perc_total_clipped=total_clipped_count/total_count*100) %>%
  mutate(perc_lf_clipped=ifelse(is.na(perc_lf_clipped), 0, perc_lf_clipped),
         perc_rf_clipped=ifelse(is.na(perc_rf_clipped), 0, perc_rf_clipped),
         perc_total_clipped=ifelse(is.na(perc_total_clipped), 0, perc_total_clipped))

# Final insertion site counts
insertion_site_counts <- all_insertions %>%
  dplyr::select(sample, biorep, contig_id, core_start,
                genome_core, target_11mer, target_14mer) %>%
  unique() %>% inner_join(repr_cores) %>% inner_join(window_counts) %>%
  arrange(desc(total_count)) %>% dplyr::select(-window_start, -window_end)

###################################################
# Handle the samples with non-wildtype/reprogrammed
# Bridge RNAS
###################################################

# Get the insertion counts and compare target with expected
reprogrammed_insertion_counts <- insertion_site_counts %>% filter(!grepl("WT", sample)) %>%
  dplyr::select(
    sample, biorep, contig_id, core_start, core_end, strand, genome_core,
    target_11mer, target_14mer, perc_total_clipped, total_count
  ) %>% mutate(label=paste(
    sample, contig_id, core_start, core_end, strand, genome_core,
    target_11mer, target_14mer, sep='___')) %>%
  dplyr::select(label, biorep, perc_total_clipped, total_count) %>%
  complete(label, biorep, fill=list(perc_total_clipped=0, total_count=0)) %>%
  separate(label, into=c(
    'sample', 'contig_id', 'core_start', 'core_end', 'strand', 'genome_core',
    'target_11mer', 'target_14mer'
  ), sep='___') %>%
  mutate(core_start=as.integer(core_start), core_end=as.integer(core_end)) %>%
  inner_join(sample_info) %>% rowwise() %>%
  mutate(target_lev_dist=stringdist(target_11mer, exp_target_11mer, method = "lv"),
         donor_lev_dist=stringdist(target_11mer, exp_donor_11mer, method = "lv")) %>%
  group_by() %>% group_by(sample, biorep) %>%
  mutate(perc_insertions=total_count/sum(total_count)*100) %>% group_by() %>%
  arrange(desc(perc_insertions))

# Merge data across bioreps and calculate insertion percentage
reprogrammed_insertion_perc.prefilter <- reprogrammed_insertion_counts %>%
  group_by(sample, contig_id, core_start, core_end, strand, genome_core,
           target_11mer, target_14mer, target_lev_dist, donor_lev_dist) %>%
  summarize(perc_insertions=mean(perc_insertions),
            perc_clipped=mean(perc_total_clipped),
            min_count=min(total_count),
            max_count=max(total_count),
            avg_count=mean(total_count),
            n=n()) %>% arrange(desc(perc_insertions)) %>%
  group_by()

# Keep sites where the % insertions exceeds 1%, the edit distance from the
# expected target is < 3, or the edit distance from the expected donor is < 3
keep_filt1 <- reprogrammed_insertion_perc.prefilter %>%
  filter(perc_insertions > 1 | target_lev_dist < 3 | donor_lev_dist < 3)

# Exclude remaining sites where the % clipped is >= 25%
exclude_filt1 <- reprogrammed_insertion_perc.prefilter %>% anti_join(keep_filt1) %>%
  filter(perc_clipped >= 25)

# Calculate the maximum insertion % across all samples for each target.
# Only keep sites where the max_perc_insertions > 1
max_insertion_rate_gt1 <- reprogrammed_insertion_perc.prefilter %>%
  group_by(contig_id, core_start, core_end) %>%
  summarize(max_perc_insertions=max(perc_insertions)) %>%
  filter(max_perc_insertions > 1) %>%
  arrange(max_perc_insertions) %>% group_by()

# Calculate the minimum edit distance from the expected target across all
# samples for each target. Only keep sites where the min_target_lev_dist < 3
min_target_lev_lt3 <- reprogrammed_insertion_perc.prefilter %>%
  group_by(contig_id, core_start, core_end) %>%
  summarize(min_target_lev_dist=min(target_lev_dist)) %>%
  filter(min_target_lev_dist < 3) %>% group_by()

# Exclude remaining sites that occur at high percentage in other samples and
# closely match an expected target in other samples.
exclude_filt2 <- reprogrammed_insertion_perc.prefilter %>% anti_join(keep_filt1) %>%
  anti_join(exclude_filt1) %>% inner_join(max_insertion_rate_gt1) %>%
  inner_join(min_target_lev_lt3)

# Get the final list of sites to keep after applying all filters
reprogrammed_insertion_perc.final_keep <- keep_filt1 %>% rbind(
  reprogrammed_insertion_perc.prefilter %>% anti_join(keep_filt1) %>%
    anti_join(exclude_filt1) %>% anti_join(exclude_filt2)
) %>% dplyr::select(sample, contig_id, core_start, core_end, strand)

# Recalculate the insertion $ for each sample after applying the filters
reprogrammed_insertions.final <- reprogrammed_insertion_counts %>%
  inner_join(reprogrammed_insertion_perc.final_keep) %>%
  dplyr::select(-perc_total_clipped, -perc_insertions) %>%
  spread(biorep, total_count) %>%
  dplyr::rename(biorep1_insertion_count=biorep1,
                biorep2_insertion_count=biorep2) %>%
  group_by(sample) %>%
  mutate(biorep1_insertion_perc=biorep1_insertion_count/sum(biorep1_insertion_count)*100,
         biorep2_insertion_perc=biorep2_insertion_count/sum(biorep2_insertion_count)*100) %>%
  group_by() %>%
  mutate(avg_insertion_perc=(biorep1_insertion_perc+biorep2_insertion_perc)/2) %>%
  arrange(desc(avg_insertion_perc)) %>%
  mutate(bridgerna_id=gsub(".*bridgeRNA-(.).*", "Bridge RNA \\1", sample)) %>%
  mutate(rtg_extension=ifelse(grepl("RTG4", sample), "IS621\nRTG", "Extended\nRTG")) %>%
  mutate(rtg_extension=factor(rtg_extension, levels=c("IS621\nRTG", "Extended\nRTG")))

# Group insertion % by target distance from expected target
reprogrammed_insertions.target_dist_bin <- reprogrammed_insertions.final %>%
  mutate(target_lev_dist_bin=ifelse(target_lev_dist==0, "0", ifelse(
    target_lev_dist==1, "1", ifelse(
      target_lev_dist==2, "2", ">2")))) %>%
  mutate(target_lev_dist_bin=factor(target_lev_dist_bin, levels=c("0", "1", "2", ">2"))) %>%
  group_by(sample, bridgerna_id, rtg_extension, target_lev_dist_bin) %>%
  summarize(num_sites=n(), biorep1_insertion_count=sum(biorep1_insertion_count),
            biorep2_insertion_count=sum(biorep2_insertion_count)) %>%
  mutate(biorep1_insertion_perc=biorep1_insertion_count/sum(biorep1_insertion_count)*100,
         biorep2_insertion_perc=biorep2_insertion_count/sum(biorep2_insertion_count)*100) %>%
  group_by() %>%
  mutate(avg_insertion_perc=(biorep1_insertion_perc+biorep2_insertion_perc)/2)

# Prepare data for plotting
df.avg_insertion_perc <- reprogrammed_insertions.target_dist_bin %>%
  dplyr::select(bridgerna_id, rtg_extension, target_lev_dist_bin, num_sites, avg_insertion_perc) %>%
  complete(bridgerna_id, rtg_extension, target_lev_dist_bin,
           fill=list(avg_insertion_perc=0, num_sites=0))
df.biorep_insertion_perc <- reprogrammed_insertions.target_dist_bin %>%
  dplyr::select(bridgerna_id, rtg_extension, target_lev_dist_bin, biorep1_insertion_perc, biorep2_insertion_perc) %>%
  complete(bridgerna_id, rtg_extension, target_lev_dist_bin,
           fill=list(biorep1_insertion_perc=0, biorep2_insertion_perc=0)) %>%
  gather(metric, val, -bridgerna_id, -rtg_extension, -target_lev_dist_bin)

# Plot insertion % by distance from expected target
p <- ggplot(df.avg_insertion_perc) +
  geom_col(aes(x=rtg_extension, y=avg_insertion_perc, fill=target_lev_dist_bin),
           position='dodge', color='black') +
  geom_text(aes(x=rtg_extension, y=avg_insertion_perc+5,
                label=num_sites, group=target_lev_dist_bin),
            position=position_dodge(width=0.9)) +
  geom_point(data=df.biorep_insertion_perc, aes(
    x=rtg_extension, y=val, group=target_lev_dist_bin
  ), position=position_dodge(width=0.9), size=0.5) +
  facet_wrap(~bridgerna_id, scales='free') +
  scale_y_continuous(
    limits=c(0, 100.0), expand=c(0, 0),
    breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
    ) +
  theme_classic()  +
  scale_fill_manual(values=c("#92ff92", "#fffc88", '#ffa500', '#bebebe'),
                    name='Lev. Distance\nfrom Target') +
  ylab("Insertion Reads (%)") +
  theme(axis.title.x=element_blank())
ggsave(paste0(outdir, "/programmed_insertion_perc_by_target_dist.pdf"),
       plot=p, width=200, height=250, units='mm')


# Prepare data for plotting
insertion_ranks.df <- reprogrammed_insertions.final %>%
  mutate(bridgerna_id=gsub(".*bridgeRNA-(.).*", "Bridge RNA \\1", sample)) %>%
  group_by(sample) %>% mutate(insertion_rank=row_number()) %>%
  mutate(target_lev_dist_bin=ifelse(target_lev_dist==0, "0", ifelse(
    target_lev_dist==1, "1", ifelse(
      target_lev_dist==2, "2", ">2")))) %>%
  mutate(target_lev_dist_bin=factor(target_lev_dist_bin, levels=c("0", "1", "2", ">2"))) %>%
  group_by()

top_sites <- insertion_ranks.df %>% filter(insertion_rank <= 3) %>%
  dplyr::select(bridgerna_id, rtg_extension, insertion_rank, avg_insertion_perc,
                target_14mer)

#Plot data by ranked insertion site specificity
p <- ggplot(insertion_ranks.df, aes(x=insertion_rank, y=avg_insertion_perc)) +
  geom_point(aes(fill=target_lev_dist_bin),
             colour="black",pch=21) +
  geom_text(data=top_sites, aes(label=target_14mer), hjust=0, size=2, nudge_x=1) +
  facet_grid(bridgerna_id~rtg_extension) +
  theme_classic() +
  theme(axis.line=element_line()) +
  scale_fill_manual(values=c("#92ff92", "#fffc88", '#ffa500', '#bebebe'),
                     name='Lev. Distance\nfrom Target') +
  scale_y_continuous(
    limits=c(0, 100.0), breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  ) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.y=element_blank()) +
  xlab("Genomic Targets Ranked by Specificity")
ggsave(paste0(outdir, "/programmed_insertion_perc_specificity.pdf"),
       plot=p, width=200, height=290, units='mm')

# Categorize the off-targets
ontargets <- insertion_ranks.df %>% filter(target_lev_dist==0) %>%
  mutate(offtarget_category='On-Target') %>%
  dplyr::select(sample, contig_id, core_start, core_end, strand, genome_core,
                target_14mer, offtarget_category)
targetlike_offtargets <- insertion_ranks.df %>% anti_join(ontargets) %>%
  filter(target_lev_dist < 3) %>%
  mutate(offtarget_category='Target-like, Lev. Dist < 3') %>%
  dplyr::select(sample, contig_id, core_start, core_end, strand, genome_core,
                target_14mer, offtarget_category)
donorlike_offtargets <- insertion_ranks.df %>% anti_join(ontargets) %>%
  anti_join(targetlike_offtargets) %>%
  filter(donor_lev_dist < 3) %>%
  mutate(offtarget_category='Donor-like, Lev. Dist < 3') %>%
  dplyr::select(sample, contig_id, core_start, core_end, strand, genome_core,
                target_14mer, offtarget_category)

shares_kmers <- insertion_ranks.df %>% anti_join(ontargets) %>%
  anti_join(targetlike_offtargets) %>%
  anti_join(donorlike_offtargets) %>% rowwise() %>%
  mutate(shares_target_kmer=ifelse(
    shares_kmer(exp_target_14mer, target_14mer, k=9), 9, ifelse(
      shares_kmer(exp_target_14mer, target_14mer, k=8), 8, ifelse(
        shares_kmer(exp_target_14mer, target_14mer, k=7), 7, ifelse(
          shares_kmer(exp_target_14mer, target_14mer, k=6), 6, ifelse(
            shares_kmer(exp_target_14mer, target_14mer, k=5), 5, ifelse(
              shares_kmer(exp_target_14mer, target_14mer, k=4), 4, 0
              ))))))) %>%
  mutate(shares_donor_kmer=ifelse(
    shares_kmer(exp_donor_14mer, target_14mer, k=9), 9, ifelse(
      shares_kmer(exp_donor_14mer, target_14mer, k=8), 8, ifelse(
        shares_kmer(exp_donor_14mer, target_14mer, k=7), 7, ifelse(
          shares_kmer(exp_donor_14mer, target_14mer, k=6), 6, ifelse(
            shares_kmer(exp_donor_14mer, target_14mer, k=5), 5, ifelse(
              shares_kmer(exp_donor_14mer, target_14mer, k=4), 4, 0
            ))))))) %>% group_by() %>%
  filter(shares_target_kmer >= 4 | shares_donor_kmer >= 4)

shares_target_kmers <- shares_kmers %>%
  filter((shares_target_kmer > shares_donor_kmer) | (shares_target_kmer == shares_donor_kmer & target_lev_dist > donor_lev_dist)) %>%
  mutate(offtarget_category=paste0('Target-like, Shares ', shares_target_kmer, "-mer")) %>%
  dplyr::select(sample, contig_id, core_start, core_end, strand, genome_core,
                target_14mer, offtarget_category)

shares_donor_kmers <- shares_kmers %>%
  filter((shares_target_kmer < shares_donor_kmer) | (shares_target_kmer == shares_donor_kmer & target_lev_dist < donor_lev_dist)) %>%
  mutate(offtarget_category=paste0('Donor-like, Shares ', shares_donor_kmer, "-mer")) %>%
  dplyr::select(sample, contig_id, core_start, core_end, strand, genome_core,
                target_14mer, offtarget_category)

offtarget_categories <- rbind(
  ontargets,
  targetlike_offtargets,
  donorlike_offtargets,
  shares_target_kmers,
  shares_donor_kmers
)

######################################
# Prepare WT Insertion Data
######################################
wt_insertion_counts <- insertion_site_counts %>% filter(grepl("WT", sample)) %>%
  inner_join(sample_info) %>%
  dplyr::select(
    sample, biorep, contig_id, core_start, core_end, strand, genome_core,
    target_11mer, target_14mer, exp_target_14mer, exp_donor_11mer,
    perc_total_clipped, total_count
  ) %>% mutate(label=paste(
    sample, contig_id, core_start, core_end, strand, genome_core,
    target_11mer, target_14mer, sep='___')) %>%
  dplyr::select(label, biorep, perc_total_clipped, total_count) %>%
  complete(label, biorep, fill=list(perc_total_clipped=0, total_count=0)) %>%
  separate(label, into=c(
    'sample', 'contig_id', 'core_start', 'core_end', 'strand', 'genome_core',
    'target_11mer', 'target_14mer'
  ), sep='___') %>%
  mutate(core_start=as.integer(core_start), core_end=as.integer(core_end)) %>%
  inner_join(sample_info) %>% rowwise() %>%
  mutate(target_lev_dist=stringdist(target_11mer, exp_target_11mer, method = "lv"),
         donor_lev_dist=stringdist(target_11mer, exp_donor_11mer, method = "lv")) %>%
  group_by() %>%
  group_by(sample, biorep) %>%
  mutate(perc_insertions=total_count/sum(total_count)*100) %>% group_by() %>%
  arrange(desc(perc_insertions)) %>% group_by()

# Get all % insertions across both bioreps
combined_perc_insertions <- reprogrammed_insertion_counts %>%
  inner_join(reprogrammed_insertion_perc.final_keep) %>%
  group_by(sample, biorep) %>%
  mutate(perc_insertions=total_count/sum(total_count)*100) %>% group_by() %>%
  rbind(wt_insertion_counts)

# Prepare for plotting
combined_perc_insertions.df <- combined_perc_insertions %>%
  dplyr::select(sample, biorep, contig_id, core_start, core_end,
                strand, perc_insertions) %>%
  spread(biorep, perc_insertions) %>%
  mutate(bridgerna_id=gsub(".*bridgeRNA-(.).*", "Bridge RNA \\1", sample)) %>%
  mutate(bridgerna_id=ifelse(bridgerna_id=='Bridge RNA W', 'Bridge RNA WT', bridgerna_id)) %>%
  mutate(rtg_extension=ifelse(grepl("RTG4", sample), "IS621 RTG", "Extended RTG")) %>%
  mutate(rtg_extension=factor(rtg_extension, levels=c("IS621 RTG", "Extended RTG")))
# Calculate correlations for % insertions across bioreps
correls <- combined_perc_insertions.df %>%
  group_by(bridgerna_id, rtg_extension) %>%
  summarize(correl=cor(biorep1, biorep2)) %>%
  mutate(correl=paste("r =", round(correl, 4)))
# Plot
p <- ggplot(combined_perc_insertions.df,
       aes(x=biorep1, y=biorep2)) +
  geom_abline(linetype='dashed', color='blue') +
  geom_point(alpha=0.5) + facet_wrap(~bridgerna_id+rtg_extension, ncol=4) +
  geom_text(data=correls, aes(x=10, y=90, label=correl), hjust=0, size=3) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 100.0)) +
  scale_x_continuous(limits=c(0, 100.0)) +
  coord_equal() +
  xlab("Repl. 1 - % insertion reads at site") +
  ylab("Repl. 2 - % insertion reads at site")
ggsave(paste0(outdir, "/insertion_perc_replicate_comparison.pdf"),
       plot=p, width=300, height=250, units='mm')


# Plot WT specificity
wt_insertion_perc <- wt_insertion_counts %>% dplyr::select(
  sample, biorep, contig_id, core_start, core_end, strand, genome_core,
  target_lev_dist, donor_lev_dist, target_11mer, target_14mer, total_count
) %>% spread(biorep, total_count) %>%
  dplyr::rename(biorep1_insertion_count=biorep1,
                biorep2_insertion_count=biorep2) %>%
  mutate(biorep1_insertion_perc=biorep1_insertion_count/sum(biorep1_insertion_count)*100,
         biorep2_insertion_perc=biorep2_insertion_count/sum(biorep2_insertion_count)*100) %>%
  mutate(avg_insertion_perc=(biorep1_insertion_perc+biorep2_insertion_perc)/2) %>%
  arrange(desc(avg_insertion_perc)) %>%
  mutate(insertion_rank=row_number()) %>%
  mutate(target_lev_dist_bin=ifelse(
    target_lev_dist==0 | target_11mer=='ATCAGGCCTAC', "0 or WT Target", ifelse(
    target_lev_dist==1, "1", ifelse(
      target_lev_dist==2, "2", ">2")))) %>%
  mutate(target_lev_dist_bin=factor(target_lev_dist_bin,
    levels=c("0 or WT Target", "1", "2", ">2")
    )) %>%
  group_by() %>%
  mutate(wt_rtg_ext3=substr(target_14mer, 12, 14) == "GCA")

#Plot WT data by ranked insertion site specificity
p <- ggplot(wt_insertion_perc, aes(x=insertion_rank, y=avg_insertion_perc)) +
  geom_point(aes(fill=target_lev_dist_bin),
             colour="black",pch=21) +
  geom_point(data=wt_insertion_perc %>% filter(wt_rtg_ext3==T), color='red') +
  scale_fill_manual(values=c("#92ff92", "#fffc88", '#ffa500', '#bebebe'),
                    name='Lev. Distance\nfrom Target') +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.y=element_blank()) +
  ylab("Insertion Reads (%)") +
  xlab("Genomic Targets Ranked by Specificity") +
  theme_classic()
ggsave(paste0(outdir, "/wt_insertion_perc_specificity.pdf"),
       plot=p, width=200, height=175, units='mm')

######################################
# Write final tables
######################################
programmed_insertion_sites.df <- insertion_ranks.df %>% inner_join(offtarget_categories) %>%
  dplyr::select(
    bridgerna_id, rtg_extension, contig_id, core_start, core_end, strand,
    genome_core, target_11mer, target_14mer, exp_target_11mer, exp_target_14mer,
    exp_donor_11mer, exp_donor_14mer, target_lev_dist, donor_lev_dist,
    offtarget_category, biorep1_insertion_count, biorep2_insertion_count,
    avg_insertion_perc, insertion_rank
    ) %>%
  mutate(rtg_extension=as.character(rtg_extension) %>% gsub("\n", " ", .))

programmed_insertion_sites.df %>% write_tsv(
  paste0(outdir, "/programmed_insertion_sites.tsv")
)

wt_insertion_sites.df <- wt_insertion_perc %>%
  inner_join(sample_info) %>%
  dplyr::select(
    contig_id, core_start, core_end, strand,
    genome_core, target_11mer, target_14mer, exp_target_11mer, exp_target_14mer,
    exp_donor_11mer, exp_donor_14mer, target_lev_dist, donor_lev_dist,
    biorep1_insertion_count, biorep2_insertion_count, avg_insertion_perc, insertion_rank
      )
wt_insertion_sites.df %>% write_tsv(
  paste0(outdir, "/wt_insertion_sites.tsv")
)