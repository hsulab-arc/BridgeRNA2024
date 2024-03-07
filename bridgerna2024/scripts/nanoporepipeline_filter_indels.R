require(tidyverse)

args = commandArgs(trailingOnly=TRUE)

df <- read_tsv(args[1])

df.counts <- df %>%
  dplyr::select(sample_id, biorep, indel_type, ref_name, ref_pos, merged_position, seq_cluster) %>%
  mutate(indel_position=paste0(ref_name, ":", ref_pos)) %>%
  group_by(sample_id, biorep, indel_type, indel_position, merged_position, seq_cluster) %>%
  summarize(total=n()) %>%
  group_by(sample_id, indel_type, merged_position, seq_cluster) %>%
  mutate(indel_position=indel_position[total==max(total)][1]) %>%
  group_by(sample_id, biorep, indel_type, merged_position, seq_cluster) %>%
  summarize(indel_position=indel_position[total==max(total)][1], total=sum(total)) %>%
  group_by()

exclude_untransfected <- df.counts %>%
  filter(grepl("Untransfected", sample_id)) %>%
  dplyr::select(indel_type, merged_position, seq_cluster) %>%
  unique()

keep <- df.counts %>%
  anti_join(exclude_untransfected) %>%
  group_by(sample_id, indel_type, merged_position, seq_cluster) %>%
  summarize(total=sum(total)) %>%
  filter(total > 1) %>% group_by() %>%
  dplyr::select(indel_type, merged_position, seq_cluster) %>%
  unique()

repseqs <- df %>% dplyr::select(seq_cluster, seq) %>% unique() %>%
  dplyr::count(seq_cluster, seq) %>%
  group_by(seq_cluster) %>% mutate(medlen=median(nchar(seq))) %>%
  filter(n==max(n)) %>% summarize(
    medseqlen=unique(medlen),
    repseq=seq[abs(nchar(seq)-medlen)==min(abs(nchar(seq)-medlen))][1]
    )

filtered_indels <- df.counts %>%
  inner_join(keep) %>%
  arrange(desc(total)) %>%
  group_by() %>% inner_join(repseqs) %>%
  spread(biorep, total, fill=0) %>%
  dplyr::rename(biorep1_total_reads=biorep1,
                biorep2_total_reads=biorep2) %>%
  dplyr::select(sample_id, indel_type, merged_position, seq_cluster, indel_position, medseqlen,
                biorep1_total_reads, biorep2_total_reads, repseq)


df.out <- df %>% inner_join(keep)

filtered_indels %>% write_tsv(args[2])
df.out %>% write_tsv(args[3])