require(tidyverse)

args = commandArgs(trailingOnly=TRUE)

indels_path <- args[1]
indels_cov_path <- args[2]
breakpoints_path <- args[3]
indel_bp_overlap_path <- args[4]
insert_overlap_path <- args[5]
out_indels_path <- args[6]
out_breakpoints_path <- args[7]

print(indels_path)
print(indels_cov_path)
print(breakpoints_path)
print(indel_bp_overlap_path)
print(insert_overlap_path)
print(out_indels_path)
print(out_breakpoints_path)

indels <- read_tsv(indels_path)
indels_cov <- read_tsv(indels_cov_path)
breakpoints <- read_tsv(breakpoints_path)


breakpoints <- breakpoints %>% filter(biorep=='biorep1') %>%
  dplyr::select(-biorep) %>%
  dplyr::rename(
    biorep1_strucvar_len=strucvar_len,
    biorep1_breakpoint_position=breakpoint_position,
    biorep1_total=total, biorep1_left_templates=left_templates,
    biorep1_right_templates=right_templates, biorep1_left_freq=left_freq,
    biorep1_right_freq=right_freq, biorep1_freq=freq
  ) %>%
  full_join(
    breakpoints %>% filter(biorep=='biorep2') %>%
      dplyr::select(-biorep) %>%
      dplyr::rename(
        biorep2_strucvar_len=strucvar_len,
        biorep2_breakpoint_position=breakpoint_position,
        biorep2_total=total, biorep2_left_templates=left_templates,
        biorep2_right_templates=right_templates, biorep2_left_freq=left_freq,
        biorep2_right_freq=right_freq, biorep2_freq=freq
      )
  ) %>%
  dplyr::select(sample_id, category, left_merged_breakpoint, right_merged_breakpoint,
                biorep1_strucvar_len, biorep2_strucvar_len,
                biorep1_breakpoint_position, biorep2_breakpoint_position,
                biorep1_total, biorep2_total,
                biorep1_left_templates, biorep2_left_templates,
                biorep1_right_templates, biorep2_right_templates) %>%
  mutate(biorep1_total=ifelse(is.na(biorep1_total), 0, biorep1_total),
         biorep2_total=ifelse(is.na(biorep2_total), 0, biorep2_total)) %>%
  mutate(biorep1_left_freq=ifelse(biorep1_total==0, 0, biorep1_total/biorep1_left_templates),
         biorep1_right_freq=ifelse(biorep1_total==0, 0, biorep1_total/biorep1_right_templates),
         biorep2_left_freq=ifelse(biorep2_total==0, 0, biorep2_total/biorep2_left_templates),
         biorep2_right_freq=ifelse(biorep2_total==0, 0, biorep2_total/biorep2_right_templates)) %>%
  mutate(biorep1_freq=(biorep1_left_freq+biorep1_right_freq)/2,
         biorep2_freq=(biorep2_left_freq+biorep2_right_freq)/2) %>%
  mutate(freq=(biorep1_freq+biorep2_freq)/2) %>%
  arrange(desc(freq))

indel_bp_overlap <- read_tsv(indel_bp_overlap_path)
insert_overlap <- read_tsv(insert_overlap_path)

indels <- indels %>%
  inner_join(indels_cov %>% spread(biorep, coverage) %>%
               dplyr::rename(biorep1_cov=biorep1, biorep2_cov=biorep2)) %>%
  mutate(biorep1_freq=biorep1_total_reads/biorep1_cov,
         biorep2_freq=biorep2_total_reads/biorep2_cov) %>%
  mutate(freq=(biorep1_freq+biorep2_freq)/2) %>%
  left_join(
  indel_bp_overlap %>% group_by(indel) %>%
    summarize(breakpoint=paste(breakpoint, collapse=';')) %>%
    mutate(merged_position=gsub(
      ":[^:]*$", "", indel
    ))) %>% dplyr::select(-indel) %>%
  dplyr::rename(overlapping_breakpoint=breakpoint) %>%
  mutate(avg_read_count=(biorep1_total_reads+biorep2_total_reads)/2) %>%
  arrange(desc(freq))

insert_overlap <- insert_overlap %>% filter(strucvar_type=='breakpoint') %>%
  group_by(strucvar_position) %>%
  summarize(insertions=paste(sort(unique(insertion_name)), collapse=';'))

breakpoints <- breakpoints %>% left_join(
  indel_bp_overlap %>% dplyr::rename(
    left_merged_breakpoint=breakpoint, left_merged_breakpoint_indel_overlap=indel)
) %>% left_join(
  indel_bp_overlap %>% dplyr::rename(
    right_merged_breakpoint=breakpoint, right_merged_breakpoint_indel_overlap=indel)
) %>% left_join(
  insert_overlap %>% dplyr::rename(
    left_merged_breakpoint=strucvar_position, left_merged_breakpoint_insert_overlap=insertions)
) %>% left_join(
  insert_overlap %>% dplyr::rename(
    right_merged_breakpoint=strucvar_position, right_merged_breakpoint_insert_overlap=insertions)
) %>% arrange(desc(freq))

breakpoints %>% write_tsv(out_breakpoints_path)
indels %>% write_tsv(out_indels_path)