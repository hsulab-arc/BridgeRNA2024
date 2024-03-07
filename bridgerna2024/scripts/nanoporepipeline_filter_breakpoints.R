require(tidyverse)

args = commandArgs(trailingOnly=TRUE)

df <- read_tsv(args[1])

df.counts <- df %>%
  mutate(strucvar_len=right_min_pos-(left_min_pos-1)) %>%
  mutate(breakpoint_position=paste0(
    left_contig, ":", left_min_pos, ":", left_strand, "|", right_contig, ":",
    right_min_pos, ":", right_strand)) %>%
  group_by(sample_id, biorep, category, left_strand, right_strand,
           left_merged_breakpoint, right_merged_breakpoint) %>%
  summarize(strucvar_len = strucvar_len[total==max(total)][1],
            breakpoint_position=breakpoint_position[total==max(total)][1],
            total=sum(total),
            left_templates=max(left_templates),
            right_templates=max(right_templates)) %>%
  group_by() %>%
  mutate(left_freq=total/left_templates, right_freq=total/right_templates)

exclude_untransfected <- df.counts %>% filter(grepl("Untransfected", sample_id)) %>%
  dplyr::select(left_merged_breakpoint, right_merged_breakpoint) %>%
  unique()

keep <- df.counts %>%
  anti_join(exclude_untransfected) %>%
  group_by(sample_id, left_merged_breakpoint, right_merged_breakpoint) %>%
  summarize(total=sum(total)) %>%
  filter(total > 1) %>% group_by() %>%
  dplyr::select(left_merged_breakpoint, right_merged_breakpoint) %>%
  unique()

filtered_breakpoints <- df.counts %>% inner_join(keep) %>%
  mutate(freq=(left_freq+right_freq)/2)

write_tsv(filtered_breakpoints, args[2])