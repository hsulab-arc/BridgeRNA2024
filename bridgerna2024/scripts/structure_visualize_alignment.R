require(tidyverse)
require(ggseqlogo)
require(patchwork)

melt_flankseq_aln <- function(foldaln){
  foldaln_nuc <- foldaln %>%
    select(recid, flankseq_aln_nogaps) %>%
    separate(flankseq_aln_nogaps,
             into=paste0('X', 1:unique(nchar(foldaln$flankseq_aln_nogaps[1]))),
             sep=c(1:unique(nchar(foldaln$flankseq_aln_nogaps[1])))) %>%
    gather(pos, nuc, -recid) %>%
    mutate(pos=as.integer(gsub("X", "", pos)))
  return(foldaln_nuc)
}

melt_linfold_aln <- function(foldaln){
  foldaln_dot <- foldaln %>%
    select(recid, linfold_aln_nogaps) %>%
    separate(linfold_aln_nogaps,
             into=paste0('X', 1:unique(nchar(foldaln$linfold_aln_nogaps[1]))),
             sep=c(1:unique(nchar(foldaln$linfold_aln_nogaps[1])))) %>%
    gather(pos, dot, -recid) %>%
    mutate(pos=as.integer(gsub("X", "", pos)))
  return(foldaln_dot)
}

melt_2struct_aln <- function(foldaln){
  foldaln_2struct <- foldaln %>%
    select(recid, linfold_2struct_aln_nogaps) %>%
    separate(linfold_2struct_aln_nogaps,
             into=paste0('X', 1:unique(nchar(foldaln$linfold_2struct_aln_nogaps[1]))),
             sep=c(1:unique(nchar(foldaln$linfold_2struct_aln_nogaps[1])))) %>%
    gather(pos, struc, -recid) %>%
    mutate(pos=as.integer(gsub("X", "", pos)))
  return(foldaln_2struct)
}

# MAIN
args = commandArgs(trailingOnly=TRUE)

foldalign_path <- args[1]
fiveprime_flank_offset <- as.integer(args[2])
threeprime_flank_offset <- as.integer(args[3])
prefix <- args[4]
outdir <- args[5]

this_tpase_id <- prefix

df <- read_tsv(foldalign_path) %>%
  mutate(left_flankseq_aln=toupper(left_flankseq_aln)) %>%
  mutate(left_flankseq_aln_nogaps=toupper(left_flankseq_aln_nogaps)) %>%
  mutate(right_flankseq_aln=toupper(right_flankseq_aln)) %>%
  mutate(right_flankseq_aln_nogaps=toupper(right_flankseq_aln_nogaps))

if (nrow(df) == 0){
  print('No sequences to visualize')
  write('No sequences to visualize', paste0(outdir, '/', prefix, '.png'))
  quit()
}

flankseqs <- df %>% dplyr::select(recid, left_flankseq_aln_nogaps, right_flankseq_aln_nogaps)
linfolds <- df %>% dplyr::select(recid, left_linfold_aln_nogaps, right_linfold_aln_nogaps)
twostructs <- df %>% dplyr::select(recid, left_linfold_2struct_aln_nogaps, right_linfold_2struct_aln_nogaps)

right_flankseq_melt <- melt_flankseq_aln(flankseqs %>% dplyr::select(recid, flankseq_aln_nogaps=right_flankseq_aln_nogaps)) %>%
  mutate(flankside='RE')
left_flankseq_melt <- melt_flankseq_aln(flankseqs %>% dplyr::select(recid, flankseq_aln_nogaps=left_flankseq_aln_nogaps)) %>%
  mutate(flankside='LE')

combined_flankseq_melt <- rbind(left_flankseq_melt, right_flankseq_melt)

exclude_seqs_gappy <- combined_flankseq_melt %>% group_by(recid) %>%
  summarize(total_missing=sum(nuc=='-'), total=n()) %>%
  filter(total_missing/total > 0.25) %>% dplyr::select(recid)

combined_flankseq_melt <- combined_flankseq_melt %>% anti_join(exclude_seqs_gappy)

right_linfold_melt <- melt_linfold_aln(linfolds %>% dplyr::select(recid, linfold_aln_nogaps=right_linfold_aln_nogaps)) %>%
  mutate(flankside='RE')
left_linfold_melt <- melt_linfold_aln(linfolds %>% dplyr::select(recid, linfold_aln_nogaps=left_linfold_aln_nogaps)) %>%
  mutate(flankside='LE')

combined_linfold_melt <- rbind(right_linfold_melt, left_linfold_melt) %>%
  anti_join(exclude_seqs_gappy)

right_2struct_melt <- melt_2struct_aln(twostructs %>% dplyr::select(recid, linfold_2struct_aln_nogaps=right_linfold_2struct_aln_nogaps)) %>%
  mutate(flankside='RE') %>% inner_join(right_linfold_melt) %>%
  mutate(struc=ifelse(struc=='S' & dot=='(', '5′ Stem',
                      ifelse(struc=='S' & dot==')', '3′ Stem',
                             ifelse(struc=='H', 'Hairpin',
                                    ifelse(struc=='-', 'Gap', 'Other'))))) %>%
  dplyr::select(-dot) %>% anti_join(exclude_seqs_gappy)

left_2struct_melt <- melt_2struct_aln(twostructs %>% dplyr::select(recid, linfold_2struct_aln_nogaps=left_linfold_2struct_aln_nogaps)) %>%
  mutate(flankside='LE') %>% inner_join(left_linfold_melt) %>%
  mutate(struc=ifelse(struc=='S' & dot=='(', '5′ Stem',
                      ifelse(struc=='S' & dot==')', '3′ Stem',
                             ifelse(struc=='H', 'Hairpin',
                                    ifelse(struc=='-', 'Gap', 'Other'))))) %>%
  dplyr::select(-dot) %>% anti_join(exclude_seqs_gappy)


combined_2struct_melt <- rbind(right_2struct_melt, left_2struct_melt)


dot_perc <- combined_linfold_melt %>% anti_join(exclude_seqs_gappy) %>%
  filter(dot!='-') %>%
  group_by(flankside, pos) %>%
  summarize(
    L=sum(dot=='(')/n()*100, R=sum(dot==')')/n()*100, n=n(),
    ) %>% dplyr::select(pos, flankside, L, R, n) %>%
  gather(dot, perc_dot, -pos, -flankside, -n)

struct_perc <- combined_2struct_melt %>% anti_join(exclude_seqs_gappy) %>%
  filter(struc!='Gap') %>%
  group_by(flankside, pos) %>%
  summarize(Hairpin=sum(struc=='Hairpin')/n()*100,
            `5′ Stem`=sum(struc=='5′ Stem')/n()*100,
            `3′ Stem`=sum(struc=='3′ Stem')/n()*100,
            n=n()) %>%
  gather(struc, perc_struc, -pos, -flankside, -n)

outlist <- list()
outlist[['flankseq']] <- combined_flankseq_melt %>%
  mutate(tpase_id=this_tpase_id) %>%
  dplyr::select(tpase_id, recid, flankside, pos, nuc) %>%
  group_by()
outlist[['linfold']] <- combined_linfold_melt %>%
  mutate(tpase_id=this_tpase_id) %>%
  dplyr::select(tpase_id, recid, flankside, pos, dot) %>%
  group_by()
outlist[['2struct']] <- combined_2struct_melt %>%
  mutate(tpase_id=this_tpase_id) %>%
  dplyr::select(tpase_id, recid, flankside, pos, struc) %>%
  group_by()
outlist[['dot_perc']] <- dot_perc %>%
  mutate(tpase_id=this_tpase_id) %>%
  dplyr::select(tpase_id, flankside, pos, dot, perc_dot, n) %>%
  group_by()
outlist[['struct_perc']] <- struct_perc %>%
  mutate(tpase_id=this_tpase_id) %>%
  dplyr::select(tpase_id, flankside, pos, struc, perc_struc, n) %>%
  group_by()

dir.create(outdir, showWarnings=F)

out_flankseq <- paste0(outdir, '/', prefix, '.flankseq.tsv')
out_linfold <- paste0(outdir, '/', prefix, '.linfold.tsv')
out_2struct <- paste0(outdir, '/', prefix, '.2struct.tsv')
out_dot_perc <- paste0(outdir, '/', prefix, '.dot_perc.tsv')
out_struct_perc <- paste0(outdir, '/', prefix, '.struct_perc.tsv')

write_tsv(outlist[['flankseq']], out_flankseq)
write_tsv(outlist[['linfold']], out_linfold)
write_tsv(outlist[['2struct']], out_2struct)
write_tsv(outlist[['dot_perc']], out_dot_perc)
write_tsv(outlist[['struct_perc']], out_struct_perc)

# Plotting

struct_colors <- c('Gap'='white', 'Other'='grey', 'Hairpin'='cyan',
                 '5′ Stem'='green', '3′ Stem'='blue')
struct_colors2 <- c('Hairpin'='cyan', '5′ Stem'='green', '3′ Stem'='blue')

this_flankseq <- outlist[['flankseq']]
this_2struct <- outlist[['2struct']]
this_struct_perc <- outlist[['struct_perc']]

this_flankseqs <- this_flankseq %>% group_by(recid, flankside) %>% arrange(pos) %>%
  summarize(flankseq=paste(nuc, collapse=''))

pident <- read_tsv(foldalign_path) %>% dplyr::select(recid, tpase_pident)

left_cds_start <- this_struct_perc %>% filter(flankside=='LE') %>% .$pos %>% max-fiveprime_flank_offset
right_cds_end <- -1*threeprime_flank_offset
left_cons_struc <- ggplot(this_struct_perc %>% filter(flankside=='LE'),
                        aes(x=pos, y=perc_struc, color=struc, group=struc)) +
    geom_line() + theme_classic() +
    scale_color_manual(values=struct_colors2) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          text=element_text(size=8)) +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 100.0)) +
    geom_hline(yintercept=50.0, linetype='dotted') +
    geom_vline(xintercept=left_cds_start, color='red') +
    ylab("% Stem")

right_cons_struc <- ggplot(this_struct_perc %>% filter(flankside=='RE'),
                         aes(x=pos, y=perc_struc, color=struc, group=struc)) +
    geom_line() + theme_classic() +
    scale_color_manual(values=struct_colors2) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          text=element_text(size=8)) +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, 100.0)) +
    geom_hline(yintercept=50.0, linetype='dotted') +
    geom_vline(xintercept=right_cds_end, color='red') +
    ylab("% Stem")

left_cons_nucs <- ggseqlogo(this_flankseqs %>% filter(flankside=='LE') %>% .$flankseq) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    theme_classic() +
    geom_vline(xintercept=left_cds_start, color='red') +
    theme(axis.title.x=element_blank(),
          text=element_text(size=8)) +
    ggtitle("LE Analysis")

right_cons_nucs <- ggseqlogo(this_flankseqs %>% filter(flankside=='RE') %>% .$flankseq) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    theme_classic() +
    geom_vline(xintercept=right_cds_end, color='red') +
    theme(axis.title.x=element_blank(),
          text=element_text(size=8)) +
    ggtitle("RE Analysis")

left_this_2strict.in <- this_2struct %>% inner_join(pident) %>%
    arrange(tpase_pident) %>%
    mutate(recid=factor(recid, unique(recid))) %>%
    filter(flankside=="LE")
left_this_2strict.in2 <- left_this_2strict.in %>% filter(pos==1)
left_this_2strict.in2 <- left_this_2strict.in2[seq(1, nrow(left_this_2strict.in2), 10), ] %>%
    mutate(tpase_pident=round(tpase_pident))
left_this_2strict.in2 <- left_this_2strict.in2[2:nrow(left_this_2strict.in2),]

right_this_2strict.in <- this_2struct %>% inner_join(pident) %>%
    arrange(tpase_pident) %>%
    mutate(recid=factor(recid, unique(recid))) %>%
    filter(flankside=="RE")
    right_this_2strict.in2 <- right_this_2strict.in %>% filter(pos==max(pos))
    right_this_2strict.in2 <- right_this_2strict.in2[seq(1, nrow(right_this_2strict.in2), 10), ] %>%
    mutate(tpase_pident=round(tpase_pident))
right_this_2strict.in2 <- right_this_2strict.in2[2:nrow(right_this_2strict.in2),]

left_tiled_2struct <- ggplot(left_this_2strict.in,
     aes(x=pos, y=recid, fill=struc)) +
    geom_tile() +
    geom_label(data=left_this_2strict.in2,
               aes(x=pos+5, y=recid, label=tpase_pident),
               size=2.25, fill='white') +
    scale_fill_manual(values=struct_colors) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text=element_text(size=8)) +
    geom_vline(xintercept=left_cds_start, color='red') +
    xlab("LE Position") +
    ylab("Homologous Sequences")

right_tiled_2struct <- ggplot(right_this_2strict.in, aes(x=pos, y=recid, fill=struc)) +
    geom_tile() +
    geom_label(data=right_this_2strict.in2,
               aes(x=pos-5, y=recid, label=tpase_pident),
               size=2.25, fill='white') +
    scale_fill_manual(values=struct_colors) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text=element_text(size=8)) +
    geom_vline(xintercept=right_cds_end, color='red') +
    xlab("RE Position") +
    ylab("Homologous Sequences")

layout <- "
AB
CD
CD
CD
CD
CD
EF
EF
EF
EF
EF
EF
EF
EF
EF
EF
EF
EF
"
mean_pident <- pident %>% .$tpase_pident %>% mean
final_plot <- left_cons_nucs + right_cons_nucs +
    left_cons_struc + right_cons_struc +
    left_tiled_2struct + right_tiled_2struct +
    plot_layout(design = layout) +
    plot_annotation(title = paste0(
        this_tpase_id, '; Avg. Tpase Identity = ', round(mean_pident, 1), '%'
    ))

out_png <- paste0(outdir, '/', prefix, '.png')
out_pdf <- paste0(outdir, '/', prefix, '.pdf')
ggsave(out_png, plot=final_plot, width=600, height=175, units='mm')
ggsave(out_pdf, plot=final_plot, width=600, height=175, units='mm')