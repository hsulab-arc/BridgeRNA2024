require(tidyverse)
require(patchwork)


# Get args
args <- commandArgs(trailingOnly = TRUE)

COV_DIR <- args[1]
PLASMID_DIR <- args[2]
ALIGNED_PLASMID_TSV <- args[3]
OUTPDF <- args[4]

mean_pairwise_identity <- function(nuc){

  matches = c()
  for (n1 in nuc){
    for (n2 in nuc){
      if (n1 == n2){
        matches <- c(matches, 1)
      }else{
        matches <- c(matches, 0)
      }
    }
  }
  return(mean(matches))
}


all_annots <- NULL
all_cov <- NULL
for (f in Sys.glob(paste0(COV_DIR, "/*depth.tsv"))){
  libname <- basename(f) %>% gsub("\\..*", "", .)
  covfile <- paste0(COV_DIR, "/", libname, '.read_depth.tsv')
  abundfile <- paste0(COV_DIR, "/", libname, '.read_abundance.tsv')
  annotbed <- read_tsv(paste0(PLASMID_DIR, "/", libname, '.bed'), col_names=F)
  annots <- rbind(
    data.frame(pos=annotbed[1,]$X2:annotbed[1,]$X3, annot=annotbed[1,]$X4),
    data.frame(pos=annotbed[2,]$X2:annotbed[2,]$X3, annot=annotbed[2,]$X4),
    data.frame(pos=annotbed[3,]$X2:annotbed[3,]$X3, annot=annotbed[3,]$X4)
  ) %>% group_by() %>%
    mutate(library=libname) %>%
    dplyr::select(library, pos, annot)
  cov <- read_tsv(covfile) %>%
    mutate(library=libname) %>%
    gather(strand, read_depth, -pos, -annot, -library)
  abund <- read_tsv(abundfile) %>%
    mutate(library=libname) %>%
    left_join(annots %>% dplyr::rename(start=pos, start_annot=annot)) %>%
    left_join(annots %>% dplyr::rename(end=pos, end_annot=annot)) %>%
    dplyr::select(library, strand, start, end, minstart,
                  maxend, start_annot, end_annot, count, seq)

  all_annots <- rbind(all_annots, annots)
  all_cov <- rbind(all_cov, cov)
}

annot_ranges <- all_annots %>% group_by(library, annot) %>%
  summarize(annot_start=min(pos), annot_end=max(pos)) %>%
  group_by()

re_le_ranges <- annot_ranges %>% filter(annot!='Core') %>%
  group_by(library) %>% summarize(
    re_start=min(annot_start), le_end = max(annot_end)
  )

all_cov.prepped <- all_cov %>%
  group_by(library, annot) %>%
  summarize(core_start=min(pos), core_end=max(pos)+1) %>%
  group_by() %>% filter(annot=='Core') %>% dplyr::select(-annot) %>%
  inner_join(all_cov) %>% group_by(library) %>% filter(grepl("fwd", strand)) %>%
  filter(strand=='fwd_unfilt') %>% inner_join(re_le_ranges) %>%
  mutate(plasmid_read_depth_perc=read_depth/max(read_depth)*100) %>%
  mutate(re_le_read_depth_perc=read_depth/max(read_depth[pos >= re_start & pos < le_end])*100) %>%
  mutate(maxdepth=paste("Max. Depth:", max(read_depth))) %>%
  group_by() %>%
  inner_join(read_tsv(ALIGNED_PLASMID_TSV) %>%
               dplyr::rename(library=plasmid, pos=orig_pos)) %>%
  group_by(library) %>%
  mutate(aln_core_start=aln_pos[pos==core_start]) %>%
  mutate(aln_pos=aln_pos-aln_core_start)

# Visualize Coverage
perc_identity <- all_cov.prepped %>% filter(aln_pos>=-20, aln_pos < 195) %>%
  dplyr::select(library, aln_pos, nuc) %>% unique() %>% group_by() %>%
  spread(library, nuc, fill='-') %>%
  gather(library, nuc, -aln_pos) %>%
  group_by(aln_pos) %>%
  summarize(mean_pairwise_ident=mean_pairwise_identity(nuc)*100)

pident <- ggplot(perc_identity, aes(x=aln_pos, y=mean_pairwise_ident)) +
  geom_col(fill='#4d7ae6', width=1.01) + theme_classic() +
  scale_x_continuous(expand=c(0, 0)) +
  xlab("Plasmid Position Relative to Core") +
  ylab("Mean Pairwise Identity")

rele_cov <- ggplot(
  all_cov.prepped %>% filter(aln_pos>=-20, aln_pos < 195),
  aes(x=aln_pos, y=plasmid_read_depth_perc)) +
  geom_area(aes(x=aln_pos, y=plasmid_read_depth_perc), fill='#4d7ae6',
            alpha=0.2, position='identity') +
  geom_line(aes(x=aln_pos, y=plasmid_read_depth_perc), color='#4d7ae6') +
  facet_grid(library~.) +
  scale_x_continuous(expand=c(0, 0)) +
  theme_classic() +
  xlab("Plasmid Position Relative to Core") +
  ylab("% of Maximum Plasmid Read Depth")

layout <- "
A
B
B
B
B
B
B
B
B
B
B
B
"

pident + rele_cov + plot_layout(design = layout)
ggsave(OUTPDF, width=183, height=150, units='mm')