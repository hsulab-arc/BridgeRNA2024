require(tidyverse)
require(ggseqlogo)
require(patchwork)

args = commandArgs(trailingOnly=TRUE)

BARCODE_DIR <- args[1]
OLIGO_INFO <- read_tsv(args[2])
OLIGO_COVERAGE_DIR <- args[3]

all_barcode_counts <- NULL
for (f in Sys.glob(paste0(BARCODE_DIR, "/*mapped.tsv"))){
  libname <- strsplit(basename(f), '\\.')[[1]][1]
  counts <- read_tsv(f) %>% dplyr::count(mapped_barcode) %>% 
    mutate(library=libname) %>% 
    dplyr::select(library, barcode=mapped_barcode, barcode_count=n)
  all_barcode_counts <- rbind(all_barcode_counts, counts)
}

all_barcode_counts <- all_barcode_counts %>% 
  mutate(library=factor(library, levels=unique(library))) %>% 
  mutate(barcode=factor(barcode, levels=OLIGO_INFO$barcode)) %>% 
  complete(library, barcode, fill=list(barcode_count=0))

all_barcode_counts %>% 
  write_tsv(paste0(OLIGO_COVERAGE_DIR, "/all_oligo_coverage.tsv"))