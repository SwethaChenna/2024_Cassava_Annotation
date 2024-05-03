# Based on 2021-11-21

library(GenomicAlignments)
library(rtracklayer)
library(tidyverse)
library(devtools)

scripts <- c("skip_duplicated_reads", "filter_ga_by_terminal_subalignments", "write_grl_as_bed12")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}

bamfiles <- list.files(pattern = "mapq.bam$")

data <- lapply(bamfiles, readGAlignments, use.names = TRUE) %>% lapply(sortSeqlevels) %>% lapply(sort) %>% lapply(skip_duplicated_reads)
names(data) <- bamfiles %>% str_replace("_mapq.bam$", "")

# Filter out unrealistic alignments:
data2 <- lapply(data, filter_ga_by_terminal_subalignments)

# Save as BED12:
data2 <- lapply(data2, grglist)
mapply(write_grl_as_bed12, data2, paste0(names(data2), ".bed"))

