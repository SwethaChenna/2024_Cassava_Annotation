library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)
library(collections)
library(TranscriptomeReconstructoR)
library(devtools)

devtools::source_url("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/batch_read_track_data.R?raw=TRUE")


write_grl_as_bed12 <- function(grl, filename, name_in_mcols = FALSE, zero_based = TRUE, default_name = ".", default_score = ".", default_itemRgb = ".") {
  # zero_based = FALSE is just for backwards compatibility with the code written before 2022-11-08;
  # From now on, always use write_grl_as_bed12(zero_based = TRUE)!
  # default_score = "." and default_itemRgb = "." are OK for IGV browser;
  # However, UCSC wants numeric values, e.g. default_score = "0" and default_itemRgb = "0,0,0"
  gr <- range(grl) %>% unlist(use.names = FALSE)
  grl_unl <- unlist(grl, use.names = FALSE)
  idx_fw <- lapply(elementNROWS(grl), function(x) { seq(1, x) }) %>% unlist() %>% unname()
  idx_rev <- lapply(elementNROWS(grl), function(x) { seq(x, 1) }) %>% unlist() %>% unname()
  idx <- ifelse(strand(grl_unl) == "+", idx_fw, idx_rev)
  mc <- grl_unl[idx == 1] %>% mcols()
  chrom <- seqnames(gr) %>% as.character()
  chromStart <- start(gr)
  chromEnd <- end(gr)
  if (isTRUE(name_in_mcols)) {
    if (is.null(mc$name)) {
      name <- default_name
    } else {
      name <- mc$name
    }
  } else {
    if (is.null(names(grl))) {
      name <- default_name
    } else {
      name <- names(grl)
    }
  }
  if (is.null(mc$score)) {
    score <- default_score
  } else {
    score <- mc$score
  }
  strand <- strand(gr) %>% as.character()
  if (is.null(mc$thick)) {
    thickStart <- chromStart # already zero-offset
    thickEnd <- chromEnd
  } else {
    thickStart <- start(mc$thick)
    thickEnd <- end(mc$thick)
  }
  if (is.null(mc$itemRgb)) {
    itemRgb <- default_itemRgb
  } else {
    itemRgb <- mc$itemRgb
  }
  blockCount <- elementNROWS(grl)
  blockSizes <- width(grl) %>% lapply(str_c, collapse = ",") %>% unlist()
  blockStarts <- start(grl) %>% `-`(start(gr)) %>% lapply(str_c, collapse = ",") %>% unlist()
  if (isTRUE(zero_based)) {
    chromStart <- chromStart - 1
    thickStart <- thickStart - 1
  }
  tbl <- tibble(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts)
  tbl <- arrange(tbl, chrom, chromStart, chromEnd)
  write_tsv(tbl, filename, col_names = FALSE)
}



# Prepare Seqinfo object ------------------------------------------------------

si <- file.path("I:/SCIENCE-PLEN-Marquardt_lab/Swetha/Genomes/Cassava/v6.54/", "Manihot_esculenta_v6.54.txt") %>% read_tsv(col_names = c("chr", "len"), col_types = "ci")
si <- Seqinfo(seqnames = si$chr, seqlengths = si$len)

##### Prepare input data -----------------------------------------------------
# (TrRR needs BAM files but only Bedgraph and BED12 were available)

# ONT data from BED12 files:
read_grl_from_bed12 <- function(bedfile, seqinfo, score = FALSE) {
  message(basename(bedfile)); flush.console()
  tbl <- read_tsv(bedfile, 
                  col_names = c("chr", "start", "end", "name", "score", "strand", "tstart", "tend", "rgb", "bcount", "bwidth", "bstart"),
                  col_types = "ciiciciicicc", na = ".")
  tbl2 <- tbl %>% dplyr::select(-tstart, -tend, -rgb, -bcount) %>% separate_rows(bwidth, bstart, sep = ",") %>% 
    mutate(bwidth = as.integer(bwidth), bstart = as.integer(bstart), nstart = start + bstart)
  gr <- GRanges(seqnames = tbl2$chr, IRanges(start = tbl2$nstart, end = tbl2$nstart + tbl2$bwidth - 1), strand = tbl2$strand, seqinfo = seqinfo)
  if (isTRUE(score)) {
    score(gr) <- tbl2$score
  }
  grl <- split(gr, tbl2$name) %>% sort()
  return(grl)
}

ont_dir <- "I:/SCIENCE-PLEN-Marquardt_lab/The_Holy_Grail/Cassava/Our_ONT_Direct_RNAseq/BED12"
ont_cold <- file.path(ont_dir, "Cold3h.bed.gz") %>% read_grl_from_bed12(seqinfo = si)
ont_wt <- file.path(ont_dir, "TM12.bed.gz") %>% read_grl_from_bed12(seqinfo = si)
saveRDS(ont_cold, "ONT_cold.RDS")
saveRDS(ont_wt, "ONT_wt.RDS")

# STRIPE-seq and Lexogen data from Bedgraph files:
tss_dir <- "I:/SCIENCE-PLEN-Marquardt_lab/The_Holy_Grail/Cassava/Our_STRIPEseq/Bedgraph/Raw"
tss_files <- list.files(tss_dir, pattern = "rep.*bedgraph.gz$")
tss_data <- batch_read_track_data(tss_files, tss_dir, "bedGraph", seqinfo = si)
names(tss_data) <- names(tss_data) %>% str_replace("_fw_rev.bedgraph.gz", "")
saveRDS(tss_data, "tss_data.RDS")

pas_dir <- "I:/SCIENCE-PLEN-Marquardt_lab/The_Holy_Grail/Cassava/Our_Lexogen_QuantSeq_REV/Bedgraph/Raw"
pas_files <- list.files(pas_dir, pattern = "rep.*bedgraph.gz$")
pas_data <- batch_read_track_data(pas_files, pas_dir, "bedGraph", seqinfo = si)
names(pas_data) <- names(pas_data) %>% str_replace("_fw_rev.bedgraph.gz", "")
saveRDS(pas_data, "pas_data.RDS")

# ssRNA-seq data (as a substitute for plaNET-seq) from Bedgraph files:
ssrna_dir <- "I:/SCIENCE-PLEN-Marquardt_lab/The_Holy_Grail/Cassava/Our_ssRNA-seq/Bedgraph/Raw"
ssrna_files <- list.files(ssrna_dir, pattern = "rep.*bedgraph.gz$")
ssrna_data <- batch_read_track_data(ssrna_files, ssrna_dir, "bedGraph", seqinfo = si)
names(ssrna_data) <- names(ssrna_data) %>% str_replace("_raw.bedgraph.gz$", "") %>% str_replace("^sR\\d_", "")
# Merge replicates:
ssrna_wt_data <- TranscriptomeReconstructoR:::merge_GRanges(ssrna_data[1:4])
ssrna_cold_data <- TranscriptomeReconstructoR:::merge_GRanges(ssrna_data[5:8])
saveRDS(ssrna_wt_data, "ssrna_wt_data.RDS")
saveRDS(ssrna_cold_data, "ssrna_cold_data.RDS")


##### Load the prepared data ----------------------------------------

ont_wt <- readRDS("ONT_wt.RDS")
ont_cold <- readRDS("ONT_cold.RDS")
tss_data <- readRDS("tss_data.RDS")
tss_cold <- tss_data[1:4]
tss_wt <- tss_data[5:8]
pas_data <- readRDS("pas_data.RDS")
pas_cold <- pas_data[1:4]
pas_wt <- pas_data[5:8]
rna_wt <- readRDS("ssrna_wt_data.RDS")
rna_cold <- readRDS("ssrna_cold_data.RDS")

# ##### Call TSS and PAS ----------------------------------------------
# 
# tss_tc_wt <- call_TCs(tss_wt, min_support = 3, min_tpm = 0.5)
# tss_tc_cold <- call_TCs(tss_cold, min_support = 3, min_tpm = 0.5)
# pas_tc_wt <- call_TCs(pas_wt, min_support = 3, min_tpm = 0.1)
# pas_tc_cold <- call_TCs(pas_cold, min_support = 3, min_tpm = 0.1)
# 
# export(tss_tc_wt, "TSS_TC_wt.bed", format = "BED")
# export(tss_tc_cold, "TSS_TC_cold.bed", format = "BED")
# export(pas_tc_wt, "PAS_TC_wt.bed", format = "BED")
# export(pas_tc_cold, "PAS_TC_cold.bed", format = "BED")
# 
# saveRDS(tss_tc_wt, "tss_tc_wt.RDS")
# saveRDS(tss_tc_cold, "tss_tc_cold.RDS")
# saveRDS(pas_tc_wt, "pas_tc_wt.RDS")
# saveRDS(pas_tc_cold, "pas_tc_cold.RDS")
# 
# tss_tc_wt <- readRDS("tss_tc_wt.RDS")
# tss_tc_cold <- readRDS("tss_tc_cold.RDS")
# pas_tc_wt <- readRDS("pas_tc_wt.RDS")
# pas_tc_cold <- readRDS("pas_tc_cold.RDS")
# 
# # Correct long reads:
# ont_wt <- extend_long_reads_to_TSS_and_PAS(ont_wt, tss_tc_wt, pas_tc_wt)
# ont_wt <- adjust_exons_of_long_reads(ont_wt)
# ont_wt <- detect_alignment_errors(ont_wt)
# 
# # Call de novo gene and transcript model:
# out_wt <- call_transcripts_and_genes(ont_wt)
# hml_genes_wt <- out_wt[[1]]
# hml_tx_wt <- out_wt[[2]]
# fusion_genes_wt <- out_wt[[3]]
# fusion_tx_wt <- out_wt[[4]]
# reads_free_wt <- out_wt[[5]]
# 
# # Refine by existing annotation (optional):
# #ref_wt <- refine_transcripts_by_annotation(hml_tx_wt, ebt, tss_tc_wt, pas_tc_wt, fusion_tx_wt)
# #hml_genes_wt <- ref_wt[[1]]
# #hml_tx_wt <- ref_wt[[2]]
# #fusion_genes_wt <- ref_wt[[3]]
# #fusion_tx_wt <- ref_wt[[4]]
# 
# # Also call transcription units from RNA-seq:
# trans_wt <- call_transcribed_intervals(rna_wt)
# transcribed_wt <- trans_wt[[1]]
# gaps_wt <- trans_wt[[2]]
# 
# export(transcribed_wt, "Transcribed_intervals_RNAseq_wt.bed", format = "BED")
# 
# results_wt <- process_nascent_intervals(hml_genes_wt, transcribed_wt, tss_tc_wt, pas_tc_wt, reads_free_wt, gaps_wt)
# hml_genes_v2_wt <- results_wt[[1]]
# lowexpr_wt <- results_wt[[3]]
# 
# rtracklayer::export(hml_genes_v2_wt, "Called_genes_wt.bed", format = "BED")
# write_grl_as_bed12(hml_tx_wt, "Called_transcripts_wt.bed")
# rtracklayer::export(fusion_genes_wt, "Fusion_genes_wt.bed", format = "BED")
# write_grl_as_bed12(fusion_tx_wt, "Fusion_transcripts_wt.bed")
# rtracklayer::export(lowexpr_wt, "Low_expressed_genes_wt.bed", format = "BED")
# 
# # Also save all results as RDS file:
# final_out_wt <- list("hml_genes_wt" = hml_genes_wt, "hml_genes_v2_wt" = hml_genes_v2_wt, "hml_tx_wt" = hml_tx_wt, "fusion_genes_wt" = fusion_genes_wt, "fusion_tx_wt" = fusion_tx_wt, "lncrna_wt" = lowexpr_wt)
# saveRDS(final_out_wt, "final_out_wt.RDS")
# 
# 
# ##### Also run TrRR on cold_3h data -----------------------------------------
# 
# # Correct long reads:
# ont_cold <- extend_long_reads_to_TSS_and_PAS(ont_cold, tss_tc_cold, pas_tc_cold)
# ont_cold <- adjust_exons_of_long_reads(ont_cold)
# ont_cold <- detect_alignment_errors(ont_cold)
# 
# # Call de novo gene and transcript model:
# out_cold <- call_transcripts_and_genes(ont_cold)
# hml_genes_cold <- out_cold[[1]]
# hml_tx_cold <- out_cold[[2]]
# fusion_genes_cold <- out_cold[[3]]
# fusion_tx_cold <- out_cold[[4]]
# reads_free_cold <- out_cold[[5]]
# 
# # Refine by existing annotation (optional):
# #ref_cold <- refine_transcripts_by_annotation(hml_tx_cold, ebt, tss_tc_cold, pas_tc_cold, fusion_tx_cold)
# #hml_genes_cold <- ref_cold[[1]]
# #hml_tx_cold <- ref_cold[[2]]
# #fusion_genes_cold <- ref_cold[[3]]
# #fusion_tx_cold <- ref_cold[[4]]
# 
# # Also call transcription units from RNA-seq:
# trans_cold <- call_transcribed_intervals(rna_cold)
# transcribed_cold <- trans_cold[[1]]
# gaps_cold <- trans_cold[[2]]
# 
# export(transcribed_cold, "Transcribed_intervals_RNAseq_cold.bed", format = "BED")
# 
# results_cold <- process_nascent_intervals(hml_genes_cold, transcribed_cold, tss_tc_cold, pas_tc_cold, reads_free_cold, gaps_cold)
# hml_genes_v2_cold <- results_cold[[1]]
# lowexpr_cold <- results_cold[[3]]
# 
# rtracklayer::export(hml_genes_v2_cold, "Called_genes_cold3h.bed", format = "BED")
# write_grl_as_bed12(hml_tx_cold, "Called_transcripts_cold3h.bed")
# rtracklayer::export(fusion_genes_cold, "Fusion_genes_cold3h.bed", format = "BED")
# write_grl_as_bed12(fusion_tx_cold, "Fusion_transcripts_cold3h.bed")
# rtracklayer::export(lowexpr_cold, "Low_expressed_genes_cold3h.bed", format = "BED")
# 
# # Also save all results as RDS file:
# final_out_cold <- list("hml_genes_cold" = hml_genes_cold, "hml_genes_v2_cold" = hml_genes_v2_cold, "hml_tx_cold" = hml_tx_cold, "fusion_genes_cold" = fusion_genes_cold, "fusion_tx_cold" = fusion_tx_cold, "lncrna_cold" = lowexpr_cold)
# saveRDS(final_out_cold, "final_out_cold.RDS")



########################### Modified Parameters #######################################################


##### Call TSS and PAS ----------------------------------------------

tss_tc_wt <- call_TCs(tss_wt, min_support = 3, min_tpm = 0.5, q_trim = 0.5) # default min_tpm = 0.5, q_trim = 0.05
tss_tc_cold <- call_TCs(tss_cold, min_support = 3, min_tpm = 0.5, q_trim = 0.5)
pas_tc_wt <- call_TCs(pas_wt, min_support = 3, min_tpm = 0.1, q_trim = 0.8) # default min_tpm = 0.1, q_trim = 0.05
pas_tc_cold <- call_TCs(pas_cold, min_support = 3, min_tpm = 0.1, q_trim = 0.8)

export(tss_tc_wt, "TSS_TC_wt_qTrim0.5.bed", format = "BED")
export(tss_tc_cold, "TSS_TC_cold_qTrim0.5.bed", format = "BED")
export(pas_tc_wt, "PAS_TC_wt_qTrim0.8.bed", format = "BED")
export(pas_tc_cold, "PAS_TC_cold_qTrim0.8.bed", format = "BED")

saveRDS(tss_tc_wt, "tss_tc_wt_qTrim0.5.RDS")
saveRDS(tss_tc_cold, "tss_tc_cold_qTrim0.5.RDS")
saveRDS(pas_tc_wt, "pas_tc_wt_qTrim0.8.RDS")
saveRDS(pas_tc_cold, "pas_tc_cold_qTrim0.8.RDS")

tss_tc_wt <- readRDS("tss_tc_wt_qTrim0.5.RDS")
tss_tc_cold <- readRDS("tss_tc_cold_qTrim0.5.RDS")
pas_tc_wt <- readRDS("pas_tc_wt_qTrim0.8.RDS")
pas_tc_cold <- readRDS("pas_tc_cold_qTrim0.8.RDS")


##### Run TrRR on wt data ---------------------------------------------

# Correct long reads:
ont_wt_data <- ont_wt
ont_wt <- extend_long_reads_to_TSS_and_PAS(ont_wt_data, tss_tc_wt, pas_tc_wt, read_flanks_down = c(-25, 25)) # default read_flanks_down = c(-50, 50)
ont_wt <- adjust_exons_of_long_reads(ont_wt)
ont_wt <- detect_alignment_errors(ont_wt, min_read_support = 5) # default min_read_support = 2

# Call de novo gene and transcript model:
out_wt_50 <- call_transcripts_and_genes(ont_wt, skip_minor_tx = 0.5) # default skip_minor_tx = 0.01

# write to file
hml_genes_wt <- out_wt_50[[1]]
hml_tx_wt <- out_wt_50[[2]]
fusion_genes_wt <- out_wt_50[[3]]
fusion_tx_wt <- out_wt_50[[4]]
reads_free_wt <- out_wt_50[[5]]

# Also call transcription units from RNA-seq:
trans_wt <- call_transcribed_intervals(rna_wt)
transcribed_wt <- trans_wt[[1]]
gaps_wt <- trans_wt[[2]]

export(transcribed_wt, "Transcribed_intervals_RNAseq_wt.bed", format = "BED")

results_wt <- process_nascent_intervals(hml_genes_wt, transcribed_wt, tss_tc_wt, pas_tc_wt, reads_free_wt, gaps_wt)
hml_genes_v2_wt <- results_wt[[1]]
lowexpr_wt <- results_wt[[3]]

rtracklayer::export(hml_genes_v2_wt, "Called_genes_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed", format = "BED")
write_grl_as_bed12(hml_tx_wt, "Called_transcripts_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed")
rtracklayer::export(fusion_genes_wt, "Fusion_genes_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed", format = "BED")
write_grl_as_bed12(fusion_tx_wt, "Fusion_transcripts_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed")
rtracklayer::export(lowexpr_wt, "Low_expressed_genes_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed", format = "BED")

# Also save all results as RDS file:
final_out_wt <- list("hml_genes_wt" = hml_genes_wt, "hml_genes_v2_wt" = hml_genes_v2_wt, "hml_tx_wt" = hml_tx_wt, "fusion_genes_wt" = fusion_genes_wt, "fusion_tx_wt" = fusion_tx_wt, "lncrna_wt" = lowexpr_wt)
saveRDS(final_out_wt, "final_out_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.RDS")


##### Run TrRR on cold data ---------------------------------------------

# Correct long reads:
ont_cold_data <- ont_cold
ont_cold <- extend_long_reads_to_TSS_and_PAS(ont_cold_data, tss_tc_cold, pas_tc_cold, read_flanks_down = c(-25, 25)) # default read_flanks_down = c(-50, 50)
ont_cold <- adjust_exons_of_long_reads(ont_cold)
ont_cold <- detect_alignment_errors(ont_cold, min_read_support = 5) # default min_read_support = 2

# Call de novo gene and transcript model:
out_cold_50 <- call_transcripts_and_genes(ont_cold, skip_minor_tx = 0.5) # default skip_minor_tx = 0.01

# write to file
hml_genes_cold <- out_cold_50[[1]]
hml_tx_cold <- out_cold_50[[2]]
fusion_genes_cold <- out_cold_50[[3]]
fusion_tx_cold <- out_cold_50[[4]]
reads_free_cold <- out_cold_50[[5]]

# Also call transcription units from RNA-seq:
trans_cold <- call_transcribed_intervals(rna_cold)
transcribed_cold <- trans_cold[[1]]
gaps_cold <- trans_cold[[2]]

export(transcribed_cold, "Transcribed_intervals_RNAseq_cold.bed", format = "BED")

results_cold <- process_nascent_intervals(hml_genes_cold, transcribed_cold, tss_tc_cold, pas_tc_cold, reads_free_cold, gaps_cold)
hml_genes_v2_cold <- results_cold[[1]]
lowexpr_cold <- results_cold[[3]]

rtracklayer::export(hml_genes_v2_cold, "Called_genes_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed", format = "BED")
write_grl_as_bed12(hml_tx_cold, "Called_transcripts_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed")
rtracklayer::export(fusion_genes_cold, "Fusion_genes_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed", format = "BED")
write_grl_as_bed12(fusion_tx_cold, "Fusion_transcripts_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed")
rtracklayer::export(lowexpr_cold, "Low_expressed_genes_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed", format = "BED")

# Also save all results as RDS file:
final_out_cold <- list("hml_genes_cold" = hml_genes_cold, "hml_genes_v2_cold" = hml_genes_v2_cold, "hml_tx_cold" = hml_tx_cold, "fusion_genes_cold" = fusion_genes_cold, "fusion_tx_cold" = fusion_tx_cold, "lncrna_cold" = lowexpr_cold)
saveRDS(final_out_cold, "final_out_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.RDS")