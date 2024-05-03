library(GenomicFeatures)
library(rtracklayer)
library(devtools)
library(tidyverse)

convert_GAlignments_to_coverage <- function(ga, split = TRUE, mode = "whole_read", merge.strands = FALSE, 
                                            flip.strands = FALSE, normalize = FALSE, norm_to = 1000000L, precision = 6) {
  require(GenomicAlignments)
  if (mode %in% c("start", "end")) {
    gr <- resize(granges(ga), width = 1, fix = mode)
  } else if (mode == "whole_read") {
    if (isTRUE(split)) {
      gr <- unlist(grglist(ga), use.names = FALSE)
    } else {
      gr <- granges(ga)
    }
  }
  if (isTRUE(merge.strands)) {
    cov_all <- coverage(gr)
    out <- bindAsGRanges(score = cov_all)
  } else {
    cov_fw <- coverage(gr[strand(gr) == "+"])
    cov_rev <- coverage(gr[strand(gr) == "-"])
    out_fw <- bindAsGRanges(score = cov_fw)
    strand(out_fw) <- "+"
    out_rev <- bindAsGRangess(score = cov_rev)
    strand(out_rev) <- "-"
    out <- c(out_fw, out_rev)
    if (isTRUE(flip.strands)) {
      strand(out) <- ifelse(strand(out) == "+", "-", "+")
    }
  }
  out <- out[score(out) > 0]
  out <- sortSeqlevels(out)
  out <- out[order(seqnames(out), start(out))]
  if (isTRUE(normalize)) {
    norm_factor <- length(ga) / norm_to
    score(out) <- round(score(out) / norm_factor, precision)
  }
  return(out)
}

RNAseq_BAM_to_normalized_bedGraph <- function(bamfile, mode = "PE", stranded = TRUE, switch_strand = TRUE, normalize = TRUE, norm_to = 1e06, skip_over = NULL) {
  stopifnot(mode %in% c("PE", "SE"))
  if (mode == "PE") {
    ga <- readGAlignmentPairs(bamfile)
  } else {
    ga <- readGAlignments(bamfile)
  }
  if (!is.null(skip_over) && class(skip_over) == "GRanges") {
    if (!isTRUE(stranded)) {
      strand(skip_over) <- "*"
    } else if (isTRUE(switch_strand)) {
      strand(skip_over) <- ifelse(strand(skip_over) == "+", "-", "+")
    }
    bad <- ga %over% skip_over
    message(sum(bad), " (", round(mean(bad) * 100, 1), "%) reads were skipped due to overlap with unwanted intervals;")
    ga <- ga[!bad]
  }
  cov <- convert_GAlignments_to_coverage(ga, merge.strands = !stranded, flip.strands = switch_strand, normalize = normalize, norm_to = norm_to)
  return(cov)
}


filenames <- list.files(".", pattern="dedup.bam$")

for (f in filenames) {
  print(f)
  cov <- RNAseq_BAM_to_normalized_bedGraph(f)
  rtracklayer::export(cov, paste0(f, ".bedgraph"), format = "bedGraph")
}
    

