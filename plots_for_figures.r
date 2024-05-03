library(tidyverse)
library(readxl)
library(R.utils)
library(GenomicFeatures)
library(rtracklayer)
library(devtools)

devtools::source_url("https://github.com/Maxim-Ivanov/Ivanov_et_al_2021/blob/main/07-Custom_functions.R?raw=TRUE")
devtools::source_url("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/merge_and_normalize_GRanges.R?raw=TRUE")

# Load the new annotation returned by TranscriptomeReconstructoR:
final_out <- readRDS("path_to_file/final_out.RDS")
hml_genes <- final_out$hml_genes
hml_genes_v2 <- final_out$hml_genes_v2
hml_tx <- final_out$hml_tx
lncrna <- final_out$lncrna

# Load EnsemblPlants genes:
txdb_Ensmbl <- makeTxDbFromGFF("path_to_gff_file/input.gff3") ### annotation file from the database to compare against TrRR annotation
genes_Ensmbl <- genes(txdb_Ensmbl, columns = c("gene_id", "tx_type"))
mcols(genes_Ensmbl)$tx_type <- mcols(genes_Ensmbl)$tx_type %>% unlist() %>% unname()
ebt_Ensmbl <- exonsBy(txdb_Ensmbl, by = "gene", use.names = FALSE)

# Load TSS-seq - STRIPE, 3'DRS-seq - Quant-Seq, plaNET-seq and pNET-seq data:
tss_data <- read_stranded_bedGraph("STRIPE-Seq_merged.bedgraph.gz", seqinfo = seqinfo(txdb_Ensmbl))
drs_data <- read_stranded_bedGraph("Quant-Seq_merged.bedgraph.gz", seqinfo = seqinfo(txdb_Ensmbl))
rnaseq_data <- read_stranded_bedGraph("ssRNA-Seq_merged.bedgraph.gz", seqinfo = seqinfo(txdb_Ensmbl))


### Analyze overlaps between called genes and known genes -----------------------------------

tbl1 <- classify_genes_by_overlap_with_ref_genes(hml_genes_v2, genes_Ensmbl) %>% as_tibble() %>% dplyr::select(type, Overlap) %>% mutate(ann = "Ensmbl")

# Fig. 2A:

tbl1$ann <- tbl1$ann %>% as_factor() %>% relevel(ref = "Ensmbl")
tbl1$type <- tbl1$type %>% factor(levels = c("HC", "MC", "LC"))
title <- "Overlaps between called and known genes"
p <- ggplot(tbl1, aes(x = type, fill = Overlap)) + geom_bar(colour = "black") + xlab(NULL) + ylab("Number of genes") + ggtitle(title) +
  scale_fill_manual(values = c(No = "white", Unique = "grey40", Multiple = "grey80")) + theme_bw() + facet_grid(. ~ ann)
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(title, ext), plot = p, width = 5, height = 5, units = "in")
}



### Find matched pairs of called/known genes ---------------------------------------------------

hml_par1 <- classify_genes_by_overlap_with_ref_genes(hml_genes_v2, genes_Ensmbl) %>% `[`(mcols(.)$Overlap == "Unique")
Ensmbl_par <- findOverlaps(hml_par1, genes_Ensmbl, select = "first") %>% genes_Ensmbl[.]

### Find matched duos of called/Ensemble genes ------------------------------------------

t1 <- tibble(gene_id = mcols(hml_par1)$name, idx1 = 1:length(hml_par1))
hml_duo <- hml_par1[t1$idx1]
mcols(hml_duo)$Ensmbl_mate <- Ensmbl_par[t1$idx1] %>% granges()

### Compare borders of called genes vs Ensembl ------------------------------------------

grp1 <- mcols(hml_par1)$type %>% factor(levels = c("HC", "MC", "LC"))
# Fig. 2C (upper):
compare_gene_borders(hml_par1, Ensmbl_par, groups = grp1, title = "Difference of gene borders (called genes vs Ensmbl)", width = 6, height = 4)
# Fig. 2B (left):
plot_percent_overlap(hml_par1, Ensmbl_par, groups = grp1, title = "Percent overlap (called genes vs Ensmbl)", width = 4, height = 6)
compute_stats_on_percent_overlap(hml_par1, Ensmbl_par, threshold = 0.9)

### Metagene plot of independent TSS and PAS signal tracks ---------------------------------------

# Subset duos to HC genes:
hc_duo <- hml_duo[mcols(hml_duo)$type == "HC"]

# Generate windows around TSS and PAS:
tss_win_hc <- granges(hc_duo) %>% generate_windows_around_tss_or_pas(mode = "start", width = 50)
tss_win_Ensmbl <- mcols(hc_duo)$Ensmbl_mate %>% generate_windows_around_tss_or_pas(mode = "start", width = 50)

pas_win_hc <- granges(hc_duo) %>% generate_windows_around_tss_or_pas(mode = "end", width = 50)
pas_win_Ensmbl <- mcols(hc_duo)$Ensmbl_mate %>% generate_windows_around_tss_or_pas(mode = "end", width = 50)

###########################################################################################################

# Draw metagenes around TSS:
tss_m1 <- metageneMatrix(tss_data, tss_win_hc, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE) # seqlevels(tss_data) <- seqlevels(tss_win_hc) ;  seqinfo(tss_data) <- seqinfo(tss_win_hc)
tss_m2 <- metageneMatrix(tss_data, tss_win_Ensmbl, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE) #  seqlevels(tss_win_Ensmbl) <- seqlevels(tss_win_hc)  ; seqinfo(tss_win_Ensmbl) <- seqinfo(tss_win_hc)
tss_matlist <- list("HC genes" = tss_m1, "Ensmbl" = tss_m2)
drawMetagenePlot(tss_matlist, x.axis = seq(-25, 24), vline = 0, linetype = "dotted", title = "STRIPE-seq around gene starts",
                 xlabel = "Windows 50 bp centered at gene starts", ylabel = "STRIPE-seq signal", width = 5, height = 8, units = "in") # Fig. 2D

###-------------- Work on ----------------------------------------------------------------------------------

# Draw metagenes around PAS:
pas_m1 <- metageneMatrix(drs_data, pas_win_hc, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE) # seqlevels(drs_data) <- seqlevels(pas_win_hc)  ;  seqinfo(drs_data) <- seqinfo(pas_win_hc)
pas_m2 <- metageneMatrix(drs_data, pas_win_Ensmbl, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE)  # seqlevels(pas_win_Ensmbl) <- seqlevels(pas_win_hc) ;  seqinfo(pas_win_Ensmbl) <- seqinfo(pas_win_hc)
pas_matlist <- list("HC genes" = pas_m1, "Ensmbl" = pas_m2)
drawMetagenePlot(pas_matlist, x.axis = seq(-25, 24), vline = 0, linetype = "dotted", title = "Quant-seq around gene ends",
                 xlabel = "Windows 50 bp centered at gene ends", ylabel = "Quant-seq signal", width = 5, height = 8, units = "in") # Fig. 2E

#############################################################################################################


##### Analyze internal exons in re-discovered known genes (only matched pairs of genes are considered) ---------------------------------

# Fig. 3B:
out1_m <- find_novel_internal_exons(hml_tx, ebt_Ensmbl, hml_par1, Ensmbl_par, xlab_nbreaks = 8, title = "Difference of matched exon borders in matched gene pairs (called vs Ensmbl)") 

# Fig. 3A:
tbl <- tibble("Ensmbl" = lapply(out1_m, length), type = names(out1_m)) %>% 
  pivot_longer(cols = c("Ensmbl"), names_to = "ann") %>% arrange(ann)
tbl$ann <- tbl$ann %>% as_factor() %>% relevel(ref = "Ensmbl")
#tbl$type <- tbl$type %>% factor(levels = c("No overlap", "Exact match", "IR", "Alt donor", "Alt acceptor", "Other"))
title <- "Classification of internal exons in matched gene pairs"
p <- ggplot(tbl, aes(x = ann, y = value, fill = type)) + geom_bar(stat = "identity", width = 0.6, colour = "white") + xlab(NULL) + 
  ylab("Number of exons") + theme_bw() + theme(legend.title = element_blank()) + ggtitle(title)
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(title, ext), plot = p, width = 4, height = 6, units = "in")
}


#################################################################################################################################

### Annotate called genes by overlap with EnsmblPlants -------------------------------------------

# Find novel genes which do not overlaps with any known gene:
mcols(hml_genes_v2)$novel <- hml_genes_v2 %outside% genes_Ensmbl
hml_novel <- hml_genes_v2[mcols(hml_genes_v2)$novel]
length(hml_novel)
table(mcols(hml_novel)$type)
mcols(lncrna)$novel <- lncrna %outside% genes_Ensmbl
lncrna_novel <- lncrna[mcols(lncrna)$novel]
length(lncrna_novel)

hml_novel_NOTas <- subset(hml_novel, as %in% c("FALSE"))
rtracklayer::export(hml_novel_NOTas, "Novel_genes.bed", format = "BED")

# How many of them are antisense to known genes?
mcols(hml_novel)$as <- overlapsAny(hml_novel, genes_Ensmbl, ignore.strand = TRUE)
sum(mcols(hml_novel)$as)
table(mcols(hml_novel)$type[mcols(hml_novel)$as])
mcols(lncrna_novel)$as <- overlapsAny(lncrna_novel, genes_Ensmbl, ignore.strand = TRUE)
sum(mcols(lncrna_novel)$as)

# Antisense to known genes only:
hml_as2 <- overlapsAny(hml_novel, genes_Ensmbl, ignore.strand = TRUE)
table(mcols(hml_novel)$type[hml_as2])
lnc_as2 <- overlapsAny(lncrna_novel, genes_Ensmbl, ignore.strand = TRUE)
sum(lnc_as2)

# Barplot of percent novel intergenic and antisense genes (Fig. 4A):
hc <- mcols(hml_novel)$type == "HC"
mc <- mcols(hml_novel)$type == "MC"
lc <- mcols(hml_novel)$type == "LC"
hml_as_Ensmbl <- overlapsAny(hml_novel, genes_Ensmbl, ignore.strand = TRUE)
lnc_as_Ensmbl <- overlapsAny(lncrna_novel, genes_Ensmbl, ignore.strand = TRUE)


tbl <- tibble(total = c(sum(hc), sum(mc), sum(lc), length(lncrna_novel)), 
              as_Ensmbl = c(sum(hc & hml_as_Ensmbl), sum(mc & hml_as_Ensmbl), sum(lc & hml_as_Ensmbl), sum(lnc_as_Ensmbl)),
              type = c("HC", "MC", "LC", "tr.RNA"))
tbl2 <- tbl %>% mutate(ig_Ensmbl = total - as_Ensmbl) %>% dplyr::select(-total) %>% 
  pivot_longer(cols = c(as_Ensmbl, ig_Ensmbl)) %>% separate(name, into = c("as", "ann"), sep = "_")
tbl2$as <- ifelse(tbl2$as == "as", "Antisense", "Intergenic") %>% as_factor() %>% relevel(ref = "Antisense")
tbl2$type <- tbl2$type %>% factor(levels = c("HC", "MC", "LC", "tr.RNA"))
tbl2$ann <- ifelse(tbl2$ann == "Ensmbl", "EnsemblPlants") %>% as_factor() %>% relevel(ref = "EnsemblPlants")
title <- "Novel genes and lncRNAs vs EnsemblPlants"
p <- ggplot(tbl2, aes(x = type, y = value, fill = as)) + geom_bar(stat = "identity", colour = "white") + facet_wrap(. ~ ann) + 
  theme_bw() + theme(legend.title = element_blank()) + xlab(NULL) + ylab("Number of genes") + ggtitle(title)
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(title, ext), plot = p, width = 5, height = 4, units = "in")
}

### Confirm novel genes by ssRNA-seq (Fig. 4B) ------------------------------------------------------------

m1 <- metagene_matrix_with_flanks(rnaseq_data, hml_novel)
m2 <- metagene_matrix_with_flanks(rnaseq_data, lncrna_novel)
ml1 <- list("Novel HC MC LC genes" = m1, "Novel SMC" = m2)
drawMetagenePlot(ml1, title = "Metagene ssRNA-seq novel genes", x.axis = seq(-19, 120), vline = c(0, 100), xlabel = "Scaled genes with 100 bp flanks", width = 7, height = 4, units = "in")



##### Metagene plot of plaNET-seq signal at donor sites of the novel exons ----------------------------------------------

out_Ensembl <- classify_exons(hml_tx, ebt_Ensmbl)
lapply(out_Ensembl, length)

# Fig. 3C (left):
ml1 <- vector("list", length(out_Ensembl))
names(ml1) <- names(out_Ensembl)
for (i in seq_along(out_Ensembl)) {
  exons <- out_Ensembl[[i]]
  windows <- resize(exons, 1, "end") %>% resize(50, "center")
  ml1[[i]] <- metageneMatrix(rnaseq_data, windows, scaling = FALSE, skip.zeros = FALSE)
}
drawMetagenePlot(ml1, title = "ssRNA-seq signal over donor sites (EnsemblPlants)", xlabel = "50 bp windows centered on donor sites", vline = 0, x.axis = seq(-24, 25), width = 5, height = 6, units = "in")
