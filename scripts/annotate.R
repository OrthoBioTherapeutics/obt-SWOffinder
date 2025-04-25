args <- commandArgs(trailingOnly = TRUE)

# Expected arguments
gtf_path <- args[1]
in_df <- args[2]
guidename <- args[3]

if (length(args) < 3) {
  stop("Usage: Rscript annotate.R <gtf_path> <in_df> <guidename>")
}

library(foreach)
library(doParallel)
library(GenomicRanges)
library(dplyr)
library(GenomicFeatures)


data <- read.csv(in_df)
out_name <- paste0(guidename, "gap_ot.csv")

# Account for offset which occurs for some reason. This is to make it comparable with guide pipeline output.
data <- data %>%
  mutate(EndPosition = case_when(
    Strand == "+" ~ EndPosition - 5,
    Strand == "-" ~ EndPosition - 17
  ))

txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")

genes_gr <- unlist(genes(txdb, single.strand.genes.only=FALSE))
exons_gr <- exons(txdb)
introns_gr <- unlist(intronsByTranscript(txdb), use.names = FALSE)

# Create a GRanges object from Chromosome and EndPosition columns
gr_ot <- GRanges(seqnames = Rle(data$Chromosome),
              ranges = IRanges(start = data$EndPosition, end = data$EndPosition),
              mcols = data[ , !(colnames(data) %in% c("Chromosome", "EndPosition"))])

# Default annotations
mcols(gr_ot)$region <- "intergenic"
mcols(gr_ot)$distance_to_gene <- NA
mcols(gr_ot)$gene <- NA
mcols(gr_ot)$nearest_gene <- NA

# Initialize columns
gr_ot$region <- "intergenic"
gr_ot$gene <- NA_character_

# Mark exonic
hits_exon <- findOverlaps(gr_ot, exons_gr)
gr_ot$region[queryHits(hits_exon)] <- "exon"

# Mark intronic (only for non-exonic regions)
non_exonic_idx <- which(gr_ot$region != "exon")
hits_intron <- findOverlaps(gr_ot[non_exonic_idx], introns_gr)
gr_ot$region[non_exonic_idx[queryHits(hits_intron)]] <- "intron"

# Get gene_id by overlapping with gene annotations
hits_gene <- findOverlaps(gr_ot, genes_gr)
gr_ot$gene[queryHits(hits_gene)] <- names(genes_gr)[subjectHits(hits_gene)]

intergenic_idx <- which(gr_ot$region == "intergenic")
gr_intergenic <- gr_ot[intergenic_idx]
nearest_hits <- nearest(gr_intergenic, genes_gr)
closest_gene_ids <- names(genes_gr)[nearest_hits]

# Get intergenic ranges
intergenic_idx <- which(gr_ot$region == "intergenic")
gr_intergenic <- gr_ot[intergenic_idx]

# Find nearest gene and compute distances
nearest_hits <- nearest(gr_intergenic, genes_gr)
closest_gene_ids <- names(genes_gr)[nearest_hits]
distances <- rep(NA_integer_, length(gr_intergenic))

# Compute distances only where nearest gene exists
valid <- !is.na(nearest_hits)
distances[valid] <- distance(gr_intergenic[valid], genes_gr[nearest_hits[valid]])

# Add to original GRanges
gr_ot$nearest_gene[intergenic_idx] <- closest_gene_ids
gr_ot$distance_to_gene[intergenic_idx] <- distances

write.csv(as.data.frame(gr_ot), paste0("../out/", out_name), row.names = FALSE)

