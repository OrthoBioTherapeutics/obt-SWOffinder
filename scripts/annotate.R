library(foreach)
library(doParallel)
library(GenomicRanges)

load("../../misc/gtf/canFam3_flattened_gtf.RData")
data <- read.csv("/mnt/new_home/obt-guide-design/SWOffinder/out/4mm4w2b_AGAACCGCAGAGCTAAATGCNGG_processed.csv")
out_name <- "OCC02S_annotated_offtargets.csv"


# Create a GRanges object from Chromosome and EndPosition columns
gr <- GRanges(seqnames = Rle(data$Chromosome),
              ranges = IRanges(start = data$EndPosition, end = data$EndPosition),
              mcols = data[ , !(colnames(data) %in% c("Chromosome", "EndPosition"))])

# Initialize the 'region' and 'distance_to_gene' columns
mcols(gr)$region <- NA
mcols(gr)$distance_to_gene <- NA
mcols(gr)$gene <- NA  # Initialize the 'gene' column with NA
mcols(gr)$nearest_gene <- NA

gene_ranges_gr <- gene_granges
gene_ranges_gr$gene_id <- names(gene_ranges_gr)

overlaps <- findOverlaps(gr, gene_ranges_gr)
overlapping_genes <- gene_ranges_gr$gene_id[subjectHits(overlaps)]
mcols(gr)$gene[queryHits(overlaps)] <- overlapping_genes[seq_along(queryHits(overlaps))]

# Update rows with overlaps
mcols(gr)$region[queryHits(overlaps)] <- "genic"
mcols(gr)$gene[queryHits(overlaps)] <- gene_ranges_gr$gene_id[subjectHits(overlaps)]

# For rows with no overlaps, mark as intergenic and calculate distance to closest gene
no_overlap <- setdiff(seq_along(gr), queryHits(overlaps))
nearest_genes <- nearest(gr[no_overlap], gene_ranges_gr)
mcols(gr)$region[no_overlap] <- "intergenic"
valid_indices <- !is.na(nearest_genes)
mcols(gr)$distance_to_gene[no_overlap[valid_indices]] <- 
distance(gr[no_overlap[valid_indices]], gene_ranges_gr[nearest_genes[valid_indices]])
mcols(gr)$nearest_gene[no_overlap[valid_indices]] <- 
mcols(gene_ranges_gr)$gene_id[nearest_genes[valid_indices]]

for (i in queryHits(overlaps)) {
  if (mcols(gr)$mcols.score[i] > 0.01) { 
    gene_name <- mcols(gr)$gene[i]
    print(gene_name)
    if (!is.na(gene_name)) {
      gene_exons <- merged_exons_list[[gene_name]]
    if (length(findOverlaps(gr[i], gene_exons))!=0) {
      mcols(gr)$region[i] <- "exon"
    } else {
      mcols(gr)$region[i] <- "intron"
    }
  }
}
}


write.csv(as.data.frame(gr), out_name, row.names = FALSE)
