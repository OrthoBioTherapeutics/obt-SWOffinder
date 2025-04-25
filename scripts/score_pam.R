library(dplyr)
library(stringr)
library(parallel)

complement <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G", "-" = "-")

# Get command-line arguments (excluding the first default ones)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No CSV files provided.")
}
csv_files <- args

calculate_mismatch_score <- function(spacer, protospacer, scores) {
  # Calculate the PAM score (for the last two nucleotides of the protospacer)
  pam_dinucleotide <- substr(protospacer, nchar(protospacer) - 1, nchar(protospacer))
  pam_score_value <- scores %>%
    filter(pos == 22, DNA == pam_dinucleotide) %>%
    dplyr::select(score) %>%
    as.numeric()
  
  # Remove the PAM part from both sequences
  spacer <- substr(spacer, 0, nchar(spacer)-3)
  protospacer <- substr(protospacer, 0, nchar(protospacer)-3)
  
  # Split sequences into characters
  spacer_split <- strsplit(spacer, "")[[1]]
  protospacer_split <- strsplit(protospacer, "")[[1]]

  # Initialize variables
  final_score <- 1
  position <- 0  # Keeps track of the last non-gap position
  
  for (i in seq_along(spacer_split)) {
    if (spacer_split[i] != "-") {
      position = position +1
    }

    # Check if there is a mismatch at this position
    if (spacer_split[i] != protospacer_split[i]) {
      dna_base <- complement[protospacer_split[i]]
      rna_base <- spacer_split[i]

      # Fetch the score for the current position
      score_value <- scores %>%
        filter(RNA == rna_base, DNA == dna_base, pos == position) %>%
        dplyr::select(score) %>%
        as.numeric()
      print(score_value)
      # Multiply the score
      final_score <- final_score * score_value
    }
  }
  
  # Multiply final score by PAM score value
  final_score <- final_score * pam_score_value
  return(final_score)
}

process_file <- function(file, scores) {
  swofftargets <- read.csv(paste0(file))
  swofftargets <- swofftargets %>%
    rowwise() %>%
    mutate(score = calculate_mismatch_score(AlignedTarget, AlignedText, scores))
  
  # Calculate the final result
  result <- 100 / (100 + sum(swofftargets$score))
  
  print(paste("Result for", file, ":", result))
  write.csv(swofftargets, paste0(file, "_processed.csv"), row.names = FALSE)
}

scores <- read.csv("../misc/gap_aware_matrix.csv")
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)

# Export necessary objects and functions to the cluster
clusterExport(cl, list("complement", "calculate_mismatch_score", "process_file", "scores"))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(stringr))

parLapply(cl, csv_files, process_file, scores)

# Stop the cluster
stopCluster(cl)
