### helpers ####
# A function that decomposes an allele into the supplied motifs given either a
# character string or DNAString of the allele
#
# allele DNAString or character string of the Tandem repeat
#
# motifs a vector of strings or a DNAStringSet that are the expected
#
# motifs of the Tandem Repeat
#
# match_score for if motif matches sequence this is by default 1
#
# indel score for if there is an insertion or deletion in sequence by
# default this is -1
#
# allowence is the number of mismatches between the motif
# and the sequence permitted to be still be later encoded as a motif. By
# default this is 0
#
# This function returns the tandem repeat allele decomposed into its
# motifs as a character vector. For example, the allele `"ATATAT"` with motif
# `"AT"` is decomposed into `c("AT", "AT", "AT")` as well as an updated motifs
# DNAStringSet consisting of novel motifs detected in the allele.
#
# If no motifs are found in the allele from the initial motif set,
# the function returns `character(0)`and prints a message to inform the user.

# DP decomposition with indels
#assisted by Chat Gpt To translate TRviz algorithm into an R implementation then
#debugged heavily to allow for proper decomposition (given a skeleton)
decomposeTR <- function(allele, motifs, match_score = 1, indel = -1, allowence = 0) {

  #---------------------------------------------------------------
  # Input normalization and validation
  #---------------------------------------------------------------

  # Convert character inputs into Biostrings objects if needed
  if (is.character(allele)) allele <- DNAString(allele)
  if (is.character(motifs)) motifs <- DNAStringSet(motifs)

  # Ensure correct classes
  if(!is(allele, "DNAString"))
    stop("allele must be DNAString or character")

  if(!is(motifs, "DNAStringSet"))
    stop("motifs must be DNAStringSet or character vector")

  #length of allele
  allele_length <- length(allele)


  #---------------------------------------------------------------
  # Dynamic programming setup
  #---------------------------------------------------------------

  # dp_table[i] = best score achievable ending at allele position i-1
  # (dp_table is indexed from 1 to L+1)
  dp_table<- rep(-Inf, allele_length + 1)

  # traceback[i] stores which motif index or gap ('-') leads to DP[i]
  traceback <- rep(NA_character_, allele_length + 1)
  dp_table[1] <- 0 # Starting position has score 0

  #---------------------------------------------------------------
  # Precompute all motif matches along the allele
  #   - For each motif, find all approximate matches (allowence mismatches)
  #   - Store start/end coordinates and motif index
  #---------------------------------------------------------------

  hits_df <- do.call(rbind, lapply(seq_along(motifs), function(m) {
    hits <- matchPattern(motifs[[m]], allele, max.mismatch = allowence)

    # If this motif does not occur anywhere, skip
    if (length(hits) == 0){
      return(NULL)}

    pos_df <- data.frame(
      start = start(hits),
      end = end(hits),
      motif_idx = m
    )

    # Filter invalid coordinates
    pos_df <- pos_df[pos_df$start >= 1 & pos_df$end <= allele_length
                     & pos_df$start <= pos_df$end, , drop = FALSE]

    if (nrow(pos_df) == 0) return(NULL)

    pos_df
  }))

  # If no motif matches found, return empty output early
  if (is.null(hits_df) || nrow(hits_df) == 0) {
    message("No motif matches found in allele.")
    return(character(0))
  }

  # Sort hits by genomic position for sequential Dp_table processing
  hits_df <- hits_df[order(hits_df$start), ]

  #---------------------------------------------------------------
  # Dynamic programming forward pass
  #   For each position in the allele:
  #      1. Consider taking a gap/indel
  #      2. Consider matching a motif starting here
  #---------------------------------------------------------------

  for(pos in seq(along = allele) ){

    #-----------------------------------------------------------
    # Option 1: Move forward by 1 with an indel penalty
    #-----------------------------------------------------------

    if(dp_table[pos] != -Inf && pos + 1 <= allele_length+ 1) {
      score_gap <-dp_table[pos] + indel

      # Update dp_table if taking the gap is better
      if(score_gap > dp_table[pos + 1]) {
        dp_table[pos + 1] <- score_gap
        traceback[pos + 1] <- "-"
      }
    }

    #-----------------------------------------------------------
    # Option 2: Extend by matching a motif starting at 'pos'
    #-----------------------------------------------------------

    # Extract all motif hits that begin at this position
    hits_here <- hits_df[hits_df$start == pos, , drop = FALSE]

    for(i in seq_len(nrow(hits_here))) {
      hit <- hits_here[i, ]

      # Extract the motif's sequence and the allele slice
      motif_seq <- as.character(motifs[[hit$motif_idx]])
      subseq_allele <- as.character(subseq(allele, hit$start, hit$end))

      # Award a 'match' score if perfectly matching, otherwise 0
      score <- score_match_custom(subseq_allele, motif_seq, match_score)

      # New cumulative DP score if we take this motif
      new_score <-dp_table[pos] + score

      # dp_table state moves to the end of the motif + 1
      if(new_score > dp_table[hit$end + 1]) {

        dp_table[hit$end + 1] <- new_score
        traceback[hit$end + 1] <- as.character(hit$motif_idx)

      }
    }
  }

  #---------------------------------------------------------------
  # Reconstruct the allele sequence from the DP traceback
  #---------------------------------------------------------------
  composition <- reconstruct(allele_length, traceback, motifs, allele)

  #---------------------------------------------------------------
  # Detect and append novel motifs
  #---------------------------------------------------------------
  new_motifs <- setdiff(composition, motifs)
  motifs <- c(motifs, DNAStringSet(new_motifs))

  #---------------------------------------------------------------
  # Return final decomposition
  #---------------------------------------------------------------
  return(list(
    composition = composition,
    motifs = motifs
  ))
}

# Score a motif match against a sequence
#
# This helper function compares a candidate sequence against a motif.
# If the sequence perfectly matches the motif, it returns a specified score.
# Otherwise, it returns 0.
#
# seq1 Character string. The sequence segment from the allele to be compared.
# seq2 Character string. The motif sequence to match against seq1.
# match Numeric. The score to return if seq1 exactly matches seq2.
#
# Returns `match_score` if `seq1` and `seq2` are identical; otherwise returns 0.
score_match_custom <- function(seq1, seq2, match_score) {
  if (identical(seq1, seq2)) {
    return(match_score)
  } else {
    return(0)
  }
}

##Used ChatGPT to debug composition reconstruction from traceback and derived
#function below
reconstruct <- function(L, traceback, motifs, allele) {
  composition <- character()
  pos <- L + 1

  while(pos > 1) {
    t <- traceback[pos]
    if (is.na(t)) break

    if (t == "-") {
      # Find previous motif position
      prev_motif_pos <- which(!is.na(traceback[1:(pos - 1)]) & traceback[1:(pos - 1)] != "-")
      next_pos <- if(length(prev_motif_pos) == 0) 1 else max(prev_motif_pos)

      if(next_pos < pos) {
        gap_seq <- as.character(subseq(allele, next_pos, pos - 1))

        composition <- c(gap_seq, composition)
      }
      pos <- next_pos

    } else {
      # Motif match
      motif_idx <- suppressWarnings(as.integer(t))
      if (!is.na(motif_idx) && motif_idx <= length(motifs)) {
        k <- width(motifs[motif_idx])
        composition <- c(as.character(motifs[motif_idx]), composition)
        pos <- pos - k
      } else {
        # fallback
        composition <- c("-", composition)
        pos <- pos - 1
      }
    }
  }

  # Trailing insertion
  if(pos <= length(allele)) {
    trailing_seq <- as.character(subseq(allele, pos, length(allele)))
  }


  return(composition)
}

###User function####
#' Decomposition of Tandem Repeats
#'
#' A function that decomposes alleles into the supplied motifs given either a
#' character string vector or DNAStringSet of the alleles
#'
#' @param alleles DNAStringSet or character string vector of the Tandem repeats
#' @param motifs a vector of strings or a DNAStringSet that are the expected
#' motifs of the Tandem Repeat
#' @param match score for if motif matches sequence this is by default 1
#' @param indel score for if there is an insertion or deletion in sequence by
#' default this is -1
#' @param allowence is the number of mismatches between the motif
#' and the sequence permitted to be still be later encoded as a motif. By
#' default this is 0
#' @return This function returns the tandem repeat alleles as a vector
#' decomposed into there motifs as character vectors.
#' For example, the allele `"ATATAT"` with motif `"AT"` is decomposed
#' into `c("AT", "AT", "AT")`. And an updated DNAStringSet of mptifs containing
#' motifs across all alleles.
#'
#' If no motifs are found in an allele, the function returns `character(0)` for
#' that allele and prints a message to inform the user.
#'
#'
#' @references
#'
#' Pagès, H., Aboyoun, P., Gentleman, R. & DebRoy, S. Biostrings: Efficient
#' manipulation of biological strings (R package version 2.77.2, 2025).
#' https://bioconductor.org/packages/Biostrings, doi:10.18129/B9.bioc.Biostrings
#'
#' Park, J., Kaufman, E., Valdmanis, P. N. & Bafna, V. TRviz: A Python Library
#' for decomposing and Visualizing Tandem Repeat Sequences. Bioinformatics
#' Advances 3, (2023).
#'
#' R Core Team. R: A Language and Environment for Statistical Computing. Vienna,
#' Austria: R Foundation for Statistical Computing (2025).
#' https://www.R-project.org/
#'
#' OpenAI. ChatGPT (GPT-5) large language model (2025).
#'  https://chat.openai.com/
#'
#' @examples
#' \dontrun{
#' ##Example 1
#'
#' library(Biostrings)
#' alleles <- DNAStringSet(c("TTTACACGTAC", "ACACGTAC", "TTACACGTAC" ))
#' motifs <- DNAStringSet(c("AC", "GT"))
#' compositions <- decomposeTRs(alleles, motifs)
#' compositions
#'
#' ###Expected:
#' $compositions
#' $compositions[[1]]
#' [1] "TTT" "AC"  "AC"  "GT"  "AC"
#'
#' $compositions[[2]]
#' [1] "AC" "AC" "GT" "AC"
#'
#' $compositions[[3]]
#' [1] "TT" "AC" "AC" "GT" "AC"
#'
#'
#' $motifs
#' DNAStringSet object of length 4:
#'   width seq
#' [1]     2 AC
#' [2]     2 GT
#' [3]     3 TTT
#' [4]     2 TT
#' }
#'
#' @import Biostrings
#' @importFrom methods is
#' @export decomposeTRs

decomposeTRs <- function(alleles, motifs, match = 1,
                         indel = -1, allowence = 0) {

  if (is.character(alleles)) alleles <- DNAStringSet(alleles) #from Biostrings
  if (is.character(motifs)) motifs <- DNAStringSet(motifs)  #from Biostrings

  if(!is(alleles, "DNAStringSet")){
    stop("Alleles should be character vector or DNAStringset")
  }

  if(!is(motifs, "DNAStringSet")){
    stop("Motifs should be character vector or DNAStringset")
  }

  TRs <- lapply(alleles, function(allele) {
    decomposeTR(allele, motifs, match, indel, allowence)
  })

  all_compositions <- lapply(TRs, `[[`, "composition")
  all_motifs <- lapply(TRs, `[[`, "motifs")
  combined_motifs <- unique(do.call(c, all_motifs))

  return(list(
    compositions = all_compositions,
    motifs = combined_motifs
  ))
}





