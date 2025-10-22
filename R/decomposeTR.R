#' Tandem Repeat Decomposition
#'
#' A function that decomposes an allele into the supplied motifs given either a
#' character string or DNAString of the allele
#'
#' @param allele DNAString or character string of the Tandem repeat
#' @param motifs a vector of strings or a DNAStringSet that are the expected
#' motifs of the Tandem Repeat
#' @param match score for if motif matches sequence this is by default 1
#' @param indel score for if there is an insertion or deletion in sequence by
#' default this is -1
#' @param allowence is the number of mismatches between the motif
#' and the sequence permitted to be still be later encoded as a motif. By
#' default this is 0
#' @return This function returns the tandem repeat allele decomposed into its motifs
#' as a character vector. For example, the allele `"ATATAT"` with motif `"AT"`
#' is decomposed into `c("AT", "AT", "AT")`.
#'
#' If no motifs are found in the allele, the function returns `character(0)`
#' and prints a message to inform the user.
#'
#' Any indels within the repeat are represented as `"-"` in the output vector.
#' Each `"-"` represents an indel the length of the preceding motif unit.
#' For example, `"ATATCAT"` with motif `"AT"` is decomposed into
#' `c("AT", "AT", "-", "AT")`.
#' @examples
#' \dontrun{
#' ##Example 1
#'
#' library(Biostrings)
#' allele <- DNAString("TTTACACGTAC")
#' motifs <- DNAStringSet(c("AC", "GT"))
#' composition <- decomposeTR(allele, motifs)
#' composition
#'
#' ###Expected: "-"  "-"  "AC" "AC" "GT" "AC"
#' }
#'
#' \dontrun{
#' ##Example 2
#'
#' allele <- DNAString("TTTACACGTCAC")
#' motifs <- DNAStringSet(c("AC", "GTC"))
#' composition <- decomposeTR(allele, motifs)
#' composition
#'
#' ###Expected "-"   "-"   "AC"  "AC"  "GTC" "AC"
#'}
#' @references OpenAI. ChatGPT (GPT-5) large language model (2025).
#' https://chat.openai.com/
#'
#'Pagès, H., Aboyoun, P., Gentleman, R. & DebRoy, S. Biostrings: Efficient
#'manipulation of biological strings (R package version 2.77.2, 2025).
#'https://bioconductor.org/packages/Biostrings, doi:10.18129/B9.bioc.Biostrings
#'
#' Park, J., Kaufman, E., Valdmanis, P. N. & Bafna, V. TRviz: A Python Library
#' for decomposing and Visualizing Tandem Repeat Sequences. Bioinformatics
#' Advances 3, (2023).
#'
#' @import Biostrings
#' @export decomposeTR

# DP decomposition with indels
#assisted by Chat Gpt To translate TRviz algorithm into an R implementation then
#debugged heavily to allow for proper decomposition (given a skeleton)
decomposeTR <- function(allele, motifs, match = 1,
                        indel = -1, allowence = 0) {
  if (is.character(allele)) allele <- DNAString(allele) #from Biostrings
  if (is.character(motifs)) motifs <- DNAStringSet(motifs) #from Biostrings

  L <- length(allele)
  M <- length(motifs)

  DP <- rep(-Inf, L + 1)
  traceback <- rep(NA_character_, L + 1)
  DP[1] <- 0

  # Precompute motif hits
  hits_df <- do.call(rbind, lapply(seq_along(motifs), function(m) {
    hits <- matchPattern(motifs[[m]], allele, max.mismatch = allowence)
    #matchPattern is from the BioStrings package
    if(length(hits) == 0) return(NULL)
    data.frame(
      start = start(hits),
      end = end(hits),
      motif_idx = m
    )
  }))

  if (is.null(hits_df) || nrow(hits_df) == 0) {
    message("No motif matches found in allele.")
    return(character(0))
  }

  hits_df <- hits_df[order(hits_df$start), ]

  # Fill DP
  for(pos in 1:L) {
    # Indel: gap of 1 base
    if(DP[pos] != -Inf && pos + 1 <= L + 1) {
      score_gap <- DP[pos] + indel
      if(score_gap > DP[pos + 1]) {
        DP[pos + 1] <- score_gap
        traceback[pos + 1] <- "-"
      }
    }

    # Motif hits starting at this position
    hits_here <- hits_df[hits_df$start == pos, , drop = FALSE]
    for(i in seq_len(nrow(hits_here))) {

      s <- hits_here[i, ]
      motif_seq <- as.character(motifs[[s$motif_idx]])
      subseq_allele <- as.character(subseq(allele, s$start, s$end))
      score <- sum(score_match_custom(subseq_allele, motif_seq, match))
      new_score <- DP[pos] + score

      if(new_score > DP[s$end + 1]) {
        DP[s$end + 1] <- new_score
        traceback[s$end + 1] <- as.character(s$motif_idx)
      }
    }
  }
  composition <- reconstruct(L, traceback, motifs)
  return(composition)
}

### helpers ####
score_match_custom <- function(seq1, seq2, match) {
  # Returns 'match' if the entire motif matches, else 0
  if (identical(seq1, seq2)) {
    return(match)
  } else {
    return(0)
  }
}

##Used ChatGPT to debug composition reconstruction from traceback and derived
#function below
reconstruct <- function(L, traceback, motifs) {
  composition <- character()
  pos <- L + 1

  while (pos > 1) {
    t <- traceback[pos]
    if (is.na(t)) break

    if (t == "-") {
      # Find next motif position before the current gap
      next_motif_pos <- which(!is.na(traceback[1:(pos - 1)]) & traceback[1:(pos - 1)] != "-")

      if (length(next_motif_pos) == 0) {
        # Leading gap: use first motif length
        k <- width(motifs[1])
        next_pos <- 1
      } else {
        next_pos <- max(next_motif_pos)
        t_next <- traceback[next_pos]
        if (is.na(t_next) || t_next == "-") {
          k <- 1
        } else {
          # Safely convert and check index
          motif_idx <- suppressWarnings(as.integer(t_next))
          if (is.na(motif_idx) || motif_idx > length(motifs)) {
            k <- 1
          } else {
            k <- width(motifs[motif_idx])
          }
        }
      }

      # Compute gap length in bases
      gap_len <- pos - next_pos
      n_units <- ceiling(gap_len / k)
      composition <- c(rep("-", n_units), composition)
      pos <- pos - gap_len

    } else {
      # Motif match
      motif_idx <- suppressWarnings(as.integer(t))
      if (!is.na(motif_idx) && motif_idx <= length(motifs)) {
        k <- width(motifs[motif_idx])
        composition <- c(as.character(motifs[motif_idx]), composition)
        pos <- pos - k
      } else {
        # Fallback for corrupted traceback entries
        composition <- c("-", composition)
        pos <- pos - 1
      }
    }
  }

  return(composition)
}







