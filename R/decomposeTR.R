#' Decomposes allele into motifs
#'
#' A function that decomposes an allele into the supplied motifs given either a
#' character string or DNAString of the allele
#'
#' @param allele DNAString or character string of the Tandem repeat
#' @param motifs a vector of strings or a DNAStringSet that are the expected
#' motifs of the Tandem Repeat
#' @param match score for if motif matches sequence this is by default 1
#' @param mis_match score for if there is a mismatch between sequence and motifs
#' by default this is -1
#' @param indel score for if there is an insertion or deletion in sequence by
#' default this is -2
#' @param allowence is the number of mismatches between the motif
#' and the sequence permitted to be still be later encoded as a motif by default
#' this is 0
#' @references OpenAI. (2025). ChatGPT (GPT-5) [Large language model]. https://chat.openai.com/
#' @import Biostrings
#' @export decomposeTR
# DP decomposition with indels
#assisted by Chat Gpt To translate TRviz algorithm into an R implementation then
#debugged heavily to allow for proper decomposition
decomposeTR <- function(allele, motifs, match = 1, mis_match = -1,
                        indel = -2, allowence = 0) {
  if (is.character(allele)) allele <- DNAString(allele)
  if (is.character(motifs)) motifs <- DNAStringSet(motifs)

  L <- length(allele)
  M <- length(motifs)

  DP <- rep(-Inf, L + 1)
  traceback <- rep(NA_character_, L + 1)
  DP[1] <- 0

  # Precompute motif hits
  hits_df <- do.call(rbind, lapply(seq_along(motifs), function(m) {
    hits <- matchPattern(motifs[[m]], allele, max.mismatch = allowence)
    if(length(hits) == 0) return(NULL)
    data.frame(
      start = start(hits),
      end = end(hits),
      motif_idx = m
    )
  }))
  hits_df <- hits_df[order(hits_df$start), ]

  # Fill DP
  for(pos in 1:L) {
    # Indel: gap of 1 base
    if(DP[pos] != -Inf && pos + 1 <= L + 1) {
      score_gap <- DP[pos] + indel
      if(score_gap > DP[pos + 1]) {
        DP[pos + 1] <- score_gap
        traceback[pos + 1] <- "--"
      }
    }

    # Motif hits starting at this position
    hits_here <- hits_df[hits_df$start == pos, , drop = FALSE]
    for(i in seq_len(nrow(hits_here))) {
      s <- hits_here[i, ]
      motif_seq <- as.character(motifs[[s$motif_idx]])
      subseq_allele <- as.character(subseq(allele, s$start, s$end))
      score <- sum(score_match_custom(subseq_allele, motif_seq, match, mis_match))
      new_score <- DP[pos] + score
      if(new_score > DP[s$end + 1]) {
        DP[s$end + 1] <- new_score
        traceback[s$end + 1] <- as.character(s$motif_idx)
      }
    }
  }

  # Reconstruct composition
  composition <- character()
  pos <- L + 1
  while(pos > 1) {
    t <- traceback[pos]
    if(is.na(t)) break

    if(t == "--") {
      # Look forward to next motif
      next_motif_pos <- which(!is.na(traceback[1:(pos-1)]) & traceback[1:(pos-1)] != "--")
      if(length(next_motif_pos) == 0) {
        # Leading gap: use first motif length
        k <- width(motifs[1])
      } else {
        # Gap before next motif
        next_pos <- max(next_motif_pos)
        k <- if(traceback[next_pos] == "--") 1 else width(motifs[as.integer(traceback[next_pos])])
      }
      gap_len <- pos - if(length(next_motif_pos) == 0) 1 else max(next_motif_pos)
      n_units <- ceiling(gap_len / k)
      composition <- c(rep("--", n_units), composition)
      pos <- pos - gap_len
    } else {
      m <- as.integer(t)
      k <- width(motifs[m])
      composition <- c(as.character(motifs[m]), composition)
      pos <- pos - k
    }
  }

  return(composition)
}

score_match_custom <- function(seq1, seq2, match = 1, mis_match = -1) {
  if(length(seq1) != length(seq2)) stop("seq1 and seq2 must be the same length")
  return(sum(ifelse(seq1 == seq2, match, mis_match)))
}





