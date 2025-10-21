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
#' this is 1
#' @import Biostrings
#' @export decomposeTR
# DP decomposition with indels

decomposeTR <- function(allele, motifs, match = 1, mis_match = -1, indel = -2, allowence = 1) {
  if (is.character(allele)) allele <- DNAString(allele)
  if (is.character(motifs)) motifs <- DNAStringSet(motifs)

  L <- length(allele)
  M <- length(motifs)

  # DP[i] = best score for prefix ending at i
  DP <- rep(-Inf, L + 1)
  traceback <- rep("", L + 1)     # motif index or "--"
  DP[1] <- 0                      # Base case: position 0

  # Flatten all hits into a data.frame: start, end, motif_idx
  hits_df <- do.call(rbind, lapply(seq_along(motifs), function(m) {
    hits <- matchPattern(motifs[[m]], allele, max.mismatch = allowence)
    if(length(hits) == 0) return(NULL)
    data.frame(
      start = start(hits),
      end = end(hits),
      motif_idx = m
    )
  }))
  hits_df <- hits_df[order(hits_df$end), ]

  # Fill DP
  for(i in seq_len(nrow(hits_df))) {
    s <- hits_df[i, ]
    motif_seq <- as.character(motifs[[s$motif_idx]])
    subseq_allele <- as.character(subseq(allele, s$start, s$end))
    score <- sum(score_match_custom(subseq_allele, motif_seq, match, mis_match))

    # Previous best score
    prev_score_noindel <- DP[s$start]           # from previous motif
    prev_score_indel <- max(DP[s$start:(s$end)]) + indel  # consider gaps
    if(prev_score_noindel >= prev_score_indel) {
      prev_score <- prev_score_noindel
      from_indel <- FALSE
    } else {
      prev_score <- prev_score_indel
      from_indel <- TRUE
    }

    new_score <- prev_score + score

    if(new_score > DP[s$end + 1]) {
      DP[s$end + 1] <- new_score
      traceback[s$end + 1] <- if(from_indel) "--" else as.character(s$motif_idx)
    }
  }

  # Reconstruct composition
  composition <- character()
  pos <- L + 1
  while(pos > 1) {
    t <- traceback[pos]
    if(t == "") break

    if(t == "--") {
      composition <- c("--", composition)
      pos <- pos - 1
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




