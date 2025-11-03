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
#' @return This function returns the tandem repeat allele decomposed into its
#' motifs as a character vector. For example, the allele `"ATATAT"` with motif
#' `"AT"` is decomposed into `c("AT", "AT", "AT")` as well as an updated motifs
#' DNAStringSet consisting of novel motifs detected in the allele.
#'
#' If no motifs are found in the allele from the initial motif set,
#' the function returns `character(0)`and prints a message to inform the user.
#'
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
#' ###Expected:
#' $composition
#' [1] "TTT" "AC"  "AC"  "GT"  "AC"
#'
#' $motifs
#' DNAStringSet object of length 3:
#'   width seq
#' [1]     2 AC
#' [2]     2 GT
#' [3]     3 TTT
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
#' ###Expected
#' $composition
#' [1] "TTT" "AC" "AC" "GTC" "AC"
#' $motifs
#' DNAStringSet object of length 3:
#'   width seq
#' [1]     2 AC
#' [2]     3 GTC
#' [3]     3 TTT
#'}
#' @references OpenAI. ChatGPT (GPT-5) large language model (2025).
#' https://chat.openai.com/
#'
#' Pagès, H., Aboyoun, P., Gentleman, R. & DebRoy, S. Biostrings: Efficient
#' manipulation of biological strings (R package version 2.77.2, 2025).
#' https://bioconductor.org/packages/Biostrings, doi:10.18129/B9.bioc.Biostrings
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
decomposeTR <- function(allele, motifs, match = 1, indel = -1, allowence = 0) {
  # Ensure correct object types
  if (is.character(allele)) allele <- DNAString(allele)
  if (is.character(motifs)) motifs <- DNAStringSet(motifs)

  if(!is(allele, "DNAString")) stop("allele must be DNAString or character")
  if(!is(motifs, "DNAStringSet")) stop("motifs must be DNAStringSet or character vector")

  L <- length(allele)
  M <- length(motifs)

  DP <- rep(-Inf, L + 1)
  traceback <- rep(NA_character_, L + 1)
  DP[1] <- 0

  # Precompute motif hits
  hits_df <- do.call(rbind, lapply(seq_along(motifs), function(m) {
    hits <- matchPattern(motifs[[m]], allele, max.mismatch = allowence)
    if (length(hits) == 0) return(NULL)

    df <- data.frame(
      start = start(hits),
      end = end(hits),
      motif_idx = m
    )

    # Filter invalid coordinates
    df <- df[df$start >= 1 & df$end <= length(allele) & df$start <= df$end, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)

    df
  }))


  if (is.null(hits_df) || nrow(hits_df) == 0) {
    message("No motif matches found in allele.")
    return(character(0))
  }

  hits_df <- hits_df[order(hits_df$start), ]

  # Fill DP table
  for(pos in 1:L) {
    # Gap/indel
    if(DP[pos] != -Inf && pos + 1 <= L + 1) {
      score_gap <- DP[pos] + indel
      if(score_gap > DP[pos + 1]) {
        DP[pos + 1] <- score_gap
        traceback[pos + 1] <- "-"
      }
    }

    # Motif hits
    hits_here <- hits_df[hits_df$start == pos, , drop = FALSE]
    for(i in seq_len(nrow(hits_here))) {
      s <- hits_here[i, ]
      motif_seq <- as.character(motifs[[s$motif_idx]])
      subseq_allele <- as.character(subseq(allele, s$start, s$end))
      score <- if(identical(subseq_allele, motif_seq)) match else 0
      new_score <- DP[pos] + score
      if(new_score > DP[s$end + 1]) {
        DP[s$end + 1] <- new_score
        traceback[s$end + 1] <- as.character(s$motif_idx)
      }
    }
  }

  # Reconstruct sequence with insertions preserved
  composition <- reconstruct(L, traceback, motifs, allele)

  new_motifs <- setdiff(composition, motifs)
  motifs <- c(motifs, DNAStringSet(new_motifs))

  return(list(
    composition = composition,
    motifs = motifs
  ))
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





