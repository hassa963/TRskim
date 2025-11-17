#' Aligning Tandem Repeats
#'
#' Aligns multiple encoded tandem repeat sequences by motif using a simple
#' progressive pairwise alignment algorithm
#'
#' @param encoded_trs A vector of encoded tandem repeats
#'
#' @returns A character matrix representing the aligned tandem repeats, where
#'   each row corresponds to an input sequence and each column to an aligned
#'   motif position.
#' @details
#' **Algorithm:**
#' The function uses progressive alignment with the following strategy:
#' 1. Sequences are sorted by length (longest first)
#' 2. The longest sequence initializes the alignment profile
#' 3. Shorter sequences are progressively aligned to the growing profile
#' 4. When gaps are inserted into the profile, all previously aligned sequences
#'    are updated to maintain column consistency
#' 5. Leading gap-only columns are removed from the final alignment
#'
#' **Alignment Method:**
#' Pairwise alignments use the Needleman-Wunsch algorithm with:
#' - Match score: +2
#' - Mismatch penalty: -1
#' - Gap penalty: -2
#'
#'  **Assumptions and Limitations:**
#' - **Length-based ordering**: Assumes longer sequences are more "complete" and
#'   should anchor the alignment. This works well for TRs
#'   where length differences reflect copy number variation.
#' - **Order-independent**: Results depend only on sequence content and lengths,
#'   not input order.
#' - **Not phylogeny-aware**: Does not consider evolutionary relationships
#' between sequences.
#' - **Best for similar sequences**: Works well when sequences share a common
#' motif structure with varying copy numbers. May perform poorly on highly
#' divergent sequences.
#' - **Small to medium datasets**: Designed for typical TR locus analysis
#'   (< 100 sequences). For very large datasets, consider specialized MSA tools.
#'
#' @examples
#' \dontrun{
#' #Example 1
#' encoded_trs <- c("AABA", "ABA", "BA")
#' alignTRs(encoded_trs)
#'
#' ###expected
#' #' [,1] [,2] [,3] [,4]
#' TR1 "A"  "A"  "B"  "A"
#' TR2 "-"  "A"  "B"  "A"
#' TR3 "-"  "-"  "B"  "A"
#' }
#'
#' @references OpenAI. ChatGPT (GPT-5) large language model (2025).
#' https://chat.openai.com/
#'
#'
#' @export alignTRs

### was given a skeleton to perform multiple sequence alignmentfrom by chat GPT
#and then filled it in
#### Progressive Alignment by Length ####

alignTRs <- function(encoded_trs) {
  seqs <- lapply(encoded_trs, function(x) strsplit(x, "")[[1]])

  # Order by length (longest first) - but keep track of original order
  lengths <- sapply(seqs, length)
  order_idx <- order(lengths, decreasing = TRUE)

  # Start with longest sequence as profile
  profile <- seqs[[order_idx[1]]]
  aligned <- vector("list", length(seqs))
  aligned[[order_idx[1]]] <- profile

  # Progressively align shorter sequences
  for (i in 2:length(order_idx)) {
    idx <- order_idx[i]

    # Align this sequence to current profile
    aln <- align_pair(profile, seqs[[idx]])

    # Update all previously aligned sequences with new gaps
    old_profile <- profile
    profile <- aln$seq1

    for (j in 1:(i-1)) {
      prev_idx <- order_idx[j]
      aligned[[prev_idx]] <- insert_gaps(aligned[[prev_idx]], old_profile,
                                         profile)
    }

    # Store newly aligned sequence
    aligned[[idx]] <- aln$seq2
  }

  # Convert to matrix
  aln_matrix <- do.call(rbind, aligned)
  rownames(aln_matrix) <- paste0("TR", seq_along(encoded_trs))

  # Remove leading all-gap columns
  has_content <- apply(aln_matrix, 2, function(col) any(col != "-"))
  if (any(has_content)) {
    first <- which(has_content)[1]
    aln_matrix <- aln_matrix[, first:ncol(aln_matrix), drop = FALSE]
  }

  return(aln_matrix)
}

#### Insert gaps based on profile change ####
insert_gaps <- function(seq, old_profile, new_profile) {

  # Quick check: if the profile lengths are the same, no gaps were added
  if (length(old_profile) == length(new_profile)) {
    return(seq)
  }

  # Initialize a new sequence vector with the length of the updated profile
  new_seq <- character(length(new_profile))
  seq_idx <- 1 # Position in the old aligned sequence
  old_idx <- 1 # Position in the old profile

  # Loop through positions in the new profile
  for (new_idx in seq_along(new_profile)) {

    if (old_idx <= length(old_profile) && new_profile[new_idx] ==
        old_profile[old_idx]) {

      # Case 1: Position matches old profile
      # Copy the corresponding element from the old sequence
      if (seq_idx <= length(seq)) {
        new_seq[new_idx] <- seq[seq_idx]
        seq_idx <- seq_idx + 1
      } else {
        # If old sequence is exhausted, insert a gap
        new_seq[new_idx] <- "-"
      }
      old_idx <- old_idx + 1 # Move to next position in old profile

    } else if (new_profile[new_idx] == "-") {
      # Case 2: New gap was inserted in the profile
      # Insert a gap in the aligned sequence
      new_seq[new_idx] <- "-"

    } else {
      # Case 3: New profile position doesn't match old profile and isn't a gap
      # This usually shouldn’t happen
      if (seq_idx <= length(seq)) {
        new_seq[new_idx] <- seq[seq_idx] # Copy next element from old sequence
        seq_idx <- seq_idx + 1
      } else {
        # If sequence is exhausted, insert a gap
        new_seq[new_idx] <- "-"
      }
    }
  }

  # Return the updated sequence aligned to the new profile
  return(new_seq)
}

#### Needleman-Wunsch alignment ####
#explained to me by Chat Gpt and broken into implementation steps
align_pair <- function(seq1, seq2) {
  len_seq1 <- length(seq1)   # Length of first sequence
  len_seq2 <- length(seq2)   # Length of second sequence

  # --- Scoring parameters ---
  match_score <- 2
  mismatch_score <- -1
  gap_penalty <- -2

  # --- Initialize DP (dynamic programming) matrix ---
  dp_matrix <- matrix(-Inf, len_seq1 + 1, len_seq2 + 1)
  dp_matrix[1, 1] <- 0  # Start point

  # Initialize first column (seq2 gaps)
  for (i in 2:(len_seq1 + 1)){
    dp_matrix[i, 1] <- dp_matrix[i - 1, 1] + gap_penalty
  }
  # Initialize first row (seq1 gaps)
  for (j in 2:(len_seq2 + 1)){
     dp_matrix[1, j] <- dp_matrix[1, j - 1] + gap_penalty
  }
  # --- Fill DP matrix ---
  for (i in 2:(len_seq1 + 1)) {
    for (j in 2:(len_seq2 + 1)) {
      current_match <- ifelse(seq1[i - 1] == seq2[j - 1],
                              match_score, mismatch_score)

      # Choose the best score among diagonal, up, and left moves
      dp_matrix[i, j] <- max(
        dp_matrix[i - 1, j - 1] + current_match,  # diagonal
        dp_matrix[i - 1, j] + gap_penalty,       # up (gap in seq2)
        dp_matrix[i, j - 1] + gap_penalty        # left (gap in seq1)
      )
    }
  }

  # --- Traceback to recover alignment ---
  i <- len_seq1 + 1
  j <- len_seq2 + 1
  aligned_seq1 <- character(0)
  aligned_seq2 <- character(0)

  while (i > 1 || j > 1) {
    if (i > 1 && j > 1) {
      current_match <- ifelse(seq1[i - 1] == seq2[j - 1], match_score,
                              mismatch_score)
      # Check if current cell came from diagonal
      if (abs(dp_matrix[i, j] - (dp_matrix[i - 1, j - 1] +
                                 current_match)) < 1e-10) {
        aligned_seq1 <- c(seq1[i - 1], aligned_seq1)
        aligned_seq2 <- c(seq2[j - 1], aligned_seq2)
        i <- i - 1
        j <- j - 1
        next
      }
    }

    # Check if current cell came from up (gap in seq2)
    if (i > 1 && abs(dp_matrix[i, j] - (dp_matrix[i - 1, j] + gap_penalty))
        < 1e-10) {
      aligned_seq1 <- c(seq1[i - 1], aligned_seq1)
      aligned_seq2 <- c("-", aligned_seq2)
      i <- i - 1
    } else if (j > 1) {
      # Otherwise, must have come from left (gap in seq1)
      aligned_seq1 <- c("-", aligned_seq1)
      aligned_seq2 <- c(seq2[j - 1], aligned_seq2)
      j <- j - 1
    }
  }

  # Return aligned sequences as a list
  return(list(seq1 = aligned_seq1, seq2 = aligned_seq2))
}
