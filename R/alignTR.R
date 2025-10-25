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
alignTRs <- function(encoded_trs){
  # Split each string into vector of symbols
  seqs <- lapply(encoded_trs, function(x) strsplit(x, "")[[1]])

  # Initialize aligned sequences with the first
  aligned_seqs <- list(seqs[[1]])

  for(k in 2:length(seqs)){
    # Align the new sequence to the consensus so far
    new_aln <- align_pair(aligned_seqs[[1]], seqs[[k]])

    # Update all previously aligned sequences to include new gaps
    for(i in seq_along(aligned_seqs)){
      aligned_seqs[[i]] <- insert_gaps(aligned_seqs[[i]], new_aln$seq1, aligned_seqs[[1]])
    }

    # Add the new aligned sequence
    aligned_seqs[[k]] <- new_aln$seq2

    # Update the reference for next round
    aligned_seqs[[1]] <- new_aln$seq1
  }

  # Convert to a matrix
  aln_matrix <- do.call(rbind, aligned_seqs)
  rownames(aln_matrix) <- paste0("TR", seq_along(encoded_trs))
  return(aln_matrix)
}

#### helpers ####
score <- function(a, b) {
  if(a == b) return(2)
  if(a == "-" || b == "-") return(-2)
  return(-1)
}

### was given a skeleton to perform multiple sequence alignmentfrom chat GPT and
#then filled it in
align_pair <- function(seq1, seq2) {
  n <- length(seq1)
  m <- length(seq2)
  S <- matrix(0, n+1, m+1)
  traceback <- matrix("", n+1, m+1)

  # Initialize
  S[1,] <- seq(0, -2*m, by=-2)
  S[,1] <- seq(0, -2*n, by=-2)
  traceback[1,] <- "left"
  traceback[,1] <- "up"
  traceback[1,1] <- "done"

  # Fill DP matrix
  for(i in 2:(n+1)) {
    for(j in 2:(m+1)) {
      scores <- c(
        S[i-1,j-1] + score(seq1[i-1], seq2[j-1]), # diag
        S[i-1,j] - 2,  # up (gap in seq2)
        S[i,j-1] - 2   # left (gap in seq1)
      )
      S[i,j] <- max(scores)
      idx <- which.max(scores)
      traceback[i,j] <- c("diag","up","left")[idx]
    }
  }

  # Traceback
  i <- n+1; j <- m+1
  aln1 <- c(); aln2 <- c()
  while(traceback[i,j] != "done") {
    dir <- traceback[i,j]
    if(dir=="diag"){
      aln1 <- c(seq1[i-1], aln1)
      aln2 <- c(seq2[j-1], aln2)
      i <- i-1; j <- j-1
    } else if(dir=="up"){
      aln1 <- c(seq1[i-1], aln1)
      aln2 <- c("-", aln2)
      i <- i-1
    } else if(dir=="left"){
      aln1 <- c("-", aln1)
      aln2 <- c(seq2[j-1], aln2)
      j <- j-1
    }
  }
  return(list(seq1 = aln1, seq2 = aln2))
}

### was given a skeleton to perform multiple sequence alignmentfrom chat GPT and
#then filled it in
# Insert gaps from new reference into old sequence
insert_gaps <- function(seq, new_ref, old_ref){
  new_seq <- c()
  j <- 1
  for(i in seq_along(new_ref)){
    if(new_ref[i] == "-"){
      new_seq <- c(new_seq, "-")
    } else {
      new_seq <- c(new_seq, seq[j])
      j <- j + 1
    }
  }
  return(new_seq)
}


