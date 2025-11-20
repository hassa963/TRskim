#' Encoding Tandem Repeats
#'
#' A function that encodes decomposed tandem repeats such that each unique
#' motif is a single character symbol. **MAXIMUM OF 66 ENCODINGS**
#'
#' @param decomposed_TRs A list or chracter vector containing the decomposed
#' tandem repeat(s) (ex."AC" "AC" "GT" "AC")
#' @param motifs  A character vector containing the motifs that may be
#' in the TR
#'
#'@returns
#' Returns encoded tandem repeats in $encoded either as a list if there is
#' multiple tandem repeats or as a single string if there is only one and the
#' mapping of the motifs to their symbols in $motif_map
#'
#' @examples
#' \dontrun{
#' decomposed_TR1 <- c("TTT", "AC", "AC", "GT", "AC")
#' decomposed_TR2 <- c("AC", "GT", "AC")
#' decomposed_TR3 <- c("CC", "AC", "GT", "AC", "GT")
#' decomposed_TRs <- list(decomposed_TR1, decomposed_TR2, decomposed_TR3)
#' motifs <- c("AC", "GT", "TTT", "CC")
#' encoded <- encodeTRs(decomposed_TRs, motifs)
#' encoded
#' $encoded
#' $encoded[[1]]
#' [1] "CAABA"
#'
#' $encoded[[2]]
#' [1] "ABA"
#'
#' $encoded[[3]]
#' [1] "DABAB"
#'
#'
#' $motif_map
#' AC  GT TTT  CC
#' "A" "B" "C" "D"
#' }
#'
#'\dontrun{
#'
#' decomposed_TR <- c("AC", "AC", "GT", "AC")
#' motifs <- c("AC", "GT")
#' encoded <- encodeTRs(decomposed_TR, motifs)
#'
#' encoded$encoded
#' [1] "AABA"
#'
#' encoded$motif_map
#'  AC  GT
#' "A" "B"
#'
#' }
#' @references
#' OpenAI. ChatGPT (GPT-5) large language model (2025).https://chat.openai.com/
#'
#' Park, J., Kaufman, E., Valdmanis, P. N. & Bafna, V. TRviz: A Python Library
#' for decomposing and Visualizing Tandem Repeat Sequences. Bioinformatics
#' Advances 3, (2023).
#'
#' R Core Team. R: A Language and Environment for Statistical Computing.
#' (R Foundation for Statistical Computing, 2025). Available at:
#' https://www.R-project.org/
#'
#'@importFrom stats setNames
#'@export encodeTRs

#checked chatgpt to figure out best way for encoding
encodeTRs <- function(decomposed_TRs, motifs){
  #------------------------------------------------------------
  # Convert DNAStringSet to character if needed
  #------------------------------------------------------------
  if (inherits(motifs, "DNAStringSet") || inherits(motifs, "XStringSet")) {
    motifs <- as.character(motifs)
  }

  #------------------------------------------------------------
  # Validate inputs first (before any cleaning)
  #------------------------------------------------------------
  if (is.null(decomposed_TRs)) {
    stop("decomposed_TRs cannot be NULL")
  }
  if (is.null(motifs) || length(motifs) == 0) {
    stop("motifs cannot be NULL or empty")
  }
  if (!is.character(motifs) || !all(is_valid_nt_string(motifs))) {
    stop("motifs must be of type character with proper nucleotide characters")
  }

  #------------------------------------------------------------
  # Clean decomposed_TRs if it is a list
  #------------------------------------------------------------
  if (is.list(decomposed_TRs)) {
    # Identify invalid or empty TR entries
    invalid_idx <- which(sapply(decomposed_TRs, function(tr) {
      length(tr) == 0 || all(!nzchar(tr)) || any(!is_valid_nt_string(tr))
    }))
    # Warn the user if any were removed
    if (length(invalid_idx) > 0) {
      warning(sprintf(
        "Removed %d invalid or empty TR entries at positions: %s",
        length(invalid_idx),
        paste(invalid_idx, collapse = ", ")
      ))
      # Keep only valid TR entries
      decomposed_TRs <- decomposed_TRs[-invalid_idx]
    }
    # If everything was removed
    if (length(decomposed_TRs) == 0) {
      stop("All decomposed_TRs entries were empty or invalid.")
    }
  }

  #------------------------------------------------------------
  # Validate decomposed_TRs after cleaning
  #------------------------------------------------------------
  # Check if it's a character vector or list of character vectors
  if (is.list(decomposed_TRs)) {
    # Validate each element in the list
    if (!all(sapply(decomposed_TRs, function(tr) {
      is.character(tr) && all(is_valid_nt_string(tr))
    }))) {
      stop("All elements in decomposed_TRs list must be character vectors with proper nucleotide characters")
    }
  } else if (is.character(decomposed_TRs)) {
    # Single character vector
    if (!all(is_valid_nt_string(decomposed_TRs))) {
      stop("decomposed_TRs must contain proper nucleotide characters")
    }
  } else {
    stop("decomposed_TRs must be a character vector or list of character vectors")
  }

  #------------------------------------------------------------
  # Encode
  #------------------------------------------------------------
  symbols <- c(LETTERS, letters, 0:9, "!", "@", "#", "$")
  if (length(motifs) > length(symbols)){
    stop("Too many motifs to encode (maximum ", length(symbols), ")")
  }
  motif_map <- setNames(symbols[seq_along(motifs)], motifs)

  # Wrap single TR in a list for uniform processing
  is_single <- is.character(decomposed_TRs) && !is.list(decomposed_TRs)
  tr_list <- if(is_single) list(decomposed_TRs) else decomposed_TRs

  # Encode each TR
  encoded <- lapply(tr_list, function(tr) paste0(motif_map[tr], collapse = ""))

  # Return single string if input was a single TR
  if(is_single) encoded <- encoded[[1]]

  return(list(encoded = encoded, motif_map = motif_map))
}
#[END]
