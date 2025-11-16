#' Encoding Tandem Repeats
#'
#' A function that encodes decomposed tandem repeats such that each unique
#' motif is a single character symbol. **MAXIMUM OF 66 ENCODINGS**
#'
#' @param decomposed_TRs A character vector containing the decomposed tandem
#' repeats (ex. "-"  "-"  "AC" "AC" "GT" "AC")
#' @param motifs  A character vector containing the motifs that may be in the TR
#'
#'@returns
#' Returns encoded tandem repeats in $encoded and the
#' mapping of the motifs to their symbols in $map
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
#' $map
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
#' encoded$map
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

encodeTRs <- function(decomposed_TRs, motifs){

  symbols <- c(LETTERS, letters, 0:9, "!", "@", "#", "$")
  motif_map <- setNames(symbols[seq_along(motifs)], motifs)

  # Wrap single TR in a list for uniform processing
  is_single <- is.character(decomposed_TRs) && !is.list(decomposed_TRs)
  tr_list <- if(is_single) list(decomposed_TRs) else decomposed_TRs

  # Encode each TR
  encoded <- lapply(tr_list, function(tr) paste0(motif_map[tr], collapse = ""))

  # Return single string if input was a single TR
  if(is_single) encoded <- encoded[[1]]

  return(list(encoded = encoded, map = motif_map))
}
