#' Tandem Repeat Encoder
#'
#' A function that encodes a decomposed tandem repeat such that each unique
#' motif is a single character symbol
#'
#' @param decomposed_TR A character vector containing the decomposed tandem
#' repeat (ex. "-"  "-"  "AC" "AC" "GT" "AC")
#' @param motifs  A character vector containing the motifs that may be in the TR
#'
#'@returns
#' Returns a list containing the encoded tandem repeat in $encoded and the
#' mapping of the motifs to their symbols in $motif_map
#'
#' @examples
#' \dontrun{
#'
#' decomposed_TR <- c("AC", "AC", "GT", "AC")
#' motifs <- c("AC", "GT")
#' encoded <- encodeTR(decomposed_TR, motifs)
#'
#' encoded$encoded
#' [1] "AABA"
#'
#' encoded$map
#'  AC  GT
#' "A" "B"
#'
#' }
#'
#' @references
#' Park, J., Kaufman, E., Valdmanis, P. N. & Bafna, V. TRviz: A Python Library
#' for decomposing and Visualizing Tandem Repeat Sequences. Bioinformatics
#' Advances 3, (2023).
#'
#' R Core Team. R: A Language and Environment for Statistical Computing.
#' (R Foundation for Statistical Computing, 2025). Available at:
#' https://www.R-project.org/
#'@importFrom stats setNames
#'@export encodeTR
#emulates how TRviz encodes TRs into motif units for later alignment
#some assistance with AI to create named chacrer vector
encodeTR <- function(decomposed_TR, motifs){
  motif_chars <- as.character(motifs)
  symbols <- c(LETTERS, letters,
               0:9, "!", "@", "#", "$", "%", "&", "*", "+", "-", "=")
  motif_map <- setNames(all_symbols[seq_along(motif_chars)], motif_chars)
  encoded_seq <- paste0(motif_map[decomposed_TR], collapse = "")

  return(list(encoded = encoded_seq,
              map = motif_map))
}

#' Encoding Tandem Repeats
#'
#' A function that encodes decomposed tandem repeats such that each unique
#' motif is a single character symbol
#'
#' @param decomposed_TRs A character vector containing the decomposed tandem
#' repeats (ex. "-"  "-"  "AC" "AC" "GT" "AC")
#' @param motifs  A character vector containing the motifs that may be in the TR
#'
#'@returns
#' Returns a list containing the encoded tandem repeats in $encoded and the
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
#' @references
#' OpenAI. ChatGPT (GPT-5) large language model (2025).https://chat.openai.com/
#'
#' Park, J., Kaufman, E., Valdmanis, P. N. & Bafna, V. TRviz: A Python Library
#' for decomposing and Visualizing Tandem Repeat Sequences. Bioinformatics
#' Advances 3, (2023).
#'
#'@export encodeTRs

encodeTRs <- function(decomposed_TRs, motifs) {

  # Create mapping from motif -> symbol
  #used ChatGpt to figure out how to expand motif alphabet
  motif_chars <- as.character(motifs)
  symbols <- c(LETTERS, letters,
               0:9, "!", "@", "#", "$", "%", "&", "*", "+", "-", "=")
  motif_map <- setNames(all_symbols[seq_along(motif_chars)], motif_chars)



  # Encode each decomposed TR as a single string in a list
  encoded <- lapply(decomposed_TRs, function(tr) {
    paste0(motif_map[tr], collapse = "")
  })

  return(list(
    encoded = encoded,
    map = motif_map
  ))
}
