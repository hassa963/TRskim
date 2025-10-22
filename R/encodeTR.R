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
#' mapping of the motifs to their synbols in $motif_map
#'
#' @examples
#' \dontrun{
#' decomposed_TR <- c("-", "-", "AC", "AC", "GT", "AC")
#' motifs <- c("AC", "GT")
#' encodeTR(decomposed_TR, motifs)
#'
#' ## expected
#' #$encoded
#' #"--AABA"
#' #$map
#' #AC  GT   -
#' #"A" "B" "-"
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
#'@export
#emulates how TRviz encodes TRs into motif units for later alignment
encodeTR <- function(decomposed_TR, motifs){

  symbols <- LETTERS[seq_along(motifs)]
  motif_map <- setNames(symbols,motifs) #from stats

  motif_map["-"] <- "-"

  encoded_seq <- paste0(motif_map[decomposed_TR], collapse = "")

  return(list(encoded = encoded_seq,
              map = motif_map))
}
