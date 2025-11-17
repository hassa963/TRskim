#' Running The TRskim Workflow
#'
#' A function that performs all the steps of the TRskim workflow
#'
#' @param alleles DNAStringSet or character string vector of the Tandem repeats
#' @param motifs a vector of strings or a DNAStringSet that are the expected
#' motifs of the Tandem Repeat
#' @param match_score score for if motif matches sequence this is by default 1
#' @param indel score for if there is an insertion or deletion in sequence by
#' default this is -1
#' @param allowance is the number of mismatches between the motif
#' and the sequence permitted to be still be later encoded as a motif. By
#' default this is 0
#' @references
#' Pagès, H., Aboyoun, P., Gentleman, R. & DebRoy, S. Biostrings: Efficient
#' manipulation of biological strings (R package version 2.77.2, 2025).
#' https://bioconductor.org/packages/Biostrings, doi:10.18129/B9.bioc.Biostrings
#'
#' @return The function returns a list containing the decompositions
#'  ($decomposition), encodings ($encoded), motif map ($motif_map)
#'  alignment matrix ($alignment), and tile plot visualization ($tile_plot) of
#'  the TRskim package work flow.
#'
#' @export run_workflow

run_workflow <- function(alleles, motifs, match_score = 1,
                         indel = -1, allowance = 0){

  #------------------------------------------------------------
  # Normalize alleles
  #------------------------------------------------------------
  #used classes from biostrings
  if (is.character(alleles)) {
    if (length(alleles) == 1) {
      alleles <- Biostrings::DNAString(alleles)      # single string
    } else {
      alleles <- Biostrings::DNAStringSet(alleles)   # vector of strings
    }
  }

  if (is.character(motifs)) {
    motifs <- Biostrings::DNAStringSet(motifs)
  }

  # Validate inputs
  if (!inherits(alleles, "DNAString") &&
     !inherits(alleles, "DNAStringSet")) {
   stop(
    "Alleles must be a character string/vector, a DNAString, or a DNAStringSet")
  }

  if (!inherits(motifs, "DNAStringSet")) {
    stop("Motifs must be a character vector or a DNAStringSet")
  }

  message("Decomposing tandem repeats...")

  ###decomposition
  decomposition <- decomposeTRs(alleles, motifs, match_score, indel, allowance)

  decomposed <- decomposition$compositions

  message("Encoding tandem repeats....")

  ####
  encoding <- encodeTRs(decomposed, decomposition$motifs)

  encoded <- encoding$encoded

  motif_map <- encoding$motif_map

  message("Aligning tandem repeats....")

  alignment <- alignTRs(encoded)

  message("Generating Tile plot visualization...")

  tile_plot <- plotTR(alignment, motif_map)

  message("Workflow complete!")

  return(
    list(
      decomposition = decomposed,
      encoded = encoded,
      alignment = alignment,
      motif_map = motif_map,
      tile_plot = tile_plot
    )
  )
}
