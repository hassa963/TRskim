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
#' @return The function returns a list containing the decompositions
#'  ($decomposition), encodings ($encoded), motif map ($motif_map)
#'  alignment matrix ($alignment), bar plot visualization ($bar_plot) and tile
#'  plot visualization ($tile_plot) of the TRskim package work flow.
#'
#' @export run_workflow

run_workflow <- function(alleles, motifs, match_score = 1,
                         indel = -1, allowance = 0){

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

  message("Generating Bar plot visualization...")

  bar_p <- plotTR(alignment, motif_map, graph_type = "bar")

  message("Workflow complete!")

  return(
    list(
      decomposition = decomposed,
      encoded = encoded,
      alignment = alignment,
      motif_map = motif_map,
      tile_plot = tile_plot,
      bar_plot = bar_p
    )
  )
}
