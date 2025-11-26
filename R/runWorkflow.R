#' Running the TRskim Workflow
#'
#' Performs all steps of the TRskim workflow, including decomposition, encoding,
#'  alignment, and visualization of tandem repeats.
#'
#' @param alleles A `DNAStringSet` or character vector of tandem repeat alleles.
#' @param motifs A vector of strings or a `DNAStringSet` representing the
#' expected motifs of the tandem repeat.
#' @param match_score Numeric. Score assigned when a motif matches the sequence.
#' Default is `1`.
#' @param indel Numeric. Score assigned for an insertion or deletion in the
#' sequence. Default is `-1`.
#' @param allowance Integer. Number of mismatches allowed between the motif and
#' the sequence for it to still be encoded as that motif. Default is `0`.
#'
#' @return A list containing:
#' \item{decomposition}{The tandem repeat alleles decomposed into their
#' constituent motifs.}
#' \item{encoded}{The encoded tandem repeats.}
#' \item{motif_map}{Mapping of motifs to single-character symbols used in the
#' encoding.}
#' \item{alignment}{Character matrix representing the aligned tandem repeats.}
#' \item{bar_plot}{Bar plot visualization of motif frequencies at each
#' position.}
#' \item{tile_plot}{Tile plot visualization of motif composition across
#' alleles.}
#'
#' @examples
#' \dontrun{
#' alleles <- c("ATATAT", "CGCG")
#' motifs <- c("AT", "CG")
#' trskim_result <- runWorkflow(alleles, motifs)
#' trskim_result$decomposition
#' trskim_result$encoded
#' trskim_result$bar_plot
#' trskim_result$tile_plot
#' }
#' @importFrom TRskim decomposeTRs encodeTRs plotTR alignTRs
#' @export runWorkflow

runWorkflow <- function(alleles, motifs, match_score = 1,
                         indel = -1, allowance = 0){

  message("Decomposing tandem repeats...")

  ###decomposition
  decomposition <- TRskim::decomposeTRs(alleles, motifs, match_score,
                                        indel, allowance)

  decomposed <- decomposition$compositions

  message("Encoding tandem repeats....")

  ####
  encoding <- TRskim::encodeTRs(decomposed, decomposition$motifs)

  encoded <- encoding$encoded

  motif_map <- encoding$motif_map

  message("Aligning tandem repeats....")

  alignment <- alignTRs(encoded)

  message("Generating Tile plot visualization...")

  tile_plot <- TRskim::plotTR(alignment, motif_map,
                              colour_palette = TRskim::default_palette)

  message("Generating Bar plot visualization...")

  bar_p <- TRskim::plotTR(alignment, motif_map, graph_type = "bar",
                          colour_palette = TRskim::default_palette )

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

#[END]
