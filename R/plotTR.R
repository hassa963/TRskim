#' Plot a Tandem Repeat Alignment
#'
#' Visualizes an aligned tandem repeat matrix as a colored tile graph.
#'
#' @param aln_matrix A character matrix returned by alignTRs()
#' @param motif_map A named character vector mapping encoded letters to motifs
#' @param colour_palette Vector of colours (default = default_palette)
#' @param graph_title Plot title
#' @param title_size Font size for title
#' @param axis_size Font size for axis titles
#' @param legend_size Font size for legend
#' @param show_motifs  Whether or not user wants actual motif or encoded motif
#' in their legend
#'
#' @examples
#' \dontrun{
#' #Example
#' motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")
#' encoded_trs <- c("AABA", "ABA", "BAC")
#'
#' alignment <- alignTRs(encoded_trs)
#' plotTR(alignment, motif_map)
#' }
#'
#' @return A ggplot object, or list(plot, legend) if large = TRUE
#'
#' @references OpenAI. ChatGPT (GPT-5) large language model (2025).
#' https://chat.openai.com/
#'
#' Wickham, H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag,
#' New York, 2016.
#'
#' Wickham, H. Reshaping data with the reshape package. J. Stat.
#' Softw. 21, 1–20 (2007).
#'
#' Wilke, C. cowplot: Streamlined Plot Theme and Plot Annotations for ggplot2.
#' R package version 1.2.0, https://CRAN.R-project.org/package=cowplot (2025)
#'
#' @import ggplot2
#' @import reshape2
#' @import cowplot
#' @export plotTR
plotTR <- function(aln_matrix,
                   motif_map,
                   colour_palette = default_palette,
                   graph_title = "Tandem Repeat Alignment",
                   title_size = 10,
                   axis_size = 8,
                   legend_size = 8,
                   show_motifs = TRUE
) {

  ## -------------------------- ##
  ## Input validation
  ## -------------------------- ##
  if (!is.matrix(aln_matrix)){
    stop("aln_matrix must be a matrix.", call. = FALSE)
  }

  if (!is.character(aln_matrix)){
    stop("aln_matrix must contain character motif codes.", call. = FALSE)
  }

  # Add column names if missing
  if (is.null(colnames(aln_matrix)) || any(colnames(aln_matrix) == "")) {
    colnames(aln_matrix) <- paste0("Pos", seq_len(ncol(aln_matrix)))
  }

  ## -------------------------- ##
  ## Data reshaping
  ## -------------------------- ##
  aln_df <- as.data.frame(aln_matrix, stringsAsFactors = FALSE)
  aln_df$TR <- rownames(aln_matrix)

  aln_long <- reshape2::melt(
    aln_df,
    id.vars = "TR",
    variable.name = "Position",
    value.name = "Motif"
  )

  aln_long$Position <- as.integer(gsub("\\D", "", aln_long$Position))

  ## -------------------------- ##
  ## Palette + legend labels
  ## -------------------------- ##
  motifs <- unique(aln_long$Motif)
  motif_colors <- setNames(colour_palette[seq_along(motifs)], motifs)
  motif_colors["-"] <- "grey90"

  # Decide what legend labels to show
  raw_labels <- sapply(names(motif_colors), function(raw_label) {
    if (raw_label %in% names(motif_map)) {
      if (show_motifs) {
        motif_map[[raw_label]]   # show full motif sequence
      } else {
        raw_label               # show encoded letter
      }
    } else {
      raw_label  # for "-"
    }
  })



  motif_labels <- wrap_label(raw_labels)

  ## -------------------------- ##
  ## Plot construction
  ## -------------------------- ##
  tr_plot <- ggplot2::ggplot(aln_long, aes(x = Position, y = TR, fill = Motif)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(
      values = motif_colors,
      breaks = names(motif_colors),
      labels = motif_labels,
      drop = FALSE,
      guide = guide_legend(
        keywidth = unit(0.5, "lines"),
        keyheight = unit(0.5, "lines")
      )
    ) +
    ggplot2::labs(
      title = graph_title,
      x = "Position",
      y = "Alleles"
    ) +
    ggplot2::theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = legend_size),
      legend.key.size = unit(0.3, "lines"),
      axis.title.x = element_text(size = axis_size, face = "bold"),
      axis.title.y = element_text(size = axis_size, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5)
    )

  return(tr_plot)
}

# Wrap long labels
wrap_label <- function(raw_labels) {
  sapply(raw_labels, function(raw_label){
    paste(strwrap(raw_label, width = 20), collapse = "\n")
    }
    )
}
