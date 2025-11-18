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
#' @param graph_type The user can select the type of visualization they want
#' for the plot so that it can either be a tile plot using the argument "tile"
#' or a bar plot using the argument "bar"
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
                   show_motifs = TRUE,
                   graph_type = "tile"
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

  if(graph_type != "tile" & graph_type != "bar" ){
    stop("Graph type must be one of (tile, bar)", call. = FALSE)
  }

  # Add column names if missing
  if (is.null(colnames(aln_matrix)) || any(colnames(aln_matrix) == "")) {
    colnames(aln_matrix) <- paste0("Pos", seq_len(ncol(aln_matrix)))
  }

  if(graph_type == "tile"){
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

  return(tr_plot) }

  if(graph_type == "bar"){

    tr_bar_plot <- bar_plot( aln_matrix = aln_matrix,
      motif_map = motif_map,
      colour_palette = colour_palette,
      graph_title = graph_title,
      title_size = title_size,
      axis_size = axis_size,
      legend_size = legend_size,
      show_motifs = show_motifs
    )

    return(tr_bar_plot)

  }
}

# Wrap long labels
wrap_label <- function(raw_labels) {
  sapply(raw_labels, function(raw_label){
    paste(strwrap(raw_label, width = 20), collapse = "\n")
    }
    )
}

bar_plot <- function( aln_matrix, motif_map, colour_palette, graph_title,
                      title_size, axis_size,legend_size, show_motifs)
  {
  ## Reshape
  aln_df <- as.data.frame(aln_matrix, stringsAsFactors = FALSE)
  aln_df$TR <- rownames(aln_matrix)

  aln_long <- reshape2::melt(
    aln_df,
    id.vars = "TR",
    variable.name = "Position",
    value.name = "Motif"
  )

  ## Remove gaps
  aln_long <- aln_long[aln_long$Motif != "-", ]

  ## Count motifs per position
  count_df <- as.data.frame(table(aln_long$Position, aln_long$Motif))
  colnames(count_df) <- c("Position", "Motif", "Count")

  ## Convert positions like "V1" → 1
  count_df$Position <- as.integer(gsub("\\D", "", count_df$Position))

  ## Palette
  motifs <- unique(count_df$Motif)
  motif_colors <- setNames(colour_palette[seq_along(motifs)], motifs)

  ## LEGEND LABELS
  motif_labels <- sapply(motifs, function(label) {
    if (show_motifs && label %in% names(motif_map)) {
      motif_map[[label]]        # actual sequence
    } else {
      label                     # encoded symbol
    }
  })

  ## Bar plot
  bar_p <- ggplot2::ggplot(count_df,
                       aes(x = Position, y = Count, fill = Motif)) +
    ggplot2::geom_col(position = "stack", color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = motif_colors,
      breaks = motifs,
      labels = motif_labels,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = graph_title,
      x = "Position",
      y = "Motif Count"
    ) +
    ggplot2::theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = legend_size),
      axis.title.x = element_text(size = axis_size, face = "bold"),
      axis.title.y = element_text(size = axis_size, face = "bold"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5)
    )

  return(bar_p)
}

#[END]
