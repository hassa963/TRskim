#' Plot a Tandem Repeat Alignment
#'
#' Visualizes an aligned tandem repeat matrix as a colored tile plot or bar
#' plot.
#'
#' @param aln_matrix A character matrix returned by `alignTRs()`. Each row
#' typically represents an allele, and each column represents a motif position.
#' @param motif_map A **fully named** character vector mapping encoded letters
#' to motifs.
#' @param colour_palette A vector of colors to use for plotting
#' (default = `default_palette`).
#' @param graph_title A character string specifying the plot title.
#' (default = "Tandem Repeat Alignment")
#' @param title_size Numeric value specifying font size for the title.
#' (default = `10`)
#' @param axis_size Numeric value specifying font size for axis labels.
#' (default = `8`)
#' @param legend_size Numeric value specifying font size for the legend.
#' (default = `8`)
#' @param show_motifs Logical. Whether to display the actual motif or the
#' encoded motif in the legend. (default = `TRUE`)
#' @param graph_type Character. Type of plot to generate: `"tile"` for a tile
#' plot showing motif composition of all alleles, or `"bar"` for a bar plot
#' showing motif frequency at each position. (default = `tile`)
#'
#' @return A ggplot object representing the tandem repeat alignment.
#' @examples
#' \dontrun{
#'
#' motif_map <- c(GCTT = "A", TCC = "B", GTG = "C")
#' encoded_trs <- c("AABA", "ABA", "BAC")
#'
#' alignment <- alignTRs(encoded_trs)
#'
#' plotTR(alignment, motif_map)
#'
#' plotTR(alignment, motif_map, graph_type = "bar")
#' }
#'
#' @references OpenAI. ChatGPT (GPT-5) large language model (2025).
#' https://chat.openai.com/
#'
#' Wickham, H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag,
#' New York, 2016.
#'
#' Wickham, H. Reshaping data with the reshape package. J. Stat.
#' Softw. 21, 1–20 (2007).
#' @import ggplot2
#' @import reshape2
#' @export plotTR
plotTR <- function(aln_matrix,
                   motif_map,
                   colour_palette = TRskim::default_palette,
                   graph_title = "Tandem Repeat Alignment",
                   title_size = 10,
                   axis_size = 8,
                   legend_size = 8,
                   show_motifs = TRUE,
                   graph_type = "tile") {

  ## -------------------------- ##
  ## Input validation
  ## -------------------------- ##
  if (!is.matrix(aln_matrix)) {
    stop("aln_matrix must be a matrix.", call. = FALSE)
  }

  if (!is.character(aln_matrix)) {
    stop("aln_matrix must contain character motif codes.", call. = FALSE)
  }

  if (graph_type != "tile" & graph_type != "bar") {
    stop("Graph type must be one of (tile, bar)", call. = FALSE)
  }

  if (!is.vector(motif_map) || is.null(names(motif_map))) {
    stop("motif_map must be named vector", call. = FALSE)
  }

  if (any(is.na(names(motif_map)) | names(motif_map) == "")) {
    stop("Vector is partially named. All elements must have names.")
  }

  # Add column names if missing
  if (is.null(colnames(aln_matrix)) || any(colnames(aln_matrix) == "")) {
    colnames(aln_matrix) <- paste0("Pos", seq_len(ncol(aln_matrix)))
  }

  ## -------------------------- ##
  ## Create reverse lookup: symbol -> motif
  ## motif_map is c(AC = "A", GT = "B")
  ## reverse_map becomes c(A = "AC", B = "GT")
  ## -------------------------- ##
  reverse_map <- setNames(names(motif_map), motif_map)

  if (graph_type == "tile") {
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
    motifs <- as.character(unique(aln_long$Motif))
    motif_colors <- setNames(colour_palette[seq_along(motifs)], motifs)
    motif_colors["-"] <- "grey90"

    # raw_label is the symbol (A, B, C, -)
    # Use reverse_map to get actual motif sequence
    raw_labels <- sapply(names(motif_colors), function(raw_label) {
      if (raw_label %in% names(reverse_map)) {
        if (show_motifs) {
          reverse_map[[raw_label]]  # returns the actual motif
        } else {
          raw_label
        }
      } else {
        raw_label  # for "-"
      }
    })

    motif_labels <- wrap_label(raw_labels)

    ## -------------------------- ##
    ## Plot construction
    ## -------------------------- ##
    tr_plot <- ggplot2::ggplot(aln_long, aes(x = Position, y = TR, fill = Motif)
                               ) +
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

  if (graph_type == "bar") {
    tr_bar_plot <- bar_plot(
      aln_matrix = aln_matrix,
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

bar_plot <- function(aln_matrix, motif_map, colour_palette, graph_title,
                     title_size, axis_size, legend_size, show_motifs) {

  # Create reverse lookup: symbol -> motif
  reverse_map <- setNames(names(motif_map), motif_map)

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

  ## Convert positions like "V1" -> 1
  count_df$Position <- as.integer(gsub("\\D", "", count_df$Position))

  ## Palette
  motifs <- as.character(unique(count_df$Motif))
  motif_colors <- setNames(colour_palette[seq_along(motifs)], motifs)

  ## LEGEND LABELS - use reverse_map
  raw_labels <- sapply(motifs, function(label) {
    if (label %in% names(reverse_map)) {
      if (show_motifs) {
        reverse_map[[label]]
      } else {
        label
      }
    } else {
      label
    }
  })

  motif_labels <- wrap_label(raw_labels)

  ## Bar plot
  bar_p <- ggplot2::ggplot(
    count_df,
    aes(x = Position, y = Count, fill = Motif)
  ) +
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

# Wrap long labels
wrap_label <- function(raw_labels) {
  sapply(raw_labels, function(raw_label) {
    paste(strwrap(raw_label, width = 20), collapse = "\n")
  })
}
# [END]

bar_plot <- function(aln_matrix, motif_map, colour_palette, graph_title,
                     title_size, axis_size, legend_size, show_motifs) {

  # Create reverse lookup: symbol -> motif
  # motif_map is c("ACGT..." = "A", "GCTT..." = "B")
  # reverse_map becomes c(A = "ACGT...", B = "GCTT...")
  reverse_map <- setNames(names(motif_map), motif_map)

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

  ## Convert positions like "V1" -> 1
  count_df$Position <- as.integer(gsub("\\D", "", count_df$Position))

  ## Palette
  motifs <- as.character(unique(count_df$Motif))
  motif_colors <- setNames(colour_palette[seq_along(motifs)], motifs)

  ## LEGEND LABELS - use reverse_map
  raw_labels <- sapply(motifs, function(label) {
    if (label %in% names(reverse_map)) {
      if (show_motifs) {
        reverse_map[[label]]
      } else {
        label
      }
    } else {
      label
    }
  })

  motif_labels <- wrap_label(raw_labels)

  ## Bar plot
  bar_p <- ggplot2::ggplot(
    count_df,
    aes(x = Position, y = Count, fill = Motif)
  ) +
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
