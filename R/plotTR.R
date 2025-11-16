#' Plot a Tandem Repeat Alignment
#'
#' Visualizes an aligned tandem repeat matrix as a colored tile graph.
#'
#' @param aln_matrix A character matrix returned by alignTRs
#' @param motif_map A named character vector where the names are different
#' motifs
#' @param large A boolean value that makes the legend separate from the plot if
#' it is a large number of motifs
#'
#' @return a tile graph showing the motif compositions/alignments of the alleles
#' containg a legend for the different motifs and grew for the indels ("-"). If
#' large is true it will return the plot in $plot and an object containing the
#' legend which can be viewed with `grid.newpage()` and `grid.draw()`
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
# Plot
#the dataframe formating was assisted by chatgpt
plotTR <- function(aln_matrix, motif_map, large = FALSE) {
    # Ensure input is a matrix
    aln_matrix <- as.matrix(aln_matrix)

    # Give columns names if missing
    if (is.null(colnames(aln_matrix)) || any(colnames(aln_matrix) == "")) {
      colnames(aln_matrix) <- paste0("Pos", seq_len(ncol(aln_matrix)))
    }

    # Convert to data frame
    df <- as.data.frame(aln_matrix, stringsAsFactors = FALSE)
    df$TR <- rownames(aln_matrix)

    # Melt into long format
    df_long <- reshape2::melt(df,
                              id.vars = "TR",
                              variable.name = "Position",
                              value.name = "Motif")

    # Convert position to numeric (extract digits)
    df_long$Position <- as.integer(gsub("\\D", "", df_long$Position))

    motifs <- unique(df_long$Motif)
    motif_colors <- setNames(grDevices::rainbow(length(motifs)), motifs)
    motif_colors["-"] <- "grey90"

    # Build motif legend labels
    motif_labels <- sapply(names(motif_colors), function(x) {
      if (x %in% motif_map) {
        names(motif_map)[motif_map == x]  # map encoded letter -> actual motif
      } else {
        x  # for gaps "-" or any unmapped symbol
      }
    })

    # Plot
    plot <- ggplot(df_long, aes(x = Position, y = TR,
                                fill = Motif,
                                label = Motif)) +
      geom_tile(color = "white", size = 0.4) +
      scale_fill_manual(
        values = motif_colors,
        breaks = names(motif_colors),
        labels = motif_labels,
        drop = FALSE,
        guide = guide_legend(ncol = 2)
      ) +
      theme_minimal(base_size = 6) +
      theme(
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_text(size = 6, face = "bold"),
        axis.title.y = element_text(size = 6, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
      ) +
      labs(title = "Tandem Repeat Alignment", x = "Position", y = "Alleles")

    if( large ){
      p_no_legend <- plot + theme(legend.position = "none")

      legend_only <- get_legend(plot)

      return(list(plot = p_no_legend, legend = legend_only))
    }

    else{
      return(plot)
    }


  }

DEFAULT <- c(
  "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0",
  "#f032e6","#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#800000",
  "#aaffc3", "#808000","#ffd8b1", "#000075", "#ff4500", "#32cd32", "#1e90ff",
  "#ff1493", "#00ced1", "#9b0000","#ff69b4", "#7cfc00", "#ff8c00", "#9400d3",
  "#00bfff", "#adff2f", "#ff6347", "#40e0d0","#ff7f50", "#6495ed", "#dc143c",
  "#ff00ff", "#228b22", "#dda0dd", "#ffa500", "#5f9ea0","#d2691e", "#ba55d3",
  "#00fa9a", "#ffdead", "#4b0082", "#ffb6c1", "#32cd32", "#ff4500","#7b68ee",
  "#00ff7f", "#ff6347", "#20b2aa", "#ff69b4", "#ff8c00", "#9acd32", "#ff1493",
  "#4682b4", "#ff6347", "#6a5acd", "#00ff7f", "#ff4500", "#00ced1", "#ff69b4",
  "#32cd32","#ff8c00", "#7cfc00"
)

