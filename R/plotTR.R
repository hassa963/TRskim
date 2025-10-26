#' Plot a Tandem Repeat Alignment
#'
#' Visualizes an aligned tandem repeat matrix as a colored tile graph.
#'
#' @param aln_matrix A character matrix returned by alignTRs
#' @param motif_map A named character vector where the names are different
#' motifs
#'
#' @return a tile graph showing the motif compositions/alignments of the alleles
#' containg a legend for the different motifs and grew for the indels ("-").
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
#' @import ggplot2
#' @import reshape2
#' @import RColorBrewer
#' @export plotTR

# Plot
#the dataframe formating was assisted by chatgpt
plotTR <- function(aln_matrix, motif_map) {
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

    # Define color palette
    motifs <- unique(df_long$Motif)
    palette <- RColorBrewer::brewer.pal(max(3, min(length(motifs), 8)), "Set2")
    motif_colors <- setNames(rep(palette, length.out = length(motifs)), motifs)
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
    ggplot(df_long, aes(x = Position, y = TR, fill = Motif,
                        label = Motif)) +
      geom_tile(color = "white", size = 0.4) +
      scale_fill_manual(
        values = motif_colors,
        breaks = names(motif_colors),
        labels = motif_labels,
        drop = FALSE
      ) +
      scale_y_discrete(limits = rev) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      ) +
      labs(title = "Tandem Repeat Alignment", x = "Position", y = "Alleles")

  }

