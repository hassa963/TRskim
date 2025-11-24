#structure is based off of Tabsets and file upload examples in the shiny gallery

# Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, et al. Shiny: Web
# Application Framework for R — Tabsets example. Shiny Gallery, RStudio (2025).
# https://shiny.posit.co/r/gallery/application-layout/tabsets

# Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, et al. Shiny: Web
# Application Framework for R — File upload example. Shiny Gallery, RStudio
# (2025). https://shiny.posit.co/r/gallery/widgets/file-upload/

library(shiny)

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel("Tandem Repeat Comparitive Analysis and Visualization"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      h4("Input Files"),

      # Input: Select a fasta file ----
      tags$p("Upload a fasta file containing the sequences of your tandem
             repeats:"),

      fileInput("tr_file", "Tandem Repeats",
                multiple = FALSE,
                accept = c(".fa", ".fasta", "application/fasta",
                           "application/x-fasta")),

      # Input: Select motif file ----

      tags$p(
        "Upload a file containing your known motifs in a single column
        (.txt, .tsv, or .csv). Each motif should be on a separate line.
        The file should be tab-separated. NOTE: current version can only support
        66 unique motifs."
        ),

      fileInput( "motif_file", "Motifs",
        multiple = FALSE,
        accept = c(".tsv", ".txt")
      ),

      # br() element to introduce extra vertical spacing ----
      br(),

      h4("Visualization Settings"),

      tags$p("The tile plot visualization displays alleles along the y-axis and
             positions along the x-axis, with colors indicating the underlying
             motifs. The bar plot shows positions on the x-axis and the
             frequency of motifs at each position on the y-axis."),

      # Input: Select the visualization type ----
      radioButtons("plot_type", "Visualization type:",
                   c("Tile Plot" = "tile",
                     "Bar Plot" = "bar"
                     )),

      radioButtons("show_motifs", "Do you want unencoded motifs in the legend?
                   (ex. 'A' -> 'ATAT'):",
                   c("Yes" = TRUE,
                     "No" = FALSE
                   ))

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      HTML(
        "TRskim decomposes tandem repeats into their motif composition, encodes
        alleles, performs motif-aware multiple sequence alignment, and generates
        visualizations in the form of bar plots or tile plots.<br><br>
        This provides motif-level decomposition and visualization of tandem
        repeats, enabling researchers to inspect the internal structure of
        repeat arrays across alleles. This level of resolution is essential for
        studying pathogenic repeat expansions and assessing repeat instability,
        as motif structure influences disease relevance, mutational dynamics,
        and secondary structure formation.<br><br>
        Example files can be found in the inst/extdata folder in the TRskim
        github repository:<br>
        - SORL1_tr_sequence.fa – contains the tandem repeat sequences<br>
        - motifs_SORL1.tsv – contains the motifs"
      ),

      # Output: Tabset w/ plot, decomposition, encoding and alignment ----
      tabsetPanel(
        type = "tabs",
        tabPanel("Plot",
                 p("Visualization of the aligned tandem repeats. Each color
                   represents a unique motif:"),
                 plotOutput("tr_plot")

        ),
        tabPanel("Alignment",
                 p("The aligned encoded sequences.
                   Gaps are represented by '-':"),
                 tableOutput("alignment")

        ),
        tabPanel("Motif Map",
                 p("Mapping between motif sequences and their
                   encoded symbols:"),
                 tableOutput("motif_map")

        )
      ))

    )
  )

# Define server logic for random distribution app ----
server <- function(input, output) {

  # ---- Load FASTA ----
  tr_data <- reactive({
    req(input$tr_file)
    Biostrings::readDNAStringSet(input$tr_file$datapath)
  })

  # ---- Load motifs (.tsv / .txt) ----
  motif_seq <- reactive({
    req(input$motif_file)

    motifs <- readLines(input$motif_file$datapath)
    motifs[nzchar(motifs)] #in case there are empty lines
  })


  # ---- decompose ----
  decomposition <- reactive({
    req(motif_seq())
    req(tr_data())

    decomposeTRs(tr_data(), motif_seq())
  })

  # ---- encode ----
  encoded <- reactive({
    req(decomposition())

    encodeTRs(decomposition()$compositions, decomposition()$motifs)
  })

  # ---- alignment ----

  alignment <- reactive({
    req(encoded())

    alignTRs(encoded()$encoded)
  })

  # ---- Motif map output ----
  output$motif_map <- renderTable({
    req(encoded())
    motif_df <- data.frame(
      Symbol = unname(encoded()$motif_map),
      Sequence = names(encoded()$motif_map),
      stringsAsFactors = FALSE
    )
    motif_df
  })

  output$alignment <- renderTable({
    req(alignment())

    # Convert alignment matrix to data frame
    alignment_df <- as.data.frame(alignment())

    # Change column names to position numbers
    colnames(alignment_df) <- 1:ncol(alignment_df)

    # Add allele names as first column if they exist
    if (!is.null(rownames(alignment_df))) {
      alignment_df <- data.frame(
        Allele = rownames(alignment_df),
        alignment_df,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    alignment_df
  },
  spacing = "xs",
  width = "auto",
  align = "l"
  )

  output$tr_plot <- renderPlot({
    req(alignment())
    req(encoded())

    plotTR(
      alignment(),
      encoded()$motif_map,
      show_motifs = input$show_motifs,
      graph_type = input$plot_type
    )
  })

}


# Create Shiny app ----
shinyApp(ui, server)

#[END]
