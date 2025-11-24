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
      # Output: Tabset w/ plot, decomposition, encoding and alignment ----
      tabsetPanel(
        type = "tabs",
        tabPanel(
          "Welcome",
          HTML('
    <div style="max-width: 800px; line-height: 1.6;">

      <h3>Welcome to the Tandem Repeat Comparative Analysis and Visualization
      App</h3>

      <h4>About TRskim</h4>

      <p>
        TRskim decomposes tandem repeats into their motif composition, encodes
        alleles, performs motif-aware multiple sequence alignment, and generates
        visualizations in the form of bar plots or tile plots.
        This enables motif-level inspection of repeat structure across alleles,
        helping researchers detect interruptions, examine heterogeneity, and
        evaluate overall repeat architecture. This level of detail is essential
        for understanding pathogenic repeat expansions, repeat instability, and
        motif-driven mutational dynamics.
      </p>

      <p>
        Example files are available in the <code>inst/extdata</code> folder of
        the TRskim GitHub repository.
        <br><strong>Note:</strong> this example dataset was originally taken
        from the TRviz Python library (Park et al., 2023).
      </p>

      <ul>
        <li><code>SORL1_tr_sequence.fa</code> – tandem repeat sequences</li>
        <li><code>motifs_SORL1.tsv</code> – motif definitions</li>
      </ul>

      <p>
        This app provides an end-to-end workflow for decomposing, encoding,
        aligning, and visualizing tandem repeats using <strong>TRskim</strong>.
        Upload your repeat sequences and motif list, and the app will generate a
        motif-aware alignment and customizable visualizations.
      </p>

      <h4>How to Use This App</h4>

      <ol>
        <li>
          <strong>Upload a FASTA file</strong> containing tandem repeat
          sequences.
          Each sequence is treated as an independent allele.
        </li>

        <li>
          <strong>Upload a motif file</strong> containing one motif per line
          (<code>.txt</code> or <code>.tsv</code>).
          <em>Current version supports up to 66 unique motifs.</em>
        </li>

        <li>
          Choose the <strong>visualization type</strong>:
          <ul>
            <li><strong>Tile Plot</strong> – each allele shown as a row,
            colored by motif.</li>
            <li><strong>Bar Plot</strong> – motif frequency at each aligned
            position.</li>
          </ul>
        </li>

        <li>
          Select whether the legend should show only motif symbols or
          <strong>symbol → motif</strong> mappings (e.g., A → ATAT).
        </li>

        <li>
          Explore your results in the tab panels:
          <ul>
            <li><strong>Plot</strong> – visualization of the aligned
            repeats</li>
            <li><strong>Alignment</strong> – encoded multiple sequence
            alignment</li>
            <li><strong>Motif Map</strong> – motif-to-symbol mapping
            table</li>
          </ul>
        </li>
      </ol>

      <h4>Motif File Requirements</h4>

      <ul>
        <li>One motif per line (no header).</li>
        <li>Motifs must be uppercase and composed of valid nucleotide characters
        only.</li>
        <li>Blank lines will be ignored.</li>
        <li>Up to 66 unique motifs.</li>
      </ul>

      <h4>Notes</h4>

      <ul>
        <li>The decomposition relies on exact matching between motif patterns
        and sequence structure.</li>
        <li>The order of motifs in the file does not affect encoding or
        alignment.</li>
        <li>Gaps introduced during alignment are shown as <code>-</code> in the
        Alignment tab.</li>
      </ul>

    </div>
  ')
        )
        ,
        tabPanel(
          "References",
          HTML('
    <div style="max-width: 800px; line-height: 1.6;">

      <h3>References</h3>

      <ul>
        <li>
          BioRender.com. <em>BioRender</em> [Online]. Available at:
          <a href="https://www.biorender.com"
          target="_blank">https://www.biorender.com</a>
          (accessed 26 October 2025).
        </li>

        <li>
          Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, et al. (2025).
          <em>Shiny: Web Application Framework for R — Tabsets example.</em>
          Shiny Gallery, RStudio.
          <a href="https://shiny.posit.co/r/gallery/application-layout/tabsets"
          target="_blank">
            Link
          </a>
        </li>

        <li>
          Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, et al. (2025).
          <em>Shiny: Web Application Framework for R — File upload example.</em>
          Shiny Gallery, RStudio.
          <a href="https://shiny.posit.co/r/gallery/widgets/file-upload/"
          target="_blank">
            Link
          </a>
        </li>

        <li>
          OpenAI. (2025). ChatGPT (GPT-5) large language model.
          <a href="https://chat.openai.com/" target="_blank">
          https://chat.openai.com/</a>
        </li>

        <li>
          Pagès H, Aboyoun P, Gentleman R & DebRoy S. (2025).
          <em>Biostrings: Efficient manipulation of biological strings</em>
          (R package version 2.77.2).
          <a href="https://bioconductor.org/packages/Biostrings"
          target="_blank">
            Bioconductor
          </a>, doi:10.18129/B9.bioc.Biostrings
        </li>

        <li>
          Park J, Kaufman E, Valdmanis PN & Bafna V. (2023).
          <em>TRviz: A Python Library for decomposing and Visualizing
          Tandem Repeat Sequences.</em>
          Bioinformatics Advances 3.
        </li>

        <li>
          R Core Team. (2025). <em>R: A Language and Environment for Statistical
          Computing.</em>
          Vienna, Austria: R Foundation for Statistical Computing.
          <a href="https://www.R-project.org/" target="_blank">
          https://www.R-project.org/</a>
        </li>

        <li>
          Wickham H. (2016). <em>ggplot2: Elegant Graphics for Data Analysis.
          </em>
          Springer-Verlag, New York.
        </li>

        <li>
          Wickham H. (2007). Reshaping data with the reshape package.
          <em>J. Stat. Softw.</em> 21, 1–20.
        </li>

      </ul>

    </div>
  ')
        ),
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
