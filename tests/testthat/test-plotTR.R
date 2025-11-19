#OpenAI. ChatGPT (GPT-5) large language model (2025).
# https://chat.openai.com/

#supported by AI to write tests

pdf(NULL) #so plots don't go to pdf

# Test suite for plotTR error handling
test_that("plotTR validates aln_matrix is a matrix", {
  motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")

  # Test with data frame
  aln_df <- data.frame(
    Pos1 = c("A", "B"),
    Pos2 = c("B", "A")
  )
  expect_error(
    plotTR(aln_df, motif_map),
    "aln_matrix must be a matrix.",
    fixed = TRUE
  )

  # Test with vector
  aln_vec <- c("A", "B", "C")
  expect_error(
    plotTR(aln_vec, motif_map),
    "aln_matrix must be a matrix.",
    fixed = TRUE
  )

  # Test with list
  aln_list <- list(c("A", "B"), c("C", "A"))
  expect_error(
    plotTR(aln_list, motif_map),
    "aln_matrix must be a matrix.",
    fixed = TRUE
  )
})

test_that("plotTR validates aln_matrix contains character data", {
  motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")

  # Test with numeric matrix
  aln_numeric <- matrix(c(1, 2, 3, 4), nrow = 2)
  expect_error(
    plotTR(aln_numeric, motif_map),
    "aln_matrix must contain character motif codes.",
    fixed = TRUE
  )

  # Test with logical matrix
  aln_logical <- matrix(c(TRUE, FALSE, TRUE, FALSE), nrow = 2)
  expect_error(
    plotTR(aln_logical, motif_map),
    "aln_matrix must contain character motif codes.",
    fixed = TRUE
  )
})

test_that("plotTR validates graph_type parameter", {
  motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")
  aln_matrix <- matrix(c("A", "B", "A", "C"), nrow = 2)
  rownames(aln_matrix) <- c("TR1", "TR2")

  # Test with invalid graph_type
  expect_error(
    plotTR(aln_matrix, motif_map, graph_type = "scatter"),
    "Graph type must be one of (tile, bar)",
    fixed = TRUE
  )

  expect_error(
    plotTR(aln_matrix, motif_map, graph_type = "line"),
    "Graph type must be one of (tile, bar)",
    fixed = TRUE
  )

  expect_error(
    plotTR(aln_matrix, motif_map, graph_type = "histogram"),
    "Graph type must be one of (tile, bar)",
    fixed = TRUE
  )

  # Test with empty string
  expect_error(
    plotTR(aln_matrix, motif_map, graph_type = ""),
    "Graph type must be one of (tile, bar)",
    fixed = TRUE
  )
})

test_that("plotTR validates motif_map is a named vector", {
  aln_matrix <- matrix(c("A", "B", "A", "C"), nrow = 2)
  rownames(aln_matrix) <- c("TR1", "TR2")

  # Test with unnamed vector
  motif_map_unnamed <- c("GCTT", "TCC", "GTG")
  expect_error(
    plotTR(aln_matrix, motif_map_unnamed),
    "motif_map must be named vector",
    fixed = TRUE
  )

  # Test with partially named vector
  motif_map_partial <- c(A = "GCTT", "TCC", C = "GTG")
  expect_error(
    plotTR(aln_matrix, motif_map_partial),
    "Vector is partially named. All elements must have names.",
    fixed = TRUE
  )
})

test_that("plotTR handles valid inputs correctly", {
  motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")
  aln_matrix <- matrix(c("A", "B", "A", "C"), nrow = 2)
  rownames(aln_matrix) <- c("TR1", "TR2")
  colnames(aln_matrix) <- c("Pos1", "Pos2")

  # Should not throw error with valid tile plot
  expect_silent(plotTR(aln_matrix, motif_map, graph_type = "tile"))

  # Should not throw error with valid bar plot
  expect_silent(plotTR(aln_matrix, motif_map, graph_type = "bar"))
})

test_that("plotTR handles matrix without column names", {
  motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")
  aln_matrix <- matrix(c("A", "B", "A", "C"), nrow = 2)
  rownames(aln_matrix) <- c("TR1", "TR2")
  # No column names set

  # Should not throw error and should handle missing column names
  expect_silent(plotTR(aln_matrix, motif_map, graph_type = "tile"))
})

test_that("plotTR handles gap characters correctly", {
  motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")
  aln_matrix <- matrix(c("A", "B", "-", "C"), nrow = 2)
  rownames(aln_matrix) <- c("TR1", "TR2")
  colnames(aln_matrix) <- c("Pos1", "Pos2")

  # Should handle gaps without error
  expect_silent(plotTR(aln_matrix, motif_map, graph_type = "tile"))
  expect_silent(plotTR(aln_matrix, motif_map, graph_type = "bar"))
})

test_that("plotTR returns correct object type", {
  motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")
  aln_matrix <- matrix(c("A", "B", "A", "C"), nrow = 2)
  rownames(aln_matrix) <- c("TR1", "TR2")
  colnames(aln_matrix) <- c("Pos1", "Pos2")

  # Check return type for tile plot
  tile_plot <- plotTR(aln_matrix, motif_map, graph_type = "tile")
  expect_s3_class(tile_plot, "ggplot")

  # Check return type for bar plot
  bar_plot <- plotTR(aln_matrix, motif_map, graph_type = "bar")
  expect_s3_class(bar_plot, "ggplot")
})

test_that("plotTR handles show_motifs parameter correctly", {
  motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")
  aln_matrix <- matrix(c("A", "B", "A", "C"), nrow = 2)
  rownames(aln_matrix) <- c("TR1", "TR2")
  colnames(aln_matrix) <- c("Pos1", "Pos2")

  # Should work with show_motifs = TRUE
  expect_silent(plotTR(aln_matrix, motif_map, show_motifs = TRUE))

  # Should work with show_motifs = FALSE
  expect_silent(plotTR(aln_matrix, motif_map, show_motifs = FALSE))
})

test_that("plotTR handles custom styling parameters", {
  motif_map <- c(A = "GCTT", B = "TCC", C = "GTG")
  aln_matrix <- matrix(c("A", "B", "A", "C"), nrow = 2)
  rownames(aln_matrix) <- c("TR1", "TR2")
  colnames(aln_matrix) <- c("Pos1", "Pos2")

  # Should handle custom parameters without error
  expect_silent(
    plotTR(
      aln_matrix,
      motif_map,
      graph_title = "Custom Title",
      title_size = 14,
      axis_size = 10,
      legend_size = 9
    )
  )
})

#[END]
