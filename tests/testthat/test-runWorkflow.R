library(testthat)

test_that("runWorkflow produces expected messages", {
  # Set up test data
  alleles <- c("AGATAGATAGATAGAT", "AGATAGATAGAT", "AGATAGAT")
  motifs <- c("AGAT")

  # Capture messages using expect_message
  expect_message(
    result <- runWorkflow(alleles, motifs),
    "Decomposing tandem repeats..."
  )

  expect_message(
    result <- runWorkflow(alleles, motifs),
    "Encoding tandem repeats...."
  )

  expect_message(
    result <- runWorkflow(alleles, motifs),
    "Aligning tandem repeats...."
  )

  expect_message(
    result <- runWorkflow(alleles, motifs),
    "Generating Tile plot visualization..."
  )

  expect_message(
    result <- runWorkflow(alleles, motifs),
    "Generating Bar plot visualization..."
  )

  expect_message(
    result <- runWorkflow(alleles, motifs),
    "Workflow complete!"
  )
})

test_that("runWorkflow produces all messages in order", {
  alleles <- c("AGATAGATAGATAGAT", "AGATAGATAGAT", "AGATAGAT")
  motifs <- c("AGAT")

  # Capture all messages at once
  messages <- character()

  result <- withCallingHandlers(
    runWorkflow(alleles, motifs),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
    }
  )

  # Check all expected messages are present in order
  expect_length(messages, 6)
  expect_match(messages[1], "Decomposing tandem repeats...")
  expect_match(messages[2], "Encoding tandem repeats....")
  expect_match(messages[3], "Aligning tandem repeats....")
  expect_match(messages[4], "Generating Tile plot visualization...")
  expect_match(messages[5], "Generating Bar plot visualization...")
  expect_match(messages[6], "Workflow complete!")
})

test_that("runWorkflow returns expected list structure", {
  alleles <- c("AGATAGATAGATAGAT", "AGATAGATAGAT", "AGATAGAT")
  motifs <- c("AGAT")

  # Suppress messages for cleaner test output
  result <- suppressMessages(runWorkflow(alleles, motifs))

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("decomposition", "encoded", "alignment",
                         "motif_map", "tile_plot", "bar_plot"))

  # Check each component exists and is correct type
  expect_true(!is.null(result$decomposition))
  expect_true(!is.null(result$encoded))
  expect_true(!is.null(result$alignment))
  expect_true(!is.null(result$motif_map))
  expect_true(!is.null(result$tile_plot))
  expect_true(!is.null(result$bar_plot))

  # Check specific types
  expect_true(is.matrix(result$alignment))
  expect_true(is.character(result$motif_map))
  expect_true(!is.null(names(result$motif_map)))  # motif_map must be named
})

test_that("runWorkflow passes parameters correctly", {
  alleles <- c("AGATAGATAGATAGAT", "AGATAGATAGAT", "AGATAGAT")
  motifs <- c("AGAT")

  # Test with custom parameters
  result <- suppressMessages(
    runWorkflow(alleles, motifs,
                 match_score = 2,
                 indel = -2,
                 allowance = 1)
  )

  # Should complete without error
  expect_true(!is.null(result))
  expect_type(result, "list")
})

test_that("runWorkflow handles edge cases", {
  # Single allele
  alleles <- c("AGATAGAT")
  motifs <- c("AGAT")

  result <- suppressMessages(runWorkflow(alleles, motifs))
  expect_true(!is.null(result))
  expect_equal(nrow(result$alignment), 1)

  # Multiple motifs
  alleles2 <- c("AGATAGACAGAT", "AGATAGAT", "AGACAGAC")
  motifs2 <- c("AGAT", "AGAC")

  result2 <- suppressMessages(runWorkflow(alleles2, motifs2))
  expect_true(!is.null(result2))
  expect_true(length(result2$motif_map) >= 2)
})

test_that("runWorkflow produces valid plots", {
  alleles <- c("AGATAGATAGATAGAT", "AGATAGATAGAT", "AGATAGAT")
  motifs <- c("AGAT")

  result <- suppressMessages(runWorkflow(alleles, motifs))

  # Check plots are ggplot objects
  expect_s3_class(result$tile_plot, "ggplot")
  expect_s3_class(result$bar_plot, "ggplot")

  # Verify plots are different (one is tile, one is bar)
  expect_false(identical(result$tile_plot, result$bar_plot))
})

# Integration test - run the full workflow and verify end-to-end
test_that("runWorkflow integration test", {
  # Realistic STR data
  alleles <- c(
    "AGATAGATAGATAGATAGATAGAT",  # 6 repeats
    "AGATAGATAGATAGATAGAT",        # 5 repeats
    "AGATAGATAGATAGAT"              # 4 repeats
  )
  motifs <- c("AGAT")

  # Run full workflow
  messages <- character()
  result <- withCallingHandlers(
    runWorkflow(alleles, motifs),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
    }
  )

  # Verify complete execution
  expect_length(messages, 6)
  expect_equal(nrow(result$alignment), 3)
  expect_type(result$encoded, "list")
  expect_length(result$encoded, 3)

  # Verify motif_map is properly named
  expect_true(!is.null(names(result$motif_map)))
  expect_true(all(nchar(names(result$motif_map)) > 0))

  # Verify plots can be rendered (don't error)
  expect_silent(print(result$tile_plot))
  expect_silent(print(result$bar_plot))
})

#[END]
