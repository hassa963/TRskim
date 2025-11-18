test_that("encodeTR correctly encodes decomposed motifs", {
  decomposed_TR <- c("AC", "AC", "GT", "AC")
  motifs <- c("AC", "GT")
  encoded <- encodeTRs(decomposed_TR, motifs)
  # Expected outputs
  expected_encoded <- "AABA"
  expected_map <- c("A", "B")
  names(expected_map) <- c("AC", "GT")
  expect_named(encoded, c("encoded", "motif_map"))
  expect_equal(encoded$encoded, expected_encoded)
  expect_equal(encoded$motif_map, expected_map)
})

test_that("encodeTRs errors when decomposed_TRs is NULL", {
  expect_error(
    encodeTRs(NULL, c("AC", "GT")),
    "decomposed_TRs cannot be NULL"
  )
})

test_that("encodeTRs errors when motifs is NULL", {
  expect_error(
    encodeTRs(c("AC"), NULL),
    "motifs cannot be NULL"
  )
})

test_that("encodeTRs errors when motifs is empty", {
  expect_error(
    encodeTRs(c("AC"), character(0)),
    "motifs cannot be NULL or empty"  # Updated error message
  )
})

test_that("encodeTRs errors when motifs contain invalid nucleotides", {
  motifs <- c("AC", "1G")  # 1 is invalid
  expect_error(
    encodeTRs(c("AC", "GT"), motifs),
    "motifs must be of type character with proper nucleotide characters"
  )
})

test_that("encodeTRs errors when decomposed_TRs contains invalid characters", {
  decomposed <- c("AC", "GT", "22")
  expect_error(
    encodeTRs(decomposed, c("AC", "GT")),
    "decomposed_TRs must contain proper nucleotide characters"
    # Updated error message
  )
})

test_that("encodeTRs correctly encodes multiple TR entries in a list", {
  decomposed_list <- list(
    c("AC", "GT"),
    c("GT", "AC")
  )
  motifs <- c("AC", "GT")
  out <- encodeTRs(decomposed_list, motifs)
  expect_equal(out$encoded[[1]], "AB")
  expect_equal(out$encoded[[2]], "BA")
  expect_equal(out$motif_map, c(AC="A", GT="B"))
})

test_that("encodeTRs warns and removes empty TR entries", {
  decomposed_list <- list(
    c("AC","GT"),
    character(0),
    c("", "")
  )
  motifs <- c("AC", "GT")
  expect_warning(
    out <- encodeTRs(decomposed_list, motifs),
    "Removed"
  )
  expect_length(out$encoded, 1)  # Changed from expect_equal
  expect_equal(out$encoded[[1]], "AB")  # List element, not bare list
})

test_that("encodeTRs produces correct warning for removed TRs", {
  decomposed <- list(
    c("AC","GT"),
    character(0),
    c("", "")
  )
  motifs <- c("AC","GT")
  expect_warning(
    encodeTRs(decomposed, motifs),
    "Removed 2 invalid or empty TR entries at positions: 2, 3"
  )
})

test_that("encodeTRs returns list when input is a list", {
  decomposed <- list(c("AC", "GT"))
  motifs <- c("AC", "GT")
  out <- encodeTRs(decomposed, motifs)
  expect_type(out$encoded, "list")
  expect_length(out$encoded, 1)
  expect_equal(out$encoded[[1]], "AB")
})

test_that("encodeTRs errors when all TR entries are invalid or empty", {
  decomposed_list <- list(
    character(0),
    c("", "")
  )
  # Should get warning first, then error
  expect_error(
    suppressWarnings(encodeTRs(decomposed_list, c("AC", "GT"))),
    "All decomposed_TRs entries were empty or invalid"
  )
})

test_that("encodeTRs errors when motif count exceeds available symbols", {
  # Create valid nucleotide motifs that exceed symbol limit
  motifs <- paste0(rep(c("A", "C", "G", "T"), length.out = 100),
                   rep(c("A", "C", "G", "T"), each = 25))
  expect_error(
    encodeTRs(c("AA"), motifs),
    "Too many motifs to encode"
  )
})

test_that("encodeTRs accepts DNAStringSet for motifs", {
  decomposed_TR <- c("AC", "AC", "GT", "AC")
  motifs <- Biostrings::DNAStringSet(c("AC", "GT"))

  encoded <- encodeTRs(decomposed_TR, motifs)

  expected_encoded <- "AABA"
  expected_map <- c("A", "B")
  names(expected_map) <- c("AC", "GT")

  expect_equal(encoded$encoded, expected_encoded)
  expect_equal(encoded$motif_map, expected_map)
})

#[END]
