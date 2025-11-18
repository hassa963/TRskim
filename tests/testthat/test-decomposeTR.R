#OpenAI. ChatGPT (GPT-5) large language model (2025).
# https://chat.openai.com/

#supported by AI to write tests

test_that("decomposeTRs correctly decomposes allele into motif composition", {
  allele <- DNAString("TTTACACGTAC")
  motifs <- DNAStringSet(c("AC", "GT"))

  result <- decomposeTRs(allele, motifs)

  expected <- c("TTT", "AC", "AC", "GT", "AC")

  expect_equal(result$compositions, expected)
})

test_that("decomposeTRs correctly decomposes mismatched motif sizes", {
  allele <-DNAString("TTTACACGTCAC")
  motifs <- DNAStringSet(c("AC", "GTC"))

  result <- decomposeTRs(allele, motifs)

  expected <- c("TTT", "AC", "AC", "GTC", "AC")

  expect_equal(result$compositions, expected)
})

test_that("decomposeTRs creates empty character", {
  allele <- DNAString("TTTACACGTCAC")
  motifs <- DNAStringSet(c("GC"))

  expect_warning(
    result <- decomposeTRs(allele, motifs),
    "No motif matches found in allele"
  )
  expected <- character(0)
  expect_equal(result$compositions, expected)

  allele <- DNAString("G")
  motifs <- DNAStringSet(c("GC"))

  expect_warning(
    result <- decomposeTRs(allele, motifs),
    "No motif matches found in allele"
  )
  expected <- character(0)
  expect_equal(result$compositions, expected)
})

test_that("decomposeTRs prioritizes longer motifs", {
  allele <-DNAString("ATATATAT")
  motifs <- DNAStringSet(c("AT", "ATAT"))

  result <- decomposeTRs(allele, motifs)

  expected <- c("ATAT", "ATAT")

  expect_equal(result$compositions, expected)
})

test_that("decomposeTRs takes in character vectors", {
  allele <-"ATATATAT"
  motifs <- c("AT", "ATAT")

  result <- decomposeTRs(allele, motifs)

  expected <- c("ATAT", "ATAT")

  expect_equal(result$compositions, expected)
})

test_that("decomposeTRs rejects invalid nucleotide characters and empty inputs",
          {

  expect_error(decomposeTRs("AC1T", c("AC")),
               "Alleles contain invalid nucleotide characters")

  expect_error(decomposeTRs(NULL, c("AC")),
               "Alleles cannot be NULL or empty")

  expect_error(decomposeTRs("ACT", NULL),
               "Motifs cannot be NULL or empty")

  expect_error(decomposeTRs("ACT", c()),
               "Motifs cannot be NULL or empty")

  expect_error(decomposeTRs(c(), c("AC")),
               "Alleles cannot be NULL or empty")

  expect_error(decomposeTRs(c("ACGT", "ACGT"), c("4C")),
               "Motifs contain invalid nucleotide characters")

  expect_error(decomposeTRs(c("ACGT", "A1GT"), c("AC")),
               "Alleles contain invalid nucleotide characters")

  expect_error(decomposeTRs(c("ACGT", ""), c("AC")),
               "Alleles contain invalid nucleotide characters")

  expect_error(decomposeTRs(c("ACGT", ""), c("")),
               "Alleles contain invalid nucleotide characters")
})
