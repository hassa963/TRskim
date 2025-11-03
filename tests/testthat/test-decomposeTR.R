test_that("decomposeTR correctly decomposes allele into motif composition", {
  allele <- DNAString("TTTACACGTAC")
  motifs <- DNAStringSet(c("AC", "GT"))

  result <- decomposeTR(allele, motifs)

  expected <- c("TTT", "AC", "AC", "GT", "AC")

  expect_equal(result$composition, expected)
})

test_that("decomposeTR correctly decomposes mismatched motif sizes", {
  allele <-DNAString("TTTACACGTCAC")
  motifs <- DNAStringSet(c("AC", "GTC"))

  result <- decomposeTR(allele, motifs)

  expected <- c("TTT", "AC", "AC", "GTC", "AC")

  expect_equal(result$composition, expected)
})

test_that("decomposeTR creates empty character", {
  allele <-DNAString("TTTACACGTCAC")
  motifs <- DNAStringSet(c("GC"))

  result <- decomposeTR(allele, motifs)

  expected <- character(0)

  expect_equal(result, expected)
})
