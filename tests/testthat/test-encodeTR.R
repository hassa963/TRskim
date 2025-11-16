
test_that("encodeTR correctly encodes decomposed motifs", {
  decomposed_TR <- c("AC", "AC", "GT", "AC")
  motifs <- c("AC", "GT")

  encoded <- encodeTRs(decomposed_TR, motifs)

  # Expected outputs
  expected_encoded <- "AABA"
  expected_map <- c("A", "B")
  names(expected_map) <- c("AC", "GT")

  expect_named(encoded, c("encoded", "map"))
  expect_equal(encoded$encoded, expected_encoded)
  expect_equal(encoded$map, expected_map)
})
