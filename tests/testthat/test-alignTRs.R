test_that("alignTRs correctly aligns encoded TR sequences", {
  encoded_trs <- c("AABA", "ABA", "BA")

  result <- alignTRs(encoded_trs)

  expected <- matrix(
    c("A", "A", "B", "A",
      "-", "A", "B", "A",
      "-", "-", "B", "A"),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      c("TR1", "TR2", "TR3")
    )
  )

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 4))
  expect_equal(result, expected)
})
