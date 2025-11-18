#OpenAI. ChatGPT (GPT-5) large language model (2025).
# https://chat.openai.com/

#supported by AI to write tests

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


test_that("alignTRs accepts list of encoded strings", {
  encoded_trs <- list("AABA", "ABA", "BA")
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

test_that("alignTRs handles single element list", {
  encoded_trs <- list("AABA")
  result <- alignTRs(encoded_trs)

  expected <- matrix(
    c("A", "A", "B", "A"),
    nrow = 1, byrow = TRUE,
    dimnames = list("TR1")
  )

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(1, 4))
  expect_equal(result, expected)
})

test_that("alignTRs preserves names from named list", {
  encoded_trs <- list(sample1 = "AAB", sample2 = "AB", sample3 = "B")
  result <- alignTRs(encoded_trs)

  expect_equal(rownames(result), c("sample1", "sample2", "sample3"))
})

test_that("alignTRs handles nested list structure", {
  # Test if someone passes the full output from encodeTRs
  encoded_output <- list(encoded = list("AABA", "ABA", "BA"))

  # This should work by passing just the encoded part
  result <- alignTRs(encoded_output$encoded)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 4))
})

test_that("alignTRs errors with list containing non-character elements", {
  encoded_trs <- list("AAB", 123, "BA")

  expect_error(
    alignTRs(encoded_trs),
    "encoded_trs must be a character vector or list of character strings"
  )
})

test_that("alignTRs errors with empty list", {
  expect_error(
    alignTRs(list()),
    "encoded_trs cannot be empty"
  )
})

test_that("alignTRs errors with list containing NA", {
  encoded_trs <- list("AAB", NA, "BA")

  expect_error(
    alignTRs(encoded_trs),
    "encoded_trs must be a character vector or list of character strings"
  )
})

test_that("alignTRs errors with list containing empty strings", {
  encoded_trs <- list("AAB", "", "BA")

  expect_error(
    alignTRs(encoded_trs),
    "encoded_trs cannot contain empty strings"
  )
})
