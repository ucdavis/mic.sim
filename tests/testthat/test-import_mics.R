test_that("multiplication works", {
  expected = tibble(
      left_bound = c(1, 3, -Inf, -Inf, NA),
      right_bound = c(2, Inf, 1, 1, NA),
      mic_column = c(4, ">8", "≤2", "=<2", NA))
  attr(expected, "source") = "imported"
  attr(expected, "lr_col") = FALSE
  attr(expected, "mic_class") = "imported_mic_column"
  attr(expected, "metadata") = FALSE
  attr(expected, "scale") = "log"

  actual = import_mics(mic_column = c(4, ">8", "≤2", "=<2", NA))
  expect_equal(actual, expected)
})
