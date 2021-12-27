test_that("check_normal_prior works", {
  # This should produce no error messages
  expect_equal(check_normal_prior(c(0, 2), "intercept"), expected=TRUE)

  # Error: non-numeric data
  expect_error(check_normal_prior(c("0", "2"), "intercept"))

  # Error: wrong number of values
  expect_error(check_normal_prior(c(0, 2, 3), "intercept"))
  expect_error(check_normal_prior(c(0), "intercept"))

  # Error: non-positive variance
  expect_error(check_normal_prior(c(0, -2), "intercept"))
  expect_error(check_normal_prior(c(0, 0), "intercept"))
})
