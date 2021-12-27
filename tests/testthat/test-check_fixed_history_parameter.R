test_that("check_fixed_history_parameter works", {
  # This should produce no error messages and return rep(1, 10)
  expect_equal(check_fixed_history_parameter("tau", 1, 10, Inf), rep(1, 10))
  expect_equal(check_fixed_history_parameter("tau", 1, 1, Inf), rep(1, 1))

  # Wrong type for value
  expect_error(check_fixed_history_parameter("tau", "1", 10, Inf))
  expect_error(check_fixed_history_parameter("tau", list(1), 10, Inf))
  expect_error(check_fixed_history_parameter("tau", as.integer(c()), 10, Inf))

  # Positive integer scalar for randomN
  expect_error(check_fixed_history_parameter("tau", 1, 10.1, Inf))
  expect_error(check_fixed_history_parameter("tau", 1, c(10, 10), Inf))
  expect_error(check_fixed_history_parameter("tau", 1, 0, Inf))

  # Wrong number of values, neither 1, nor randomN
  expect_error(check_fixed_history_parameter("tau", c(1, 2), 10, Inf))
  expect_error(check_fixed_history_parameter("tau", as.numeric(c()), 10, Inf))

  # Value above upper limit
  expect_error(check_fixed_history_parameter("tau", 2, 10, 1))
})
