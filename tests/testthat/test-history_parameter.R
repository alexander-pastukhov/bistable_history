test_that("evaluate_history_option works", {
  # type and size
  expect_type(evaluate_history_option("tau", 1, 1, Inf), "list")
  expect_equal(length(evaluate_history_option("tau", 1, 1, Inf)), 5)

  # constant value
  expect_equal(evaluate_history_option("tau", 0.5, 1, Inf)[["fixed_tau"]], 0.5)
  expect_equal(evaluate_history_option("tau", 0.5, 1, Inf)[["tau_option"]], 1)
  expect_equal(evaluate_history_option("tau", 0.5, 1, Inf)[["tau_mu_size"]], 0)
  expect_equal(evaluate_history_option("tau", 0.5, 1, Inf)[["tau_sigma_size"]], 0)
  expect_equal(evaluate_history_option("tau", 0.5, 1, Inf)[["tau_rnd_size"]], 0)

  # single value
  expect_equal(evaluate_history_option("tau", NULL, 1, Inf)[["fixed_tau"]], 0)
  expect_equal(evaluate_history_option("tau", NULL, 1, Inf)[["tau_option"]], 2)
  expect_equal(evaluate_history_option("tau", NULL, 1, Inf)[["tau_mu_size"]], 1)
  expect_equal(evaluate_history_option("tau", NULL, 1, Inf)[["tau_sigma_size"]], 0)
  expect_equal(evaluate_history_option("tau", NULL, 1, Inf)[["tau_rnd_size"]], 0)

  # single value despite declaring a random factor
  expect_equal(evaluate_history_option("tau", "random", 1, Inf)[["tau_option"]], 2)
  expect_equal(evaluate_history_option("tau", "1|random", 1, Inf)[["tau_option"]], 2)

  # independent random
  randomN <- 10
  expect_equal(evaluate_history_option("tau", "random", randomN, Inf)[["fixed_tau"]], 0)
  expect_equal(evaluate_history_option("tau", "random", randomN, Inf)[["tau_option"]], 3)
  expect_equal(evaluate_history_option("tau", "random", randomN, Inf)[["tau_mu_size"]], 0)
  expect_equal(evaluate_history_option("tau", "random", randomN, Inf)[["tau_sigma_size"]], 0)
  expect_equal(evaluate_history_option("tau", "random", randomN, Inf)[["tau_rnd_size"]], randomN)

  # pooled random
  randomN <- 10
  expect_equal(evaluate_history_option("tau", "1|random", randomN, Inf)[["fixed_tau"]], 0)
  expect_equal(evaluate_history_option("tau", "1|random", randomN, Inf)[["tau_option"]], 4)
  expect_equal(evaluate_history_option("tau", "1|random", randomN, Inf)[["tau_mu_size"]], 1)
  expect_equal(evaluate_history_option("tau", "1|random", randomN, Inf)[["tau_sigma_size"]], 1)
  expect_equal(evaluate_history_option("tau", "1|random", randomN, Inf)[["tau_rnd_size"]], randomN)
})

test_that("evaluate_history_option fails", {
  # unknown option
  expect_error(evaluate_history_option("Random", -0.5, 1, Inf))

  # supplied value outside of valid range
  expect_error(evaluate_history_option("tau", -0.5, 1, Inf))
  expect_error(evaluate_history_option("history_mix", -0.5, 1, 1))
  expect_error(evaluate_history_option("history_mix", 1.1, 1, 1))
})


test_that("evaluate_history_init works", {
  # returns numeric vector with two elements
  expect_vector(evaluate_history_init(0.5), ptype=double(), size=2)
  expect_vector(evaluate_history_init(c(0, 0)), ptype=double(), size=2)
})

test_that("evaluate_history_init fails", {
  # too few or too many elements
  expect_error(evaluate_history_init(c()))
  expect_error(evaluate_history_init(c(0, 0, 0)))

  # wrong type
  expect_error(evaluate_history_init(c("0", "0")))
  expect_error(evaluate_history_init(c(TRUE, FALSE)))

  # wrong range
  expect_error(evaluate_history_init(1.1))
  expect_error(evaluate_history_init(-0.1))
  expect_error(evaluate_history_init(c(-0.1, 0.5)))
  expect_error(evaluate_history_init(c(-0.1, 1.1)))
})
