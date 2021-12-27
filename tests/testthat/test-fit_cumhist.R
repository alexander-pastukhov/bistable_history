test_that("test the fit_cumhist fails", {
  # should work without sampling, returning cumhist s3 class
  expect_s3_class(fit_cumhist(br_singleblock, state="State", duration="Duration", chains = 0), "cumhist")

  # should fail for data parameter
  expect_error(fit_cumhist(data=NULL, state="State", duration="Duration", chains = 0))
  expect_error(fit_cumhist(list(br_singleblock), state="State", duration="Duration", chains = 0))

  # should fail as state was not specified
  expect_error(fit_cumhist(br_singleblock, duration="Duration", chains = 0))

  # should fail as neither duration nor onset were specified
  expect_error(fit_cumhist(br_singleblock, state="State", chains = 0))

  # should fail because of bad tau parameter
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", tau = -1, chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", tau = c(0.5, 0.5), chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", tau = "RANDOM", chains = 0))

  # should fail because of bad history_mix parameter
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_mix = -1, chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_mix = 2, chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_mix = c(0.5, 0.5), chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_mix = "RANDOM", chains = 0))

  # should fail because of bad mixed_state parameter
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", mixed_state = -1, chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", mixed_state = 2, chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", mixed_state = c(0.5, 0.5), chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", mixed_state = "RANDOM", chains = 0))

  # should fail because of bad history_init values
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_init = -1, chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_init = 1.1, chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_init = c(0.2, 0.2, 0.2), chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_init = "one-half", chains = 0))

  # should fail because of bad prior
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_priors = c(1, 1), chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_priors = list(c(1, 1)), chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", history_priors = list("tau" = c(1, 1, 2)), chains = 0))

  # should fail because of bad intercept prior
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", intercept_priors = list(c(1, 1)), chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", intercept_priors = c(1, 1), chains = 0))
  expect_error(fit_cumhist(br_singleblock, state="State", duration="Duration", intercept_priors = c(1, 1, 1, 1), family = "normal", chains = 0))
})
