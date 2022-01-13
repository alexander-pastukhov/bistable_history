test_that("test fit and post fit utilities", {
  fit <- readRDS("fit.RDS")
  sampler_params <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)
  samplesN <- dim(sampler_params[[1]])[1] * length(sampler_params)

  # check that the fit itself is what we think it is
  expect_s3_class(fit, "cumhist")

  # checking that extraction methods work
  expect_s3_class(bayes_R2(fit), "data.frame")

  # history term
  expect_s3_class(bistablehistory::historyef(fit), "data.frame")
  expect_equal(nrow(bistablehistory::historyef(fit, summary = FALSE)), 2 * samplesN)

  # tau
  expect_s3_class(bistablehistory::history_tau(fit), "data.frame")
  expect_equal(nrow(bistablehistory::history_tau(fit, summary = FALSE)), samplesN)
  expect_equal(ncol(bistablehistory::history_tau(fit, summary = FALSE)), 1)

  # mixed state
  expect_s3_class(bistablehistory::history_mixed_state(fit), "data.frame")

  # mean history
  expect_equal(nrow(extract_history(fit)), fit$data$rowsN)
  expect_equal(ncol(extract_history(fit)), 2)

  # predictions
  expect_equal(length(predict(fit)), fit$data$rowsN)
  expect_equal(length(predict(fit, full_length = FALSE)), fit$data$clearN)
  expect_equal(nrow(predict(fit, summary = FALSE)), samplesN)
  expect_equal(ncol(predict(fit, summary = FALSE)), fit$data$rowsN)
  expect_equal(ncol(predict(fit, summary = FALSE, full_length = FALSE)), fit$data$clearN)
  expect_equal(nrow(predict(fit, summary = FALSE)), samplesN)

  # not fixed effect, should be null
  expect_message(fixef(fit))
  expect_equal(fixef(fit), NULL)

  # fixed effects
  fit <- readRDS("fit-fixed.RDS")
  sampler_params <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)
  samplesN <- dim(sampler_params[[1]])[1] * length(sampler_params)
  expect_s3_class(fixef(fit), "data.frame")
  expect_equal(nrow(fixef(fit)), 2)
  expect_equal(nrow(fixef(fit, summary=FALSE)), 2 * samplesN)
  expect_equal(ncol(fixef(fit, summary=FALSE)), 3)
})
