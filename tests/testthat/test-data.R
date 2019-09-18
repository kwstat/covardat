# test-data.R

test_that("data accessible", {
  expect_silent(data(tumor))
  expect_silent(data(alzheimers))
})
