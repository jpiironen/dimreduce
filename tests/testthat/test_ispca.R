context("ispca")
source("helpers.R")

# tests for ispca-function


# get some data
set.seed(130934)
n <- 200
dat1 <- get.test.data(n, nfeat_per_latent = 1)
dat2 <- get.test.data(n, nfeat_per_latent = 2)




test_that("arguments to ispca have expected effect", {
  dat <- dat1
  x <- dat$x
  y <- dat$y

  # verbosity
  expect_silent(ispca(x, y, verbose = F))

  # nctot
  for (nctot in 1:3) {
    dr <- ispca(x, y, nctot = nctot, verbose = F)
    z <- predict(dr, x)
    expect_equal(ncol(z), nctot)
  }
  expect_error(ispca(x, y, nctot = ncol(x) + 1))

  # ncsup
  for (ncsup in 1:3) {
    dr <- ispca(x, y, ncsup = ncsup, verbose = F, permtest = F)
    z <- predict(dr, x)
    expect_equal(ncol(z), ncsup)
  }
})





test_that("result seems plausible", {


  ## data 1
  dat <- dat1
  x <- dat$x
  y <- dat$y

  dr <- ispca(x, y, verbose = F)
  z <- predict(dr, x)

  expect_equal(ncol(z), 2)
  expect_equal(dr$v[, 1], c(1, 0, 0))
  expect_equal(dr$v[, 2], c(0, 1, 0))


  ## data 2
  dat <- dat2
  x <- dat$x
  y <- dat$y

  dr <- ispca(x, y, verbose = F)
  z <- predict(dr, x)

  expect_equal(ncol(z), 2)
  expect_equal(dr$v[, 1] != 0, c(T, T, F, F, F, F))
  expect_equal(dr$v[, 2] != 0, c(F, F, T, T, F, F))
})
