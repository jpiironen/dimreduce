context("spca")
source("helpers.R")

# tests for spca-function


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
  expect_silent(spca(x, y, verbose = F))

  # nctot
  for (nctot in 1:3) {
    dr <- spca(x, y, nctot = nctot, verbose = F)
    z <- predict(dr, x)
    expect_equal(ncol(z), nctot)
  }
})




test_that("result seems plausible", {


  ## data 1
  dat <- dat1
  x <- dat$x
  y <- dat$y

  dr <- spca(x, y, verbose = F)
  z <- predict(dr, x)

  expect_equal(ncol(z), 3)
  expect_equal(dr$v[, 1] != 0, c(T, T, F))
  expect_equal(dr$v[, 2] != 0, c(T, T, F))


  ## data 2
  dat <- dat2
  x <- dat$x
  y <- dat$y

  dr <- spca(x, y, verbose = F)
  z <- predict(dr, x)

  expect_equal(ncol(z), 6)
  expect_equal(dr$v[, 1] != 0, c(T, T, T, T, F, F))
  expect_equal(dr$v[, 2] != 0, c(T, T, T, T, F, F))
  expect_equal(dr$v[, 3] != 0, c(T, T, T, T, F, F))
  expect_equal(dr$v[, 4] != 0, c(T, T, T, T, F, F))
})
