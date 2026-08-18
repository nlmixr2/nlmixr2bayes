## linCmt() support.  rxode2 will not give sensitivities for the closed form --
## calcSens creates the rx__sens_* compartments but nothing writes to them, so
## they solve to all zeros -- and its real linCmt Jacobians live in a separate
## internal path that rxSolve does not expose.  The generator therefore
## replaces linCmt() with the equivalent ODE system and checks the replacement
## against linCmt() itself.
library(testthat)

oralData <- function(nsub = 3L, cmt = "depot") {
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- rxode2::et(amt = 100, cmt = cmt)
    e <- rxode2::et(e, c(0.25, 1, 2, 4, 8, 12, 24))
    d <- as.data.frame(e)
    d$id <- i
    d$dv <- 1
    d
  }))
}

lin1Oral <- function() {
  ini({
    tcl <- 1.386
    tv <- 3.401
    tka <- 0.0953
    eta.cl ~ 0.0625
    add.sd <- 0.3
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    ka <- exp(tka)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

lin2Oral <- function() {
  ini({
    tcl <- 1.386
    tv <- 3.401
    tq <- 0.7
    tvp <- 3.0
    tka <- 0.0953
    add.sd <- 0.3
  })
  model({
    cl <- exp(tcl)
    v <- exp(tv)
    q <- exp(tq)
    vp <- exp(tvp)
    ka <- exp(tka)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

lin1Iv <- function() {
  ini({
    tcl <- 1.386
    tv <- 3.401
    add.sd <- 0.3
  })
  model({
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

test_that("rxode2 gives NO usable sensitivities for linCmt", {
  ## Documents why the expansion exists.  If this ever starts failing, rxode2
  ## has grown real linCmt sensitivities and the expansion can be revisited.
  m <- rxode2::rxode2("cl <- Cl\nv <- V\nka <- Ka\ncp <- linCmt()",
                      calcSens = c("Cl", "V", "Ka"))
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, c(0.5, 2, 6))
  s <- rxode2::rxSolve(m, params = c(Cl = 4, V = 30, Ka = 1.1), events = ev,
                       returnType = "data.frame", cores = 1L,
                       atol = 1e-12, rtol = 1e-12)

  sens <- grep("^rx__sens_", names(s), value = TRUE)
  expect_gt(length(sens), 0)
  expect_true(all(vapply(sens, function(k) all(s[[k]] == 0), logical(1))))

  ## And the columns are not even named the way ODE models name them.
  expect_false(any(grepl("__$", sens)))
})

test_that("the expansion reproduces linCmt() itself", {
  skip_if_not_installed("nlmixr2")

  for (mod in list(lin1Oral, lin2Oral)) {
    ui <- rxode2::rxode2(mod)
    info <- nlmixr2bayes:::.rxsLinCmtInfo(ui)
    txt <- paste(c(grep("linCmt\\(", strsplit(rxode2::rxNorm(ui), "\n")[[1]],
                        invert = TRUE, value = TRUE),
                   nlmixr2bayes:::.rxsLinCmtOdes(info),
                   sprintf("cp = central / %s;", info$v)), collapse = "\n")
    ## Errors if the expansion and the closed form disagree.
    d <- nlmixr2bayes:::.rxsCheckLinCmt(ui, txt, info)
    expect_lt(d, 1e-8)
  }
})

test_that("compartment order matches linCmt's, so cmt numbers still line up", {
  skip_if_not_installed("nlmixr2")
  ui <- rxode2::rxode2(lin2Oral)
  info <- nlmixr2bayes:::.rxsLinCmtInfo(ui)
  expect_equal(info$states, c("depot", "central", "peripheral1"))

  txt <- paste(c("cl <- 4", "v <- 30", "q <- 2", "vp <- 20", "ka <- 1.1",
                 nlmixr2bayes:::.rxsLinCmtOdes(info)), collapse = "\n")
  expect_equal(rxode2::rxState(rxode2::rxode2(txt)),
               c("depot", "central", "peripheral1"))
})

test_that("a linCmt model generates and its states come from the expansion", {
  skip_if_not_installed("nlmixr2")

  gen <- rxsStanFromUi(lin1Oral, oralData(2L))
  on.exit(rxsRelease(gen$handle))

  expect_equal(gen$states, c("depot", "central"))
  expect_match(gen$code, "real cp = (central / v);", fixed = TRUE)
  ## 3 thetas + 1 eta per subject
  expect_equal(attr(gen$handle, "nBlock"), 4L)
})

test_that("an IV linCmt model works too", {
  skip_if_not_installed("nlmixr2")
  gen <- rxsStanFromUi(lin1Iv, oralData(2L, cmt = "central"))
  on.exit(rxsRelease(gen$handle))
  expect_equal(gen$states, "central")
})

test_that("a parameterization the expansion cannot cover is refused", {
  skip_if_not_installed("nlmixr2")
  ## Rate-constant parameterization: no cl/v for the expansion to use.
  micro <- function() {
    ini({
      tkel <- -2
      tv <- 3.401
      tka <- 0.0953
      add.sd <- 0.3
    })
    model({
      kel <- exp(tkel)
      v <- exp(tv)
      ka <- exp(tka)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  expect_error(rxsStanFromUi(micro, oralData(2L)),
               "clearance parameterization|disagrees with linCmt")
})

test_that("a linCmt model has correct gradients end to end", {
  skipUnlessStan()

  gen <- rxsStanFromUi(lin1Oral, oralData(4L))
  on.exit(rxsRelease(gen$handle))

  sm <- stanModelFor(gen$code, "rxstan_lincmt")
  fit <- rstan::sampling(sm, data = gen$standata, chains = 0)

  set.seed(31)
  for (trial in 1:3) {
    u <- stats::runif(rstan::get_num_upars(fit), -0.5, 0.5)
    chk <- rxsCheckGradient(fit, u)
    expect_true(all(chk$relDiff < 1e-4),
                info = paste(utils::capture.output(
                  print(chk[order(-chk$relDiff), ][1:3, ])), collapse = "\n"))
  }
})
