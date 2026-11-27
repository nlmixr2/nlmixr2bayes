# Extracted from test-map-bayes.R:36

# prequel ----------------------------------------------------------------------
.mapBayesMod <- function() {
  ini({
    tcl <- fix(log(0.5))
    tv <- fix(log(5))
    eta.cl + eta.v ~ fix(0.09,
                         0.03, 0.04)
    add.sd <- fix(0.2)
    prop.sd <- fix(0.1)
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(centr) <- -cl / v * centr
    cp <- centr / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}
.mapBayesData <- function() {
  data.frame(
    id = 1,
    time = c(0, 2, 11.9, 12, 24, 26),
    amt = c(100, NA, NA, 100, 100, NA),
    dv = c(NA, 17.5, 6.1, NA, NA, 13.9),
    evid = c(1, 0, 0, 1, 1, 0)
  )
}

# test -------------------------------------------------------------------------
.assert <- function(f) .stanAssert(rxode2::rxode2(f))
