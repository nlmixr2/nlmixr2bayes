# Shared fixture: a small analytic (non-ODE) model so the link tests compile
# fast and the conditional density is hand-computable; same shape as
# nlmixr2est's own foceiLik fixtures.

.linkMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- 0.5
    eta.cl ~ 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

.linkData <- function() {
  set.seed(42)
  do.call(rbind, lapply(1:4, function(id) {
    tt <- c(0.5, 1, 2, 4, 8)
    data.frame(ID = id, TIME = tt,
               DV = 5 * exp(-0.05 * tt) + stats::rnorm(length(tt), 0, 0.5),
               AMT = 0, EVID = 0)
  }))
}
