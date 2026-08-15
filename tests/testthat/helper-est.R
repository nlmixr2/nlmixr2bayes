# A prior-carrying analytic model for the est="stan" tests
.estMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- c(0, 0.5)
    eta.cl ~ 0.1
    prior(tcl) ~ dnorm(1, 2)
    prior(tv) ~ dnorm(3, 2)
    prior(add.sd) ~ dcauchy(0, 2.5)
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}
