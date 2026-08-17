# Prior emission (issue #1, Spec 4): turn rxode2::rxUiPriors() rows into Stan
# declarations and sampling statements via lotri::lotriPriorDists()$stanName.
# Pure R -- no Stan, no C.  Four load-bearing rules:
#
# 1. SUPPORT PROMOTION.  lotri validates a prior's support only against
#    *finite* parameter bounds, so dlnorm() on an unbounded parameter passes
#    validation and would yield -Inf at init; the declaration is tightened by
#    the prior's support first.
# 2. TRUNCATION.  With effective bounds (L, U): both infinite -> plain
#    statement; finite with every argument a numeric literal -> declare the
#    constraint and OMIT T[L,U] (the dropped normalizer is a constant -- the
#    half-Cauchy idiom); finite with any argument referencing a parameter ->
#    emit T[L,U] (the normalizer is parameter-dependent; dropping it is a real
#    bug), erroring if Stan has no _lcdf for the distribution.
# 3. MULTIVARIATE DEDUP.  lotri stores a multivariate prior byte-identically
#    on every member row; dedupe by (prior string, sorted member set), with
#    member ORDER from the embedded lotri() expression's dimnames.
# 4. DISCRETE REFUSAL.  Every iniDf parameter is real-valued; the 8 discrete
#    distributions are refused with an R error naming the parameter, rather
#    than a stanc type error later.

# Stan distributions in the lotri catalogue without a _lcdf/_lccdf (cannot be
# truncated with parameter-dependent bounds/arguments)
.stanNoLcdf <- c("wiener")

#' Mangle an nlmixr2 parameter name into a Stan identifier
#' @noRd
.stanParName <- function(x) {
  .x <- gsub(".", "_", x, fixed = TRUE)
  .x <- gsub("[^A-Za-z0-9_]", "_", .x)
  ifelse(grepl("^[0-9_]", .x), paste0("p_", .x), .x)
}

#' Format a number for Stan source
#' @noRd
.stanNum <- function(x) {
  vapply(x, function(v) {
    if (is.infinite(v)) stop("cannot emit an infinite value into Stan source",
                             call. = FALSE) # nocov
    format(v, digits = 15, trim = TRUE, scientific = FALSE)
  }, character(1))
}

#' The bounds declaration for a Stan parameter, "" when unbounded
#' @noRd
.stanConstraint <- function(lower, upper) {
  .l <- is.finite(lower)
  .u <- is.finite(upper)
  if (.l && .u) return(paste0("<lower=", .stanNum(lower), ",upper=", .stanNum(upper), ">"))
  if (.l) return(paste0("<lower=", .stanNum(lower), ">"))
  if (.u) return(paste0("<upper=", .stanNum(upper), ">"))
  ""
}

#' Is a prior argument a numeric literal (or an expression of literals)?
#'
#' `dcauchy(0, 5)` -> both TRUE; `dcauchy(0, tauScale)` -> the second is
#' FALSE.  Evaluated against baseenv() only, so a parameter name can never
#' accidentally resolve.
#' @noRd
.stanArgIsLiteral <- function(arg) {
  if (is.numeric(arg)) return(TRUE)
  .v <- tryCatch(eval(arg, envir = baseenv()), error = function(e) NULL)
  is.numeric(.v) && length(.v) == 1L
}

#' Deparse a prior argument into Stan source, mangling parameter references
#' @noRd
.stanArgDeparse <- function(arg) {
  if (.stanArgIsLiteral(arg)) {
    return(.stanNum(eval(arg, envir = baseenv())))
  }
  if (is.name(arg)) return(.stanParName(as.character(arg)))
  stop("cannot translate prior argument '", deparse1(arg),
       "' to Stan: arguments must be numeric literals or parameter names",
       call. = FALSE)
}

#' Look a prior call up in the lotri distribution catalogue
#'
#' @param priorStr the canonical deparsed prior string from `iniDf$prior`
#' @return list(name=, stanName=, kind=, support=, args=list of language,
#'   dist=the catalogue row)
#' @noRd
.stanPriorLookup <- function(priorStr) {
  .lang <- str2lang(priorStr)
  if (is.name(.lang)) .lang <- as.call(list(.lang)) # bare stdNormal
  .nm <- as.character(.lang[[1]])
  .tbl <- lotri::lotriPriorDists()
  .w <- which(.tbl$name == .nm | .tbl$camelName == .nm | .tbl$stanName == .nm)
  if (length(.w) != 1L) {
    stop("unknown prior distribution '", .nm, "'", call. = FALSE) # nocov
  }
  .row <- .tbl[.w, ]
  list(name = .nm, stanName = .row$stanName, kind = .row$kind,
       support = .row$support, args = as.list(.lang)[-1], dist = .row)
}

#' Apply the support-promotion rule
#'
#' @return c(lower=, upper=) effective bounds
#' @noRd
.stanPromoteSupport <- function(support, lower, upper, name) {
  if (support %in% c("positive", "nonneg")) {
    lower <- max(lower, 0)
  } else if (support == "unit") {
    lower <- max(lower, 0)
    upper <- min(upper, 1)
  } else if (support == "circular") {
    if ((is.finite(lower) && (lower < -pi || lower > pi)) ||
          (is.finite(upper) && (upper > pi || upper < -pi))) {
      stop("parameter '", name, "' has bounds outside [-pi, pi] but a circular",
           " prior", call. = FALSE)
    }
    if (!is.finite(lower)) lower <- -pi
    if (!is.finite(upper)) upper <- pi
  }
  c(lower = lower, upper = upper)
}

#' Emit one univariate prior
#'
#' @return list(constraint=, statement=, truncated=, lower=, upper=)
#' @noRd
.stanPriorUnivariate <- function(name, lk, lower, upper) {
  .b <- .stanPromoteSupport(lk$support, lower, upper, name)
  .lit <- vapply(lk$args, .stanArgIsLiteral, logical(1))
  .args <- vapply(lk$args, .stanArgDeparse, character(1))
  .argStr <- paste(.args, collapse = ", ")
  .par <- .stanParName(name)
  .bounded <- is.finite(.b[["lower"]]) || is.finite(.b[["upper"]])
  .truncated <- .bounded && !all(.lit)
  if (.truncated && lk$stanName %in% .stanNoLcdf) {
    stop("prior '", lk$name, "' on '", name, "' has parameter-dependent ",
         "arguments and finite bounds, but Stan has no ", lk$stanName,
         "_lcdf to normalize the truncation", call. = FALSE)
  }
  .stmt <- paste0(.par, " ~ ", lk$stanName, "(", .argStr, ")")
  if (.truncated) {
    .tl <- if (is.finite(.b[["lower"]])) .stanNum(.b[["lower"]]) else ""
    .tu <- if (is.finite(.b[["upper"]])) .stanNum(.b[["upper"]]) else ""
    .stmt <- paste0(.stmt, " T[", .tl, ", ", .tu, "]")
  }
  list(constraint = .stanConstraint(.b[["lower"]], .b[["upper"]]),
       statement = paste0(.stmt, ";"),
       truncated = .truncated,
       lower = .b[["lower"]], upper = .b[["upper"]])
}

#' Evaluate a multivariate prior argument by value (c() vectors and lotri()
#' matrices only; nothing user-parameter-dependent is allowed here)
#' @noRd
.stanEvalMvArg <- function(arg) {
  .env <- new.env(parent = baseenv())
  assign("lotri", lotri::lotri, envir = .env)
  .v <- tryCatch(eval(arg, envir = .env), error = function(e) NULL)
  if (is.null(.v) || !is.numeric(.v)) {
    stop("cannot evaluate multivariate prior argument '", deparse1(arg),
         "': it must be a numeric vector or a lotri() matrix",
         call. = FALSE)
  }
  .v
}

#' Stan source for a numeric vector / matrix literal
#' @noRd
.stanVecLit <- function(v) paste0("[", paste(.stanNum(v), collapse = ", "), "]'")
.stanMatLit <- function(m) {
  paste0("[", paste(apply(m, 1, function(r) {
    paste0("[", paste(.stanNum(r), collapse = ", "), "]")
  }), collapse = ", "), "]")
}

#' Emit one multivariate (normal-family) prior over a member set
#' @noRd
.stanPriorMultivariate <- function(members, lk) {
  .pars <- .stanParName(members)
  .lhs <- paste0("to_vector({", paste(.pars, collapse = ", "), "})")
  .args <- vapply(lk$args, function(a) {
    .v <- .stanEvalMvArg(a)
    if (is.matrix(.v)) .stanMatLit(.v) else .stanVecLit(.v)
  }, character(1))
  list(statement = paste0("target += ", lk$stanName, "_lpdf(", .lhs, " | ",
                          paste(.args, collapse = ", "), ");"))
}

#' Member set (in order) of a multivariate prior from its embedded lotri()
#' @noRd
.stanMvMembers <- function(lk, fallback) {
  for (a in lk$args) {
    .v <- tryCatch(.stanEvalMvArg(a), error = function(e) NULL)
    if (is.matrix(.v) && !is.null(dimnames(.v)[[1]])) return(dimnames(.v)[[1]])
  }
  fallback
}

#' Translate a model's `ini({})` priors into Stan
#'
#' Consumes [rxode2::rxUiPriors()] and [lotri::lotriPriorDists()] and returns
#' the Stan-side pieces: per-parameter bound declarations (support-promoted)
#' and sampling statements for population parameters, and the classified
#' omega-block priors for the model generator.
#'
#' @param ui an rxode2 ui (or model function)
#' @return a list:
#' \describe{
#' \item{pop}{data.frame over population-parameter priors: `name` (nlmixr2),
#'   `par` (Stan identifier), `prior`, `stanName`, `kind`, `lower`, `upper`
#'   (effective, support-promoted), `constraint`, `truncated`, `members`
#'   (comma-joined for a multivariate group, `NA` otherwise), `statement`}
#' \item{omega}{data.frame over omega-block priors: `name`, `neta1`, `neta2`,
#'   `prior`, `stanName`, `kind`}
#' }
#' @export
#' @author Matthew L Fidler
stanPriors <- function(ui) {
  .ui <- rxode2::assertRxUi(ui)
  .pri <- rxode2::rxUiPriors(.ui)
  .popRow <- function(name, par, prior, stanName, kind, lower, upper,
                      constraint, truncated, members, statement) {
    data.frame(name = name, par = par, prior = prior, stanName = stanName,
               kind = kind, lower = lower, upper = upper,
               constraint = constraint, truncated = truncated,
               members = members, statement = statement,
               stringsAsFactors = FALSE)
  }
  .pop <- .popRow(character(0), character(0), character(0), character(0),
                  character(0), numeric(0), numeric(0), character(0),
                  logical(0), character(0), character(0))
  .om <- data.frame(name = character(0), neta1 = numeric(0),
                    neta2 = numeric(0), prior = character(0),
                    stanName = character(0), kind = character(0),
                    stringsAsFactors = FALSE)
  if (nrow(.pri) == 0L) return(list(pop = .pop, omega = .om))
  .mvSeen <- character(0)
  for (.i in seq_len(nrow(.pri))) {
    .r <- .pri[.i, ]
    .lk <- .stanPriorLookup(.r$prior)
    if (.lk$kind == "discrete") {
      stop("prior '", .lk$name, "' on '", .r$name, "' is a discrete ",
           "distribution; every nlmixr2 parameter is real-valued",
           call. = FALSE)
    }
    if (!is.na(.r$neta1)) {
      # omega-block prior: classified here, emitted by the model generator
      # (the declared Stan parameter is whatever the prior is written on, so
      # Stan's own constraining transform supplies the only Jacobian needed)
      .om <- rbind(.om, data.frame(name = .r$name, neta1 = .r$neta1,
                                   neta2 = .r$neta2, prior = .r$prior,
                                   stanName = .lk$stanName, kind = .lk$kind,
                                   stringsAsFactors = FALSE))
      next
    }
    if (.lk$kind == "matrix") {
      stop("matrix prior '", .lk$name, "' on population parameter '",
           .r$name, "' is not supported", call. = FALSE) # nocov
    }
    if (.lk$kind == "multivariate") {
      # duplicated byte-identically on every member row: dedupe, keep order
      # from the embedded lotri() dimnames
      .members <- .stanMvMembers(.lk, .pri$name[.pri$prior == .r$prior])
      .key <- digest::digest(list(.r$prior, sort(.members)))
      if (.key %in% .mvSeen) next
      .mvSeen <- c(.mvSeen, .key)
      .mv <- .stanPriorMultivariate(.members, .lk)
      .pop <- rbind(.pop, .popRow(.r$name, .stanParName(.r$name), .r$prior,
                                  .lk$stanName, .lk$kind, .r$lower, .r$upper,
                                  .stanConstraint(.r$lower, .r$upper), FALSE,
                                  paste(.members, collapse = ","),
                                  .mv$statement))
      next
    }
    .u <- .stanPriorUnivariate(.r$name, .lk, .r$lower, .r$upper)
    .pop <- rbind(.pop, .popRow(.r$name, .stanParName(.r$name), .r$prior,
                                .lk$stanName, .lk$kind, .u$lower, .u$upper,
                                .u$constraint, .u$truncated, NA_character_,
                                .u$statement))
  }
  list(pop = .pop, omega = .om)
}
