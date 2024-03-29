
#' Options for influence calculation
#'
#' @param ... Dispatched internally to \code{\link{set_compute}}.
#' @param p_max,n_max Numeric scalar with the maximum percentage / number of
#' removals to consider.
#' @param sm_re Integer scalar with the number of steps until inverses are
#' recalculated exactly (instead of using the Sherman-Morrison update).
#' @param fwl Optional integer scalar with the position of variables to
#' marginalise out.
#' @param fwl_re Optional integer scalar with the steps until variables are
#' marginalised out again.
#' @param adaptive Logical scalar setting adaptive set identificaiton.
#'
#' @return Returns a list with options.
#'
#' @export
set_options <- function(...,
  p_max = NULL, n_max = NULL,
  sm_re = 1L, fwl = 0L, fwl_re = 1L,
  adaptive = TRUE) {

  out <- set_compute(...)

  if(is.null(p_max)) {
    if(is.null(n_max)) {
      p_max <- .5
      n_max <- 100L
    } else {
      p_max = 1
    }
  } else if(is.null(n_max)) {
    n_max <- Inf
  }

  out$p_max = num_check(p_max, 0, 1, msg = "Choose a valid maximum percentage.")
  out$n_max = int_check(n_max, 1L, Inf, msg = "Choose a valid maximum number.")
  out$sm_re = int_check(sm_re, 0L, 1e6L,
    msg = "Choose a valid step size for Sherman-Morrison updates.")
  out$fwl = vapply(fwl, int_check, integer(1L), 0L, Inf,
    msg = "Choose valid indices for variables to retain after FWL.")
  out$fwl_re = int_check(fwl_re, 1L, 1e6L,
    msg = "Choose a valid step size for recalculating FWL.")
  out$adaptive <- isTRUE(adaptive)

  return(out)
}


#' Options for influence computation
#'
#' @param x Character scalar with short-hand options.
#' @param hat,beta,sigma,se,tstat,cooksd,dffits,rstudent,covration,cluster
#' Logical scalars determining which quantities to compute (exactly).
#'
#' @return Returns a list with options.
#'
#' @export
set_compute <- function(x = c("some", "none", "all"),
  hat, beta, sigma, se, tstat, cooksd, dffits, rstudent, covratio, cluster) {

  x <- match.arg(x)

  out <- if(x == "all") {
    list("hat" = TRUE, "beta" = TRUE, "sigma" = TRUE,
      "se" = TRUE, "tstat" = TRUE,
      "cooksd" = TRUE, "dffits" = TRUE, "rstudent" = TRUE, "covratio" = TRUE,
      "cluster" = TRUE)
  } else if(x == "some") {
    list("hat" = TRUE, "beta" = TRUE, "sigma" = TRUE,
      "se" = FALSE, "tstat" = FALSE, "cooksd" = FALSE, "dffits" = FALSE,
      "rstudent" = FALSE, "covratio" = FALSE, "cluster" = FALSE)
  } else if(x == "none") {
    list("hat" = FALSE, "beta" = FALSE, "sigma" = FALSE,
      "se" = FALSE, "tstat" = FALSE, "cooksd" = FALSE, "dffits" = FALSE,
      "rstudent" = FALSE, "covratio" = FALSE, "cluster" = FALSE)
  }

  if(!missing(hat)) {out$hat <- isTRUE(hat)}
  if(!missing(beta)) {out$beta <- isTRUE(beta)}
  if(!missing(sigma)) {out$sigma <- isTRUE(sigma)}
  if(!missing(se)) {out$se <- isTRUE(se)}
  if(!missing(tstat)) {out$tstat <- isTRUE(tstat)}
  if(!missing(cooksd)) {out$cooksd <- isTRUE(cooksd)}
  if(!missing(dffits)) {out$dffits <- isTRUE(dffits)}
  if(!missing(rstudent)) {out$rstudent <- isTRUE(rstudent)}
  if(!missing(covratio)) {out$covratio <- isTRUE(covratio)}
  if(!missing(cluster)) {out$cluster <- isTRUE(cluster)}

  if(out$se) {out$beta <- out$sigma <- TRUE}
  if(out$tstat) {out$se <- out$beta <- out$sigma <- TRUE}
  if(out$cooksd || out$dffits || out$rstudent || out$covratio) {out$hat <- TRUE}
  if((out$dffits || out$rstudent || out$covratio )&& !out$sigma) {
    warning("Consider setting `sigma=TRUE` for exact computation.")
  }

  return(out)
}
