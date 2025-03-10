# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compute the generalized binomial test statistic
#' 
#' @noRd
compute_GBT <- function(s, p) {
    .Call(`_aberrance_compute_GBT`, s, p)
}

#' Compute the M4 statistic
#' 
#' @noRd
compute_M4 <- function(s_1, s_0, p_1, p_0, p_n) {
    .Call(`_aberrance_compute_M4`, s_1, s_0, p_1, p_0, p_n)
}

