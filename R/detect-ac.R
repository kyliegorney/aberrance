#' Detect answer copying
#'
#' @description Detect answer copying for all possible source-copier pairs.
#'
#' @param method The answer copying statistic(s) to compute. Options for
#'   score-based statistics are:
#'   - `"OMG_S"` for the conditional \eqn{\omega} statistic (Wollack, 1997).
#'   - `"GBT_S"` for the conditional \eqn{GBT} statistic (van der Linden &
#'     Sotaridona, 2006).
#'
#'   Options for response-based statistics are:
#'   - `"OMG_R"` for the conditional \eqn{\omega} statistic (Wollack, 1997).
#'   - `"GBT_R"` for the conditional \eqn{GBT} statistic (van der Linden &
#'     Sotaridona, 2006).
#'
#' @param x,r Matrices of raw data. `x` is for the item scores and `r` the item
#'   responses.
#'
#' @inheritParams detect_pm
#'
#' @returns A list is returned with the following elements:
#' \item{stat}{A matrix of answer copying statistics.}
#' \item{pval}{A matrix of *p*-values.}
#' \item{flag}{An array of flagging results. The first dimension corresponds to
#'   source-copier pairs, the second dimension to methods, and the third
#'   dimension to significance levels.}
#'
#' @references
#' van der Linden, W. J., & Sotaridona, L. (2006). Detecting answer copying when
#' the regular response process follows a known response model. *Journal of
#' Educational and Behavioral Statistics*, *31*(3), 283--304.
#'
#' Wollack, J. A. (1997). A nominal response model approach for detecting answer
#' copying. *Applied Psychological Measurement*, *21*(4), 307--320.
#'
#' @seealso [detect_as()] to detect answer similarity.
#'
#' @examples
#' # Setup for Examples 1 and 2 ------------------------------------------------
#'
#' # Settings
#' set.seed(0)     # seed for reproducibility
#' N <- 50         # number of persons
#' n <- 40         # number of items
#'
#' # Randomly select 10% sources and 10% copiers
#' s <- sample(1:N, size = N * 0.10)
#' c <- sample(setdiff(1:N, s), size = N * 0.10)
#'
#' # Create vector of indicators (1 = copying pair, 0 = non-copying pair)
#' pair <- t(combn(N, 2))
#' pair <- rbind(pair, pair[, 2:1])
#' ind <- ifelse(1:nrow(pair) %in% apply(
#'   rbind(cbind(s, c), cbind(c, s)), 1, function(p)
#'   which(pair[, 1] == p[1] & pair[, 2] == p[2])), 1, 0)
#' names(ind) <- paste(pair[, 1], pair[, 2], sep = "-")
#'
#' # Example 1: Item Scores ----------------------------------------------------
#'
#' # Generate person parameters for the 3PL model
#' xi <- cbind(theta = rnorm(N, mean = 0.00, sd = 1.00))
#'
#' # Generate item parameters for the 3PL model
#' psi <- cbind(
#'   a = rlnorm(n, meanlog = 0.00, sdlog = 0.25),
#'   b = rnorm(n, mean = 0.00, sd = 1.00),
#'   c = runif(n, min = 0.05, max = 0.30)
#' )
#'
#' # Simulate uncontaminated data
#' x <- sim(psi, xi)$x
#'
#' # Modify contaminated data by replacing 40% of the copier scores with source
#' # scores
#' for (v in 1:length(c)) {
#'   ci <- sample(1:n, size = n * 0.40)
#'   x[c[v], ci] <- x[s[v], ci]
#' }
#'
#' # Detect answer copying
#' out <- detect_ac(
#'   method = c("OMG_S", "GBT_S"),
#'   psi = psi,
#'   x = x
#' )
#'
#' # Example 2: Item Responses -------------------------------------------------
#'
#' # Generate person parameters for the nominal response model
#' xi <- cbind(eta = rnorm(N, mean = 0.00, sd = 1.00))
#'
#' # Generate item parameters for the nominal response model
#' psi <- cbind(
#'   lambda1 = rnorm(n, mean = -0.50, sd = 0.50),
#'   lambda2 = rnorm(n, mean = -0.50, sd = 0.50),
#'   lambda3 = rnorm(n, mean = -0.50, sd = 0.50),
#'   lambda4 = rnorm(n, mean = 1.50, sd = 0.50),
#'   zeta1 = rnorm(n, mean = -0.50, sd = 0.50),
#'   zeta2 = rnorm(n, mean = -0.50, sd = 0.50),
#'   zeta3 = rnorm(n, mean = -0.50, sd = 0.50),
#'   zeta4 = rnorm(n, mean = 1.50, sd = 0.50)
#' )
#'
#' # Simulate uncontaminated data
#' r <- sim(psi, xi)$r
#'
#' # Modify contaminated data by replacing 40% of the copier responses with
#' # source responses
#' for (v in 1:length(c)) {
#'   ci <- sample(1:n, size = n * 0.40)
#'   r[c[v], ci] <- r[s[v], ci]
#' }
#'
#' # Detect answer copying
#' out <- detect_ac(
#'   method = c("OMG_R", "GBT_R"),
#'   psi = psi,
#'   r = r
#' )
#' @export

detect_ac <- function(method,
                      psi,
                      xi = NULL,
                      x = NULL, r = NULL,
                      interval = c(-4, 4),
                      alpha = 0.05) {

  # Checks
  if (any("S" %in% extract(method, 2)) && ("R" %in% extract(method, 2))) {
    stop("`method` may contain either score-based statistics or ",
         "response-based statistics, but not both.", call. = FALSE)
  }
  if (any("S" %in% extract(method, 2))) {
    check_par("x", psi, xi)
  } else if ("R" %in% extract(method, 2)) {
    check_par("r", psi, xi)
  }
  method <- match.arg(
    arg = unique(method),
    choices = c("OMG_S", "GBT_S", "OMG_R", "GBT_R"),
    several.ok = TRUE
  )
  check_data(x, r)

  # Setup
  N <- max(nrow(x), nrow(r))
  n <- max(ncol(x), ncol(r))
  pair <- t(combn(N, 2))
  NN <- nrow(pair)
  pair <- rbind(pair, pair[, 2:1])
  stat <- pval <- matrix(
    nrow = NN * 2, ncol = length(method),
    dimnames = list(
      pair = paste(pair[, 1], pair[, 2], sep = "-"),
      method = method
    )
  )
  flag <- array(
    dim = c(NN * 2, length(method), length(alpha)),
    dimnames = list(
      pair = row.names(stat),
      method = method,
      alpha = alpha
    )
  )

  # Estimate person parameters
  if (is.null(xi)) {
    xi <- est(interval, psi, x = x, r = r)
  }

  # Compute score-based statistics
  if (any(c("OMG_S", "GBT_S") %in% method)) {
    m <- count(psi, ignore = "lambda1")
    p_mat <- irt_p(m, psi, xi, ignore = "lambda1")
    for (v in 1:NN) {
      s <- as.integer(x[pair[v, 1], ] == x[pair[v, 2], ])
      if (all(!is.na(s))) {
        p <- q <- rep(NA, times = n)
        for (i in 1:n) {
          p[i] <- p_mat[pair[v, 2], i, x[pair[v, 1], i] + 1]
          q[i] <- p_mat[pair[v, 1], i, x[pair[v, 2], i] + 1]
        }
        if ("OMG_S" %in% method) {
          stat[v, "OMG_S"] <- compute_OMG(s, p)
          stat[v + NN, "OMG_S"] <- compute_OMG(s, q)
        }
        if ("GBT_S" %in% method) {
          stat[v, "GBT_S"] <- compute_GBT(s, p)
          stat[v + NN, "GBT_S"] <- compute_GBT(s, q)
        }
      }
    }
  }

  # Compute response-based statistics
  if (any(c("OMG_R", "GBT_R") %in% method)) {
    m <- count(psi)
    p_mat <- irt_p(m, psi, xi)
    for (v in 1:NN) {
      s <- as.integer(r[pair[v, 1], ] == r[pair[v, 2], ])
      if (all(!is.na(s))) {
        p <- q <- rep(NA, times = n)
        for (i in 1:n) {
          p[i] <- p_mat[pair[v, 2], i, r[pair[v, 1], i]]
          q[i] <- p_mat[pair[v, 1], i, r[pair[v, 2], i]]
        }
        if ("OMG_R" %in% method) {
          stat[v, "OMG_R"] <- compute_OMG(s, p)
          stat[v + NN, "OMG_R"] <- compute_OMG(s, q)
        }
        if ("GBT_R" %in% method) {
          stat[v, "GBT_R"] <- compute_GBT(s, p)
          stat[v + NN, "GBT_R"] <- compute_GBT(s, q)
        }
      }
    }
  }

  # Compute p-values
  pval[, grep("OMG", colnames(pval))] <-
    pnorm(stat[, grep("OMG", colnames(stat))], lower.tail = FALSE)
  pval[, grep("GBT", colnames(pval))] <-
    stat[, grep("GBT", colnames(stat))]

  # Compute flagging rates
  for (a in 1:length(alpha)) {
    flag[, , a] <- pval <= alpha[a]
  }

  # Output
  list(stat = stat, pval = pval, flag = flag)
}
