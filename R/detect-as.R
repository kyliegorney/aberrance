#' Detect answer similarity
#'
#' @description Detect answer similarity for all possible pairs.
#'
#' @param method The answer similarity statistic(s) to compute. Options for
#'   score-based statistics are:
#'   - `"OMG_S"` for the unconditional \eqn{\omega} statistic (Romero et al.,
#'     2015).
#'   - `"GBT_S"` for the unconditional \eqn{GBT} statistic (van der Linden &
#'     Sotaridona, 2006).
#'   - `"M4_S"` for the \eqn{M4} statistic (Maynes, 2014).
#'
#'   Options for score and distractor-based statistics are:
#'   - `"OMG_SD"` for the unconditional \eqn{\omega} statistic (Romero et al.,
#'     2015).
#'   - `"GBT_SD"` for the unconditional \eqn{GBT} statistic (van der Linden &
#'     Sotaridona, 2006).
#'   - `"M4_SD"` for the \eqn{M4} statistic (Maynes, 2014).
#'
#'   Options for response-based statistics are:
#'   - `"OMG_R"` for the unconditional \eqn{\omega} statistic (Romero et al.,
#'     2015).
#'   - `"GBT_R"` for the unconditional \eqn{GBT} statistic (van der Linden &
#'     Sotaridona, 2006).
#'   - `"M4_R"` for the \eqn{M4} statistic (Maynes, 2014).
#'
#'   Options for score and response time-based statistics are:
#'   - `"OMG_ST"` for the unconditional \eqn{\omega} statistic (Gorney &
#'     Wollack, 2024).
#'   - `"GBT_ST"` for the unconditional \eqn{GBT} statistic (Gorney & Wollack,
#'     2024).
#'
#'   Options for score, distractor, and response time-based statistics are:
#'   - `"OMG_SDT"` for the unconditional \eqn{\omega} statistic (Gorney &
#'     Wollack, 2024).
#'   - `"GBT_SDT"` for the unconditional \eqn{GBT} statistic (Gorney & Wollack,
#'     2024).
#'
#'   Options for response and response time-based statistics are:
#'   - `"OMG_RT"` for the unconditional \eqn{\omega} statistic (Gorney &
#'     Wollack, 2024).
#'   - `"GBT_RT"` for the unconditional \eqn{GBT} statistic (Gorney & Wollack,
#'     2024).
#'
#' @inheritParams detect_pm
#'
#' @returns A list is returned with the following elements:
#' \item{stat}{A matrix of answer similarity statistics.}
#' \item{pval}{A matrix of *p*-values.}
#' \item{flag}{An array of flagging results. The first dimension corresponds to
#'   pairs, the second dimension to methods, and the third dimension to
#'   significance levels.}
#'
#' @references
#' Gorney, K., & Wollack, J. A. (2024). Using response times in answer
#' similarity analysis. *Journal of Educational and Behavioral Statistics*.
#' Advance online publication.
#'
#' Maynes, D. (2014). Detection of non-independent test taking by similarity
#' analysis. In N. M. Kingston & A. K. Clark (Eds.), *Test fraud: Statistical
#' detection and methodology* (pp. 53--80). Routledge.
#'
#' Romero, M., Riascos, Á., & Jara, D. (2015). On the optimality of
#' answer-copying indices: Theory and practice. *Journal of Educational and
#' Behavioral Statistics*, *40*(5), 435--453.
#'
#' van der Linden, W. J., & Sotaridona, L. (2006). Detecting answer copying when
#' the regular response process follows a known response model. *Journal of
#' Educational and Behavioral Statistics*, *31*(3), 283--304.
#'
#' @seealso
#' [detect_ac()] to detect answer copying.
#'
#' [detect_pk()] to detect preknowledge.
#'
#' @examples
#' # Setup for Examples 1 and 2 ------------------------------------------------
#'
#' # Settings
#' set.seed(0)     # seed for reproducibility
#' N <- 50         # number of persons
#' n <- 40         # number of items
#'
#' # Randomly select 10% examinees with preknowledge and 40% compromised items
#' cv <- sample(1:N, size = N * 0.10)
#' ci <- sample(1:n, size = n * 0.40)
#'
#' # Create vector of indicators (1 = similar pair, 0 = non-similar pair)
#' pair <- t(combn(N, 2))
#' ind <- ifelse((pair[, 1] %in% cv) & (pair[, 2] %in% cv), 1, 0)
#' names(ind) <- paste(pair[, 1], pair[, 2], sep = "-")
#'
#' # Example 1: Item Scores and Response Times ---------------------------------
#'
#' # Generate person parameters for the 3PL model and lognormal model
#' xi <- MASS::mvrnorm(
#'   N,
#'   mu = c(theta = 0.00, tau = 0.00),
#'   Sigma = matrix(c(1.00, 0.25, 0.25, 0.25), ncol = 2)
#' )
#'
#' # Generate item parameters for the 3PL model and lognormal model
#' psi <- cbind(
#'   a = rlnorm(n, meanlog = 0.00, sdlog = 0.25),
#'   b = NA,
#'   c = runif(n, min = 0.05, max = 0.30),
#'   alpha = runif(n, min = 1.50, max = 2.50),
#'   beta = NA
#' )
#'
#' # Generate positively correlated difficulty and time intensity parameters
#' psi[, c("b", "beta")] <- MASS::mvrnorm(
#'   n,
#'   mu = c(b = 0.00, beta = 3.50),
#'   Sigma = matrix(c(1.00, 0.20, 0.20, 0.15), ncol = 2)
#' )
#'
#' # Simulate uncontaminated data
#' dat <- sim(psi, xi)
#' x <- dat$x
#' y <- dat$y
#'
#' # Modify contaminated data by changing the item scores and reducing the log
#' # response times
#' x[cv, ci] <- rbinom(length(cv) * length(ci), size = 1, prob = 0.90)
#' y[cv, ci] <- y[cv, ci] * 0.75
#'
#' # Detect answer similarity
#' out <- detect_as(
#'   method = c("OMG_S", "GBT_S", "OMG_ST", "GBT_ST"),
#'   psi = psi,
#'   x = x,
#'   y = y
#' )
#'
#' # Example 2: Polytomous Item Scores -----------------------------------------
#'
#' # Generate person parameters for the generalized partial credit model
#' xi <- cbind(theta = rnorm(N, mean = 0.00, sd = 1.00))
#'
#' # Generate item parameters for the generalized partial credit model
#' psi <- cbind(
#'   a = rlnorm(n, meanlog = 0.00, sdlog = 0.25),
#'   c0 = 0,
#'   c1 = rnorm(n, mean = -1.00, sd = 0.50),
#'   c2 = rnorm(n, mean = 0.00, sd = 0.50),
#'   c3 = rnorm(n, mean = 1.00, sd = 0.50)
#' )
#'
#' # Simulate uncontaminated data
#' x <- sim(psi, xi)$x
#'
#' # Modify contaminated data by changing the item scores to the maximum score
#' x[cv, ci] <- 3
#'
#' # Detect answer similarity
#' out <- detect_as(
#'   method = c("OMG_S", "GBT_S"),
#'   psi = psi,
#'   x = x
#' )
#'
#' # Setup for Examples 3 and 4 ------------------------------------------------
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
#' # Create vector of indicators (1 = similar pair, 0 = non-similar pair)
#' pair <- t(combn(N, 2))
#' ind <- ifelse(1:nrow(pair) %in% apply(
#'   rbind(cbind(s, c), cbind(c, s)), 1, function(p)
#'   which(pair[, 1] == p[1] & pair[, 2] == p[2])), 1, 0)
#' names(ind) <- paste(pair[, 1], pair[, 2], sep = "-")
#'
#' # Example 3: Item Scores and Distractors ------------------------------------
#'
#' # Generate person parameters for the nested logit model
#' xi <- MASS::mvrnorm(
#'   N,
#'   mu = c(theta = 0.00, eta = 0.00),
#'   Sigma = matrix(c(1.00, 0.80, 0.80, 1.00), ncol = 2)
#' )
#'
#' # Generate item parameters for the nested logit model
#' psi <- cbind(
#'   a = rlnorm(n, meanlog = 0.00, sdlog = 0.25),
#'   b = rnorm(n, mean = 0.00, sd = 1.00),
#'   c = runif(n, min = 0.05, max = 0.30),
#'   lambda1 = rnorm(n, mean = 0.00, sd = 1.00),
#'   lambda2 = rnorm(n, mean = 0.00, sd = 1.00),
#'   lambda3 = rnorm(n, mean = 0.00, sd = 1.00),
#'   zeta1 = rnorm(n, mean = 0.00, sd = 1.00),
#'   zeta2 = rnorm(n, mean = 0.00, sd = 1.00),
#'   zeta3 = rnorm(n, mean = 0.00, sd = 1.00)
#' )
#'
#' # Simulate uncontaminated data
#' dat <- sim(psi, xi)
#' x <- dat$x
#' d <- dat$d
#'
#' # Modify contaminated data by replacing 40% of the copier scores and
#' # distractors with source scores and distractors
#' for (v in 1:length(c)) {
#'   ci <- sample(1:n, size = n * 0.40)
#'   x[c[v], ci] <- x[s[v], ci]
#'   d[c[v], ci] <- d[s[v], ci]
#' }
#'
#' # Detect answer similarity
#' out <- detect_as(
#'   method = c("OMG_S", "GBT_S", "OMG_SD", "GBT_SD"),
#'   psi = psi,
#'   x = x,
#'   d = d
#' )
#'
#' # Example 4: Item Responses -------------------------------------------------
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
#' # Detect answer similarity
#' out <- detect_as(
#'   method = c("OMG_R", "GBT_R"),
#'   psi = psi,
#'   r = r
#' )
#' @export

detect_as <- function(method,
                      psi,
                      xi = NULL,
                      x = NULL, d = NULL, r = NULL, y = NULL,
                      interval = c(-4, 4),
                      alpha = 0.05) {

  # Checks
  if (any(c("S", "SD", "ST", "SDT") %in% extract(method, 2)) &&
      any(c("R", "RT") %in% extract(method, 2))) {
    stop("`method` may contain either score-based statistics or ",
         "response-based statistics, but not both.", call. = FALSE)
  }
  if (any(c("S", "SD", "ST", "SDT") %in% extract(method, 2))) {
    check_par("x", psi, xi)
    if (any(c("SD", "SDT") %in% extract(method, 2))) {
      check_par("d", psi, xi)
    }
  } else if (any(c("R", "RT") %in% extract(method, 2))) {
    check_par("r", psi, xi)
  }
  if (any(c("ST", "SDT", "RT") %in% extract(method, 2))) {
    check_par("y", psi)
  }
  method <- match.arg(
    arg = unique(method),
    choices = c("OMG_S", "GBT_S", "M4_S",
                "OMG_SD", "GBT_SD", "M4_SD",
                "OMG_R", "GBT_R", "M4_R",
                "OMG_ST", "GBT_ST",
                "OMG_SDT", "GBT_SDT",
                "OMG_RT", "GBT_RT"),
    several.ok = TRUE
  )
  check_data(x, d, r, y)

  # Setup
  N <- max(nrow(x), nrow(d), nrow(r), nrow(y))
  n <- max(ncol(x), ncol(d), ncol(r), ncol(y))
  pair <- t(combn(N, 2))
  NN <- nrow(pair)
  stat <- pval <- matrix(
    nrow = NN, ncol = length(method),
    dimnames = list(
      pair = paste(pair[, 1], pair[, 2], sep = "-"),
      method = method
    )
  )
  flag <- array(
    dim = c(NN, length(method), length(alpha)),
    dimnames = list(
      pair = row.names(stat),
      method = method,
      alpha = alpha
    )
  )

  # Estimate person parameters
  if (is.null(xi)) {
    xi <- est(interval, psi, x = x, d = d, r = r, y = y)
  }

  # Compute score-based statistics
  if (any(c("OMG_S", "GBT_S", "M4_S") %in% method)) {
    m <- count(psi, ignore = "lambda1")
    p_mat <- irt_p(m, psi, xi, ignore = "lambda1")
    for (v in 1:NN) {
      s <- sum(x[pair[v, 1], ] == x[pair[v, 2], ])
      p <- rowSums(p_mat[pair[v, 1], , ] * p_mat[pair[v, 2], , ], na.rm = TRUE)
      if ("OMG_S" %in% method) {
        stat[v, "OMG_S"] <- compute_OMG(s, p)
      }
      if ("GBT_S" %in% method) {
        stat[v, "GBT_S"] <- compute_GBT(s, p)
      }
      if ("M4_S" %in% method) {
        if ("b" %in% colnames(psi)) {
          s_1 <- sum((x[pair[v, 1], ] == 1) & (x[pair[v, 2], ] == 1))
          p_1 <- p_mat[pair[v, 1], , 2] * p_mat[pair[v, 2], , 2]
        } else {
          s_1 <- sum((x[pair[v, 1], ] == (m - 1)) &
                       (x[pair[v, 2], ] == (m - 1)))
          p_1 <- rep(NA, times = n)
          for (i in 1:n) {
            p_1[i] <- p_mat[pair[v, 1], i, m[i]] * p_mat[pair[v, 2], i, m[i]]
          }
        }
        s_0 <- s - s_1
        p_0 <- p - p_1
        stat[v, "M4_S"] <- compute_M4(c(s_1, s_0), cbind(p_1, p_0, 1 - p))
      }
    }
  }

  # Compute score and distractor-based statistics
  if (any(c("OMG_SD", "GBT_SD", "M4_SD") %in% method)) {
    m <- count(psi)
    p_mat <- irt_p(m, psi, xi)
    for (v in 1:NN) {
      s_1 <- sum((x[pair[v, 1], ] == 1) & (x[pair[v, 2], ] == 1))
      s_0 <- sum(d[pair[v, 1], ] == d[pair[v, 2], ], na.rm = TRUE)
      s <- s_1 + s_0
      p <- rowSums(p_mat[pair[v, 1], , ] * p_mat[pair[v, 2], , ], na.rm = TRUE)
      if ("OMG_SD" %in% method) {
        stat[v, "OMG_SD"] <- compute_OMG(s, p)
      }
      if ("GBT_SD" %in% method) {
        stat[v, "GBT_SD"] <- compute_GBT(s, p)
      }
      if ("M4_SD" %in% method) {
        p_1 <- rep(NA, times = n)
        for (i in 1:n) {
          p_1[i] <- p_mat[pair[v, 1], i, m[i]] * p_mat[pair[v, 2], i, m[i]]
        }
        p_0 <- p - p_1
        stat[v, "M4_SD"] <- compute_M4(c(s_1, s_0), cbind(p_1, p_0, 1 - p))
      }
    }
  }

  # Compute response-based statistics
  if (any(c("OMG_R", "GBT_R", "M4_R") %in% method)) {
    m <- count(psi)
    p_mat <- irt_p(m, psi, xi)
    for (v in 1:NN) {
      s <- sum(r[pair[v, 1], ] == r[pair[v, 2], ])
      p <- rowSums(p_mat[pair[v, 1], , ] * p_mat[pair[v, 2], , ], na.rm = TRUE)
      if ("OMG_R" %in% method) {
        stat[v, "OMG_R"] <- compute_OMG(s, p)
      }
      if ("GBT_R" %in% method) {
        stat[v, "GBT_R"] <- compute_GBT(s, p)
      }
      if ("M4_R" %in% method) {
        s_1 <- sum((r[pair[v, 1], ] == m) & (r[pair[v, 2], ] == m))
        p_1 <- rep(NA, times = n)
        for (i in 1:n) {
          p_1[i] <- p_mat[pair[v, 1], i, m[i]] * p_mat[pair[v, 2], i, m[i]]
        }
        s_0 <- s - s_1
        p_0 <- p - p_1
        stat[v, "M4_R"] <- compute_M4(c(s_1, s_0), cbind(p_1, p_0, 1 - p))
      }
    }
  }

  # Compute score and response time-based statistics
  if (any(c("OMG_ST", "GBT_ST") %in% method)) {
    m <- count(psi, ignore = "lambda1")
    p_mat <- irt_p(m, psi, xi, ignore = "lambda1")
    mu <- t(outer(psi[, "beta"], xi[, "tau"], "-"))
    for (v in 1:NN) {
      s <- sum((x[pair[v, 1], ] == x[pair[v, 2], ]) &
                 (y[pair[v, 1], ] < mu[pair[v, 1], ]) &
                 (y[pair[v, 2], ] < mu[pair[v, 2], ]))
      p <- 0.25 *
        rowSums(p_mat[pair[v, 1], , ] * p_mat[pair[v, 2], , ], na.rm = TRUE)
      if ("OMG_ST" %in% method) {
        stat[v, "OMG_ST"] <- compute_OMG(s, p)
      }
      if ("GBT_ST" %in% method) {
        stat[v, "GBT_ST"] <- compute_GBT(s, p)
      }
    }
  }

  # Compute score, distractor, and response time-based statistics
  if (any(c("OMG_SDT", "GBT_SDT") %in% method)) {
    m <- count(psi)
    p_mat <- irt_p(m, psi, xi)
    mu <- t(outer(psi[, "beta"], xi[, "tau"], "-"))
    for (v in 1:NN) {
      s_1 <- sum((x[pair[v, 1], ] == 1) & (x[pair[v, 2], ] == 1) &
                   (y[pair[v, 1], ] < mu[pair[v, 1], ]) &
                   (y[pair[v, 2], ] < mu[pair[v, 2], ]))
      s_0 <- sum(d[pair[v, 1], ] == d[pair[v, 2], ] &
                   (y[pair[v, 1], ] < mu[pair[v, 1], ]) &
                   (y[pair[v, 2], ] < mu[pair[v, 2], ]), na.rm = TRUE)
      s <- s_1 + s_0
      p <- 0.25 *
        rowSums(p_mat[pair[v, 1], , ] * p_mat[pair[v, 2], , ], na.rm = TRUE)
      if ("OMG_SDT" %in% method) {
        stat[v, "OMG_SDT"] <- compute_OMG(s, p)
      }
      if ("GBT_SDT" %in% method) {
        stat[v, "GBT_SDT"] <- compute_GBT(s, p)
      }
    }
  }

  # Compute response and response time-based statistics
  if (any(c("OMG_RT", "GBT_RT") %in% method)) {
    m <- count(psi)
    p_mat <- irt_p(m, psi, xi)
    mu <- t(outer(psi[, "beta"], xi[, "tau"], "-"))
    for (v in 1:NN) {
      s <- sum((r[pair[v, 1], ] == r[pair[v, 2], ]) &
                 (y[pair[v, 1], ] < mu[pair[v, 1], ]) &
                 (y[pair[v, 2], ] < mu[pair[v, 2], ]))
      p <- 0.25 *
        rowSums(p_mat[pair[v, 1], , ] * p_mat[pair[v, 2], , ], na.rm = TRUE)
      if ("OMG_RT" %in% method) {
        stat[v, "OMG_RT"] <- compute_OMG(s, p)
      }
      if ("GBT_RT" %in% method) {
        stat[v, "GBT_RT"] <- compute_GBT(s, p)
      }
    }
  }

  # Compute p-values
  pval[, grep("OMG", colnames(pval))] <-
    pnorm(stat[, grep("OMG", colnames(stat))], lower.tail = FALSE)
  pval[, grep("GBT", colnames(pval))] <-
    stat[, grep("GBT", colnames(stat))]
  pval[, grep("M4", colnames(pval))] <-
    stat[, grep("M4", colnames(stat))]

  # Compute flagging rates
  for (a in 1:length(alpha)) {
    flag[, , a] <- pval <= alpha[a]
  }

  # Output
  list(stat = stat, pval = pval, flag = flag)
}
