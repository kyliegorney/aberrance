#' Detect preknowledge
#'
#' @description Detect preknowledge under the assumption that the set of
#'   compromised items is known.
#'
#' @param method The preknowledge detection statistic(s) to compute. Options for
#'   score-based statistics are:
#'   - `"L_S"` for the signed likelihood ratio test statistic (Sinharay, 2017).
#'   - `"ML_S"` for the modified signed likelihood ratio test statistic
#'     (Sinharay & Jensen, 2019). For numerical stability, an absolute cutoff
#'     value can be specified using `cutoff`. *Note:* This statistic cannot be
#'     computed under the 3PL model or the graded response model.
#'   - `"LR_S"` for the Lugannani-Rice approximation (Sinharay & Jensen, 2019).
#'     For numerical stability, an absolute cutoff value can be specified using
#'     `cutoff`. *Note:* This statistic cannot be computed under the 3PL model
#'     or the graded response model.
#'   - `"S_S"` for the signed score test statistic (Sinharay, 2017).
#'   - `"W_S"` for the Wald test statistic (Sinharay & Jensen, 2019).
#'
#'   Options for response time-based statistics are:
#'   - `"L_T"` for the signed likelihood ratio test statistic, or equivalently,
#'     `"W_T"` for the Wald test statistic (Sinharay, 2020).
#'
#'   Options for score and response time-based statistics are:
#'   - `"L_ST"` for the constrained likelihood ratio test statistic (Sinharay &
#'     Johnson, 2020).
#'
#' @param ci A vector of compromised item positions. All other items are
#'   presumed secure.
#'
#' @param xi,xi_c,xi_s Matrices of person parameters. `xi` is based on all
#'   items, `xi_c` is based on the compromised items, and `xi_s` is based on the
#'   secure items. If `NULL` (default), person parameters are estimated using
#'   maximum likelihood estimation.
#'
#' @param x,y Matrices of raw data. `x` is for the item scores and `y` the item
#'   log response times.
#'
#' @param cutoff Use with the modified signed likelihood ratio test statistic
#'   and the Lugannani-Rice approximation. If the absolute value of the signed
#'   likelihood ratio test statistic is less than the cutoff (default is
#'   `0.05`), then the modified signed likelihood ratio test statistic is
#'   replaced with the signed likelihood ratio test statistic and the
#'   Lugannani-Rice approximation is replaced with the \eqn{p}-value of the
#'   signed likelihood ratio test statistic.
#'
#' @inheritParams detect_pm
#'
#' @returns A list is returned with the following elements:
#' \item{stat}{A matrix of preknowledge detection statistics.}
#' \item{pval}{A matrix of *p*-values.}
#' \item{flag}{An array of flagging results. The first dimension corresponds to
#'   persons, the second dimension to methods, and the third dimension to
#'   significance levels.}
#'
#' @references
#' Sinharay, S. (2017). Detection of item preknowledge using likelihood ratio
#' test and score test. *Journal of Educational and Behavioral Statistics*,
#' *42*(1), 46--68.
#'
#' Sinharay, S. (2020). Detection of item preknowledge using response times.
#' *Applied Psychological Measurement*, *44*(5), 376--392.
#'
#' Sinharay, S., & Jensen, J. L. (2019). Higher-order asymptotics and its
#' application to testing the equality of the examinee ability over two sets of
#' items. *Psychometrika*, *84*(2), 484--510.
#'
#' Sinharay, S., & Johnson, M. S. (2020). The use of item scores and response
#' times to detect examinees who may have benefited from item preknowledge.
#' *British Journal of Mathematical and Statistical Psychology*, *73*(3),
#' 397--419.
#'
#' @seealso [detect_as()] to detect answer similarity.
#'
#' @examples
#' # Setup for Examples 1 and 2 ------------------------------------------------
#'
#' # Settings
#' set.seed(0)     # seed for reproducibility
#' N <- 500        # number of persons
#' n <- 40         # number of items
#'
#' # Randomly select 10% examinees with preknowledge and 40% compromised items
#' cv <- sample(1:N, size = N * 0.10)
#' ci <- sample(1:n, size = n * 0.40)
#'
#' # Create vector of indicators (1 = preknowledge, 0 = no preknowledge)
#' ind <- ifelse(1:N %in% cv, 1, 0)
#'
#' # Example 1: Item Scores and Response Times ---------------------------------
#'
#' # Generate person parameters for the 2PL model and lognormal model
#' xi <- MASS::mvrnorm(
#'   N,
#'   mu = c(theta = 0.00, tau = 0.00),
#'   Sigma = matrix(c(1.00, 0.25, 0.25, 0.25), ncol = 2)
#' )
#'
#' # Generate item parameters for the 2PL model and lognormal model
#' psi <- cbind(
#'   a = rlnorm(n, meanlog = 0.00, sdlog = 0.25),
#'   b = NA,
#'   c = 0,
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
#' # Detect preknowledge
#' out <- detect_pk(
#'   method = c("L_S", "ML_S", "LR_S", "S_S", "W_S", "L_T", "L_ST"),
#'   ci = ci,
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
#' # Detect preknowledge
#' out <- detect_pk(
#'   method = c("L_S", "ML_S", "LR_S", "S_S", "W_S"),
#'   ci = ci,
#'   psi = psi,
#'   x = x
#' )
#' @export

detect_pk <- function(method,
                      ci,
                      psi,
                      xi = NULL, xi_c = NULL, xi_s = NULL,
                      x = NULL, y = NULL,
                      interval = c(-4, 4),
                      alpha = 0.05,
                      cutoff = 0.05) {

  # Checks
  si <- setdiff(1:nrow(psi), ci)
  if (length(ci) < 1) {
    stop("At least one item must be specified as compromised.", call. = FALSE)
  } else if (length(si) < 1) {
    stop("At least one item must be specified as secure.", call. = FALSE)
  }
  if (any(c("L_S", "ML_S", "LR_S", "S_S", "W_S", "L_ST") %in% method)) {
    check_par("x", psi, xi)
    if (any(c("ML_S", "LR_S") %in% method)) {
      if ((("b" %in% colnames(psi)) && any(psi[, "c"] > 0)) ||
          ("b1" %in% colnames(psi))) {
        method <- setdiff(method, c("ML_S", "LR_S"))
        warning("The ML_S and LR_S statistics cannot be computed under the ",
                "3PL model or the graded response model.", call. = FALSE)
      }
    }
  }
  if (any(c("L_T", "W_T", "L_ST") %in% method)) {
    check_par("y", psi, xi)
  }
  method <- match.arg(
    arg = unique(method),
    choices = c("L_S", "ML_S", "LR_S", "S_S", "W_S", "L_T", "W_T", "L_ST"),
    several.ok = TRUE
  )
  check_data(x, y)

  # Setup
  N <- max(nrow(x), nrow(y))
  n <- max(ncol(x), ncol(y))
  stat <- pval <- matrix(
    nrow = N, ncol = length(method),
    dimnames = list(
      person = 1:N,
      method = method
    )
  )
  flag <- array(
    dim = c(N, length(method), length(alpha)),
    dimnames = list(
      person = 1:N,
      method = method,
      alpha = alpha
    )
  )

  # Estimate person parameters
  if (is.null(xi)) {
    xi <- est(interval, psi, x = x, y = y)
  }
  if (is.null(xi_c)) {
    xi_c <- est(
      interval,
      psi[ci, , drop = FALSE],
      x = x[, ci, drop = FALSE],
      y = y[, ci, drop = FALSE]
    )
  }
  if (is.null(xi_s)) {
    xi_s <- est(
      interval,
      psi[si, , drop = FALSE],
      x = x[, si, drop = FALSE],
      y = y[, si, drop = FALSE]
    )
  }

  # Compute score-based statistics
  if (any(c("L_S", "ML_S", "LR_S", "S_S", "W_S", "L_ST") %in% method)) {
    m <- count(psi, ignore = "lambda1")
    p_0 <- irt_p(m, psi, xi, ignore = "lambda1")
    p1_0 <- irt_p1(p_0, m, psi, xi, ignore = "lambda1")
    if (any(c("L_S", "ML_S", "LR_S", "L_ST") %in% method)) {
      p_1 <- p_0
      p_1[, ci, ] <- irt_p(m[ci], psi[ci, , drop = FALSE], xi_c,
                           ignore = "lambda1")
      p_1[, si, ] <- irt_p(m[si], psi[si, , drop = FALSE], xi_s,
                           ignore = "lambda1")
      L_S <- compute_L(x, p_0, p_1, xi_c[, "theta"], xi_s[, "theta"])
      if ("L_S" %in% method) {
        stat[, "L_S"] <- L_S
      }
      if (any(c("ML_S", "LR_S") %in% method)) {
        p1_1 <- p1_0
        p1_1[, ci, ] <- irt_p1(p_1[, ci, , drop = FALSE], m[ci],
                               psi[ci, , drop = FALSE], xi_c,
                               ignore = "lambda1")
        p1_1[, si, ] <- irt_p1(p_1[, si, , drop = FALSE], m[si],
                               psi[si, , drop = FALSE], xi_s,
                               ignore = "lambda1")
        num <- irt_info(p_1[, ci, , drop = FALSE], p1_1[, ci, , drop = FALSE]) *
          irt_info(p_1[, si, , drop = FALSE], p1_1[, si, , drop = FALSE])
        den <- irt_info(p_0[, ci, , drop = FALSE], p1_0[, ci, , drop = FALSE]) +
          irt_info(p_0[, si, , drop = FALSE], p1_0[, si, , drop = FALSE])
        q <- (xi_c[, "theta"] - xi_s[, "theta"]) * sqrt(num / den)
        if ("ML_S" %in% method) {
          stat[, "ML_S"] <- ifelse(abs(L_S) < cutoff, L_S,
                                   L_S + log(q / L_S) / L_S)
        }
        if ("LR_S" %in% method) {
          stat[, "LR_S"] <- ifelse(abs(L_S) < cutoff, pnorm(L_S),
                                   pnorm(L_S) + (1 / L_S - 1 / q) * dnorm(L_S))
        }
      }
    }
    if (any(c("S_S", "W_S") %in% method)) {
      info_c <- irt_info(p_0[, ci, , drop = FALSE], p1_0[, ci, , drop = FALSE])
      info_s <- irt_info(p_0[, si, , drop = FALSE], p1_0[, si, , drop = FALSE])
      if ("S_S" %in% method) {
        l1_c <- irt_l1(
          x[, ci, drop = FALSE],
          p_0[, ci, , drop = FALSE],
          p1_0[, ci, , drop = FALSE]
        )
        l1_s <- irt_l1(
          x[, si, drop = FALSE],
          p_0[, si, , drop = FALSE],
          p1_0[, si, , drop = FALSE]
        )
        S <- rowSums(l1_c)^2 / info_c + rowSums(l1_s)^2 / info_s
        stat[, "S_S"] <- sign(xi_c[, "theta"] - xi_s[, "theta"]) * sqrt(S)
      }
      if ("W_S" %in% method) {
        stat[, "W_S"] <- (xi_c[, "theta"] - xi_s[, "theta"]) /
          sqrt(1 / info_c + 1 / info_s)
      }
    }
  }

  # Compute response time-based statistics
  if (any(c("L_T", "W_T", "L_ST") %in% method)) {
    Lambda_T <- xi_c[, "tau"]^2 * sum(psi[ci, "alpha"]^2) +
      xi_s[, "tau"]^2 * sum(psi[si, "alpha"]^2) -
      xi[, "tau"]^2 * sum(psi[, "alpha"]^2)
    Lambda_T[Lambda_T < 0] <- 0
    L_T <- sign(xi_c[, "tau"] - xi_s[, "tau"]) * sqrt(Lambda_T)
    if ("L_T" %in% method) {
      stat[, "L_T"] <- L_T
    }
    if ("W_T" %in% method) {
      stat[, "W_T"] <- L_T
    }
  }

  # Compute score and response time-based statistics
  if ("L_ST" %in% method) {
    stat[, "L_ST"] <- pmax(L_S, 0)^2 + pmax(L_T, 0)^2
  }

  # Compute p-values
  pval[, colnames(pval) %in% c("L_S", "ML_S", "S_S", "W_S", "L_T", "W_T")] <-
    pnorm(
      stat[, colnames(stat) %in% c("L_S", "ML_S", "S_S", "W_S", "L_T", "W_T")],
      lower.tail = FALSE
    )
  pval[, colnames(pval) == "LR_S"] <-
    1 - stat[, colnames(stat) == "LR_S"]
  pval[, colnames(pval) == "L_ST"] <-
    0.25 * pchisq(stat[, colnames(stat) == "L_ST"], df = 2,
                  lower.tail = FALSE) +
    0.50 * pchisq(stat[, colnames(stat) == "L_ST"], df = 1,
                  lower.tail = FALSE) +
    0.25 * pchisq(stat[, colnames(stat) == "L_ST"], df = 0,
                  lower.tail = FALSE)

  # Compute flagging rates
  for (a in 1:length(alpha)) {
    flag[, , a] <- pval <= alpha[a]
  }

  # Output
  list(stat = stat, pval = pval, flag = flag)
}
