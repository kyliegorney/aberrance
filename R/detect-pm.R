#' Detect parametric misfit
#'
#' @description Detect parametric misfit using person-fit statistics.
#'
#' @param method The person-fit statistic(s) to compute. Options for score-based
#'   statistics are:
#'   - `"ECI2_S_*"` for the second standardized extended caution index, also
#'     known as \eqn{\zeta_1} (Tatsuoka, 1984; see also Sinharay, 2018b).
#'   - `"ECI4_S_*"` for the fourth standardized extended caution index, also
#'     known as \eqn{\zeta_2} (Tatsuoka, 1984; see also Sinharay, 2018b).
#'   - `"L_S_*"` for the standardized log-likelihood statistic (Drasgow et al.,
#'      1985).
#'
#'   Options for distractor-based statistics are:
#'   - `"L_D_*"` for the standardized log-likelihood statistic (Gorney &
#'     Wollack, 2023).
#'
#'   Options for score and distractor-based statistics are:
#'   - `"L_SD_*"` for the log-likelihood statistic (Gorney & Wollack, 2023).
#'
#'   Options for response-based statistics are:
#'   - `"L_R_*"` for the standardized log-likelihood statistic (Drasgow et al.,
#'      1985).
#'
#'   Options for response time-based statistics are:
#'   - `"L_T"` for the log-likelihood statistic (Sinharay, 2018a).
#'
#'   Options for score and response time-based statistics are:
#'   - `"Q_ST_*"` for the log-likelihood statistic that combines `"L_S_*"` and
#'     `"L_T"` (Gorney, Sinharay, & Liu, 2024).
#'   - `"L_ST_*"` for the standardized log-likelihood statistic (Gorney,
#'     Sinharay, & Liu, 2024). *Note:* This statistic cannot be computed using
#'     the `"CF"`, `"CS"`, `"EW"`, `"TSCF"`, `"TSCS"`, or `"TSEW"` corrections.
#'
#'   Options for response and response time-based statistics are:
#'   - `"Q_RT_*"` for the log-likelihood statistic that combines `"L_R_*"` and
#'     `"L_T"` (Gorney, Sinharay, & Liu, 2024).
#'   - `"L_RT_*"` for the standardized log-likelihood statistic (Gorney,
#'     Sinharay, & Liu, 2024). *Note:* This statistic cannot be computed using
#'     the `"CF"`, `"CS"`, `"EW"`, `"TSCF"`, `"TSCS"`, or `"TSEW"` corrections.
#'
#'   Statistics ending in `"*"` can be computed using various corrections.
#'   Options are:
#'   - `"*"` for all possible corrections.
#'   - `"NO"` for no correction.
#'   - `"CF"` for the Cornish-Fisher expansion (Molenaar & Hoijtink, 1990).
#'   - `"CS"` for the chi-squared approximation (Molenaar & Hoijtink, 1990).
#'   - `"EW"` for the Edgeworth expansion (Bedrick, 1997).
#'   - `"TS"` for the Taylor series expansion (Snijders, 2001; see also
#'     Sinharay, 2016a, 2016b).
#'   - `"TSCF"` for the Taylor series expansion and Cornish-Fisher expansion
#'     (Gorney, Sinharay, & Eckerly, 2024; see also Gorney, 2024).
#'   - `"TSCS"` for the Taylor series expansion and chi-squared approximation
#'     (Gorney Sinharay, & Eckerly, 2024; see also Gorney, 2024).
#'   - `"TSEW"` for the Taylor series expansion and Edgeworth expansion (Gorney
#'     Sinharay, & Eckerly, 2024; see also Gorney, 2024).
#'
#' @param psi A matrix of item parameters.
#'
#' @param xi A matrix of person parameters. If `NULL` (default), person
#'   parameters are estimated using maximum likelihood estimation.
#'
#' @param x,d,r,y Matrices of raw data. `x` is for the item scores, `d` the item
#'    distractors, `r` the item responses, and `y` the item log response times.
#'
#' @param interval The interval to search for the person parameters. Default is
#'   `c(-4, 4)`.
#'
#' @param alpha Value(s) between 0 and 1 indicating the significance level(s)
#'   used for flagging. Default is `0.05`.
#'
#' @returns A list is returned with the following elements:
#' \item{stat}{A matrix of parametric person-fit statistics.}
#' \item{pval}{A matrix of *p*-values.}
#' \item{flag}{An array of flagging results. The first dimension corresponds to
#'   persons, the second dimension to methods, and the third dimension to
#'   significance levels.}
#'
#' @references
#' Bedrick, E. J. (1997). Approximating the conditional distribution of person
#' fit indexes for checking the Rasch model. *Psychometrika*, *62*(2), 191--199.
#'
#' Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness
#' measurement with polychotomous item response models and standardized indices.
#' *British Journal of Mathematical and Statistical Psychology*, *38*(1),
#' 67--86.
#'
#' Gorney, K. (2024). Three new corrections for standardized person-fit
#' statistics for tests with polytomous items. *British Journal of Mathematical
#' and Statistical Psychology*. Advance online publication.
#'
#' Gorney, K., Sinharay, S., & Eckerly, C. (2024). Efficient corrections for
#' standardized person-fit statistics. *Psychometrika*, *89*(2), 569--591.
#'
#' Gorney, K., Sinharay, S., & Liu, X. (2024). Using item scores and response
#' times in person-fit assessment. *British Journal of Mathematical and
#' Statistical Psychology*, *77*(1), 151--168.
#'
#' Gorney, K., & Wollack, J. A. (2023). Using item scores and distractors in
#' person-fit assessment. *Journal of Educational Measurement*, *60*(1), 3--27.
#'
#' Molenaar, I. W., & Hoijtink, H. (1990). The many null distributions of person
#' fit indices. *Psychometrika*, *55*(1), 75--106.
#'
#' Sinharay, S. (2016a). Asymptotic corrections of standardized extended caution
#' indices. *Applied Psychological Measurement*, *40*(6), 418--433.
#'
#' Sinharay, S. (2016b). Asymptotically correct standardization of person-fit
#' statistics beyond dichotomous items. *Psychometrika*, *81*(4), 992--1013.
#'
#' Sinharay, S. (2018a). A new person-fit statistic for the lognormal model for
#' response times. *Journal of Educational Measurement*, *55*(4), 457--476.
#'
#' Sinharay, S. (2018b). Extension of caution indices to mixed-format tests.
#' *British Journal of Mathematical and Statistical Psychology*, *71*(2),
#' 363--386.
#'
#' Snijders, T. A. B. (2001). Asymptotic null distribution of person fit
#' statistics with estimated person parameter. *Psychometrika*, *66*(3),
#' 331--342.
#'
#' Tatsuoka, K. K. (1984). Caution indices based on item response theory.
#' *Psychometrika*, *49*(1), 95--110.
#'
#' @seealso [detect_nm()] to detect nonparametric misfit.
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
#' # Create vector of indicators (1 = misfitting, 0 = fitting)
#' ind <- ifelse(1:N %in% cv, 1, 0)
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
#' # Detect parametric misfit
#' out <- detect_pm(
#'   method = c("L_S_TS", "L_T", "Q_ST_TS", "L_ST_TS"),
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
#' # Detect parametric misfit
#' out <- detect_pm(
#'   method = c("ECI2_S_TSCF", "ECI4_S_TSCF", "L_S_TSCF"),
#'   psi = psi,
#'   x = x
#' )
#' @export

detect_pm <- function(method,
                      psi,
                      xi = NULL,
                      x = NULL, d = NULL, r = NULL, y = NULL,
                      interval = c(-4, 4),
                      alpha = 0.05) {

  # Checks
  if (any(extract(method, 3) == "*")) {
    tmp <- NULL
    for (m in method) {
      if (extract(m, 3) == "*") {
        if (extract(m, 1:2) %in% c("L_ST", "L_RT")) {
          tmp <- c(tmp, paste(
            extract(m, 1:2),
            c("NO", "TS"),
            sep = "_")
          )
        } else {
          tmp <- c(tmp, paste(
            extract(m, 1:2),
            c("NO", "CF", "CS", "EW", "TS", "TSCF", "TSCS", "TSEW"),
            sep = "_"
          ))
        }
      } else {
        tmp <- c(tmp, m)
      }
    }
    method <- tmp
  }
  if (any(c("S", "D", "SD", "ST") %in% extract(method, 2)) &&
      any(c("R", "RT") %in% extract(method, 2))) {
    stop("`method` may contain either score and distractor-based statistics ",
         "or response-based statistics, but not both.", call. = FALSE)
  }
  if (any(c("S", "SD", "ST") %in% extract(method, 2))) {
    check_par("x", psi, xi)
  }
  if (any(c("D", "SD") %in% extract(method, 2))) {
    check_par("d", psi, xi)
  }
  if (any(c("R", "RT") %in% extract(method, 2))) {
    check_par("r", psi, xi)
  }
  if (any(c("T", "ST", "RT") %in% extract(method, 2))) {
    check_par("y", psi, xi)
    if (any(paste("L_ST", c("CF", "CS", "EW", "TSCF", "TSCS", "TSEW"),
                  sep = "_") %in% method)) {
      warning("The L_ST statistic cannot be computed using the CF, CS, EW, ",
              "TSCF, TSCS, or TSEW corrections.", call. = FALSE)
    }
    if (any(paste("L_RT", c("CF", "CS", "EW", "TSCF", "TSCS", "TSEW"),
                  sep = "_") %in% method)) {
      warning("The L_RT statistic cannot be computed using the CF, CS, EW, ",
              "TSCF, TSCS, or TSEW corrections.", call. = FALSE)
    }
  }
  mdc <- match.arg(
    arg = unique(method),
    choices = c(
      t(outer(
        c("ECI2_S", "ECI4_S", "L_S", "L_D", "L_SD", "L_R", "Q_ST", "Q_RT"),
        c("NO", "CF", "CS", "EW", "TS", "TSCF", "TSCS", "TSEW"),
        paste, sep = "_"
      )),
      "L_T",
      t(outer(
        c("L_ST", "L_RT"),
        c("NO", "TS"),
        paste, sep = "_"
      ))
    ),
    several.ok = TRUE
  )
  check_data(x, d, r, y)

  # Setup
  N <- max(nrow(x), nrow(d), nrow(r), nrow(y))
  n <- max(ncol(x), ncol(d), ncol(r), ncol(y))
  stat <- pval <- matrix(
    nrow = N, ncol = length(mdc),
    dimnames = list(
      person = 1:N,
      method = mdc
    )
  )
  flag <- array(
    dim = c(N, length(mdc), length(alpha)),
    dimnames = list(
      person = 1:N,
      method = mdc,
      alpha = alpha
    )
  )
  md <- extract(mdc, 1:2)
  c <- extract(mdc, 3)

  # Estimate person parameters
  if (is.null(xi)) {
    xi <- est(interval, psi, x = x, d = d, r = r, y = y)
  }

  # Compute score-based statistics
  if (any(c("ECI2_S", "ECI4_S", "L_S", "L_SD", "Q_ST") %in% md)) {
    m <- count(psi, ignore = "lambda1")
    p <- irt_p(m, psi, xi, ignore = "lambda1")
    p1 <- irt_p1(p, m, psi, xi, ignore = "lambda1")
    if (any(c("ECI2_S", "ECI4_S") %in% md)) {
      tmp <- mdc[md %in% c("ECI2_S", "ECI4_S")]
      stat[, tmp] <- compute_SPF_S(tmp, x, p, p1)
    }
    if (any(c("L_S", "L_SD", "Q_ST") %in% md)) {
      tmp <- unique(c(
        mdc[md == "L_S"],
        paste("L_S", c[md %in% c("L_SD", "Q_ST")], sep = "_")
      ))
      spf <- compute_SPF_S(tmp, x, p, p1)
      if ("L_S" %in% md) {
        tmp <- mdc[md == "L_S"]
        stat[, tmp] <- spf[, tmp]
      }
      if ("L_SD" %in% md) {
        tmp <- paste("L_S", c[md == "L_SD"], sep = "_")
        L_S_in_L_SD <- spf[, tmp]
      }
      if ("Q_ST" %in% md) {
        tmp <- paste("L_S", c[md == "Q_ST"], sep = "_")
        L_S_in_Q_ST <- spf[, tmp]
      }
    }
  }

  # Compute distractor-based statistics
  if (any(c("L_D", "L_SD") %in% md)) {
    m <- count(psi, ignore = "b")
    p <- irt_p(m, psi, xi, ignore = "b")
    p1 <- irt_p1(p, m, psi, xi, ignore = "b")
    tmp <- unique(c(
      mdc[md == "L_D"],
      paste("L_D", c[md == "L_SD"], sep = "_")
    ))
    spf <- compute_SPF_S(tmp, d - 1, p, p1)
    if ("L_D" %in% md) {
      tmp <- mdc[md == "L_D"]
      stat[, tmp] <- spf[, tmp]
    }
    if ("L_SD" %in% md) {
      tmp <- paste("L_D", c[md == "L_SD"], sep = "_")
      L_D_in_L_SD <- spf[, tmp]
    }
  }

  # Compute score and distractor-based statistics
  if ("L_SD" %in% md) {
    tmp <- mdc[md == "L_SD"]
    stat[, tmp] <- pmin(L_S_in_L_SD, 0)^2 + pmin(L_D_in_L_SD, 0)^2
  }

  # Compute response-based statistics
  if (any(c("L_R", "Q_RT") %in% md)) {
    m <- count(psi)
    p <- irt_p(m, psi, xi)
    p1 <- irt_p1(p, m, psi, xi)
    tmp <- unique(c(
      mdc[md == "L_R"],
      paste("L_R", c[md == "Q_RT"], sep = "_")
    ))
    spf <- compute_SPF_S(tmp, r - 1, p, p1)
    if ("L_R" %in% md) {
      tmp <- mdc[md == "L_R"]
      stat[, tmp] <- spf[, tmp]
    }
    if ("Q_RT" %in% md) {
      tmp <- paste("L_R", c[md == "Q_RT"], sep = "_")
      L_R_in_Q_RT <- spf[, tmp]
    }
  }

  # Compute response time-based statistics
  if (any(c("L_T", "Q_ST", "Q_RT") %in% md)) {
    L_T <- rowSums(
      t(psi[, "alpha"] * (t(y) - outer(psi[, "beta"], xi[, "tau"], "-")))^2)
    if ("L_T" %in% md) {
      stat[, "L_T"] <- L_T
    }
  }

  # Compute score and response time-based statistics
  if ("Q_ST" %in% md) {
    tmp <- mdc[md == "Q_ST"]
    stat[, tmp] <-
      qchisq(pnorm(L_S_in_Q_ST, lower.tail = TRUE),
             df = 1, lower.tail = FALSE) +
      qchisq(pchisq(L_T, df = n - 1, lower.tail = FALSE),
             df = 1, lower.tail = FALSE)
  }
  if ("L_ST" %in% md) {
    m <- count(psi, ignore = "lambda1")
    p <- irt_p(m, psi, xi, ignore = "lambda1")
    p1 <- irt_p1(p, m, psi, xi, ignore = "lambda1")
    m <- t(outer(psi[, "beta"], xi[, "tau"], "-"))
    s <- matrix(1 / psi[, "alpha"], nrow = N, ncol = n, byrow = TRUE)
    tmp <- mdc[md == "L_ST"]
    stat[, tmp] <- compute_SPF_ST(tmp, x, y, p, p1, m, s)
  }

  # Compute response and response time-based statistics
  if ("Q_RT" %in% md) {
    tmp <- mdc[md == "Q_RT"]
    stat[, tmp] <-
      qchisq(pnorm(L_R_in_Q_RT, lower.tail = TRUE),
             df = 1, lower.tail = FALSE) +
      qchisq(pchisq(L_T, df = n - 1, lower.tail = FALSE),
             df = 1, lower.tail = FALSE)
  }
  if ("L_RT" %in% md) {
    m <- count(psi)
    p <- irt_p(m, psi, xi)
    p1 <- irt_p1(p, m, psi, xi)
    m <- t(outer(psi[, "beta"], xi[, "tau"], "-"))
    s <- matrix(1 / psi[, "alpha"], nrow = N, ncol = n, byrow = TRUE)
    tmp <- mdc[md == "L_RT"]
    stat[, tmp] <- compute_SPF_ST(tmp, r - 1, y, p, p1, m, s)
  }

  # Compute p-values
  pval[, md %in% c("ECI2_S", "ECI4_S")] <-
    pnorm(stat[, md %in% c("ECI2_S", "ECI4_S")], lower.tail = FALSE)
  pval[, md %in% c("L_S", "L_D", "L_R", "L_ST", "L_RT")] <-
    pnorm(stat[, md %in% c("L_S", "L_D", "L_R", "L_ST", "L_RT")],
          lower.tail = TRUE)
  pval[, md == "L_SD"] <-
    0.25 * pchisq(stat[, md == "L_SD"], df = 2, lower.tail = FALSE) +
    0.50 * pchisq(stat[, md == "L_SD"], df = 1, lower.tail = FALSE) +
    0.25 * pchisq(stat[, md == "L_SD"], df = 0, lower.tail = FALSE)
  pval[, md == "L_T"] <-
    pchisq(stat[, md == "L_T"], df = n - 1, lower.tail = FALSE)
  pval[, md %in% c("Q_ST", "Q_RT")] <-
    pchisq(stat[, md %in% c("Q_ST", "Q_RT")], df = 2, lower.tail = FALSE)

  # Compute flagging rates
  for (a in 1:length(alpha)) {
    flag[, , a] <- pval <= alpha[a]
  }

  # Output
  list(stat = stat, pval = pval, flag = flag)
}
