#' Detect change point
#' 
#' @description Detect a single change point.
#' 
#' @param method The change point analysis statistic(s) to compute. Options for
#'   score-based statistics are:
#'   - `"L_S_*"` for the likelihood ratio test-based statistic (Shao et al.,
#'     2016; Sinharay, 2016; Tu et al., 2023).
#'   - `"S_S_*"` for the score test-based statistic (Sinharay, 2016; Tu et al.,
#'     2023).
#'   - `"W_S_*"` for the Wald test-based statistic (Sinharay, 2016; Tu et al.,
#'     2023).
#'   
#'   Options for response time-based statistics are:
#'   - `"L_T_*"` for the likelihood ratio test-based statistic (Cheng & Shao,
#'     2022).
#'   - `"W_T_*"` for the Wald test-based statistic (Cheng & Shao, 2022).
#'
#'   Statistics ending in `"*"` can be computed in different ways. Options are:
#'   - `"*"` for all possible ways.
#'   - `"MAX2"` for the maximum of all two-sided statistics in the change point
#'     interval.
#'   - `"MAX1"` for the maximum of all one-sided statistics in the change point
#'     interval.
#'   - `"MIN1"` for the minimum of all one-sided statistics in the change point
#'     interval.
#' 
#' @param cpi The interval to search for the change point. The lower endpoint
#'   must be greater than or equal to 1 and the upper endpoint must be less than
#'   the number of items in the test.
#' 
#' @param xi_c,xi_s Arrays of person parameters. The first dimension corresponds
#'   to persons, the second dimension to parameters, and the third dimension to
#'   change point locations. `xi_c` is based on the items before or at the
#'   change point and `xi_s` is based on the items after the change point. If
#'   `NULL` (default), person parameters are estimated using maximum likelihood
#'   estimation.
#'
#' @param x,y Matrices of raw data. Rows correspond to persons and columns to
#'   items. `x` is for the item scores and `y` the item log response times.
#'
#' @inheritParams detect_pm
#'
#' @returns A list is returned with the following elements:
#' \item{stat}{A matrix of change point analysis statistics.}
#' \item{cp}{A matrix of estimated change points.}
#'
#' @references
#' Cheng, Y., & Shao, C. (2022). Application of change point analysis of
#' response time data to detect test speededness. *Educational and Psychological
#' Measurement*, *82*(5), 1031--1062.
#' 
#' Shao, C., Li, J., & Cheng, Y. (2016). Detection of test speededness using
#' change-point analysis. *Psychometrika*, *81*(4), 1118--1141.
#' 
#' Sinharay, S. (2016). Person fit analysis in computerized adaptive testing
#' using tests for a change point. *Journal of Educational and Behavioral
#' Statistics*, *41*(5), 521--549.
#' 
#' Tu, D., Li, Y., & Cai, Y. (2023). A new perspective on detecting performance
#' decline: A change-point analysis based on Jensen-Shannon divergence.
#' *Behavior Research Methods*, *55*(3), 963--980.
#' 
#' @examples
#' # Setup for Examples 1 and 2 ------------------------------------------------
#'
#' # Settings
#' set.seed(0)     # seed for reproducibility
#' N <- 50         # number of persons
#' n <- 40         # number of items
#' 
#' # Randomly select 10% speeded examinees
#' cv <- sample(1:N, size = N * 0.10)
#' 
#' # Assign change point corresponding to 10% speeded items
#' cp <- n * 0.90
#' ci <- (cp + 1):n
#' 
#' # Create vector of indicators (1 = speeded, 0 = non-speeded)
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
#' # Modify contaminated data by changing the item scores and response times
#' x[cv, ci] <- rbinom(length(cv) * length(ci), size = 1, prob = 0.25)
#' y[cv, ci] <- runif(length(cv) * length(ci), min = log(1), max = log(10))
#' 
#' # Detect change point
#' out <- detect_cp(
#'   method = c("L_S_MAX1", "L_T_MAX1"),
#'   cpi = c(1, n - 1),
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
#' # Modify contaminated data by changing the item scores
#' x[cv, ci] <- rbinom(length(cv) * length(ci), size = 1, prob = 0.25)
#' 
#' # Detect change point in second half of the test
#' out <- detect_cp(
#'   method = c("L_S_MAX1", "S_S_MAX1", "W_S_MAX1"),
#'   cpi = c(20, n - 1),
#'   psi = psi,
#'   x = x
#' )
#' @export

detect_cp <- function(method,
                      cpi,
                      psi,
                      xi = NULL,
                      xi_c = NULL,
                      xi_s = NULL,
                      x = NULL,
                      y = NULL,
                      interval = c(-4, 4)) {
  
  # Checks
  if (cpi[1] > cpi[2]) {
    stop("The lower endpoint of the change point interval must be less than ",
         "or equal to the upper endpoint.", call. = FALSE)
  } else if (cpi[1] < 1) {
    stop("The lower endpoint of the change point interval must be greater ",
         "than or equal to 1.", call. = FALSE)
  } else if (cpi[2] >= nrow(psi)) {
    stop("The upper endpoint of the change point interval must be less than ",
         "the number of items in the test.", call. = FALSE)
  }
  if (any(extract(method, 3) == "*")) {
    tmp <- NULL
    for (m in method) {
      if (extract(m, 3) == "*") {
        tmp <- c(
          tmp,
          paste(
            extract(m, 1:2),
            c("MAX2", "MAX1", "MIN1"),
            sep = "_"
          )
        )
      } else {
        tmp <- c(tmp, m)
      }
    }
    method <- tmp
  }
  if ("S" %in% extract(method, 2)) {
    check_par("x", psi, xi)
  }
  if ("T" %in% extract(method, 2)) {
    check_par("y", psi, xi)
  }
  method <- match.arg(
    arg = unique(method),
    choices = t(outer(
      c("L_S", "S_S", "W_S", "L_T", "W_T"),
      c("MAX2", "MAX1", "MIN1"),
      paste, sep = "_"
    )),
    several.ok = TRUE
  )
  check_data(x, y)
  
  # Setup
  N <- max(nrow(x), nrow(y))
  n <- max(ncol(x), ncol(y))
  stat <- cp <- matrix(
    nrow = N,
    ncol = length(method),
    dimnames = list(
      person = 1:N,
      method = method
    )
  )
  
  # Estimate person parameters
  if (is.null(xi)) {
    xi <- est(interval, psi, x = x, y = y)
  }
  if (is.null(xi_c)) {
    xi_c <- array(
      dim = c(N, ncol(xi), n - 1),
      dimnames = list(NULL, colnames(xi), NULL)
    )
    for (i in cpi[1]:cpi[2]) {
      xi_c[, , i] <- est(
        interval,
        psi[1:i, , drop = FALSE],
        x = x[, 1:i, drop = FALSE],
        y = y[, 1:i, drop = FALSE]
      )
    }
  }
  if (is.null(xi_s)) {
    xi_s <- array(
      dim = c(N, ncol(xi), n - 1),
      dimnames = list(NULL, colnames(xi), NULL)
    )
    for (i in cpi[1]:cpi[2]) {
      xi_s[, , i] <- est(
        interval,
        psi[(i+1):n, , drop = FALSE],
        x = x[, (i+1):n, drop = FALSE],
        y = y[, (i+1):n, drop = FALSE]
      )
    }
  }
  
  # Compute score-based statistics
  if ("S" %in% extract(method, 2)) {
    m <- count(psi, ignore = "lambda1")
    p_0 <- irt_p(m, psi, xi, ignore = "lambda1")
    p1_0 <- irt_p1(p_0, m, psi, xi, ignore = "lambda1")
    if ("L_S" %in% extract(method, 1:2)) {
      p_1 <- p_0
      L_S <- matrix(nrow = N, ncol = n - 1)
      for (i in cpi[1]:cpi[2]) {
        p_1[, 1:i, 1:max(m[1:i])] <- irt_p(
          m[1:i],
          psi[1:i, , drop = FALSE],
          cbind(theta = xi_c[, "theta", i]),
          ignore = "lambda1"
        )
        p_1[, (i+1):n, 1:max(m[(i+1):n])] <- irt_p(
          m[(i+1):n],
          psi[(i+1):n, , drop = FALSE],
          cbind(theta = xi_s[, "theta", i]),
          ignore = "lambda1"
        )
        L_S[, i] <- compute_L_S(
          x,
          p_0,
          p_1,
          xi_c[, "theta", i],
          xi_s[, "theta", i],
          signed = FALSE
        )
      }
      if ("L_S_MAX2" %in% method) {
        stat[, "L_S_MAX2"] <- apply(L_S, 1, max, na.rm = TRUE)
        cp[, "L_S_MAX2"] <- apply(L_S, 1, which.max)
      }
      if (any(c("L_S_MAX1", "L_S_MIN1") %in% method)) {
        L_S[L_S < 0] <- 0
        L_S <- sign(xi_c[, "theta", ] - xi_s[, "theta", ]) * sqrt(L_S)
        if ("L_S_MAX1" %in% method) {
          stat[, "L_S_MAX1"] <- apply(L_S, 1, max, na.rm = TRUE)
          cp[, "L_S_MAX1"] <- apply(L_S, 1, which.max)
        }
        if ("L_S_MIN1" %in% method) {
          stat[, "L_S_MIN1"] <- apply(L_S, 1, min, na.rm = TRUE)
          cp[, "L_S_MIN1"] <- apply(L_S, 1, which.min)
        }
      }
    }
    if ("S_S" %in% extract(method, 1:2)) {
      S_S <- matrix(nrow = N, ncol = n - 1)
      for (i in cpi[1]:cpi[2]) {
        S_S[, i] <- compute_S_S(
          1:i,
          (i+1):n,
          x,
          p_0,
          p1_0,
          xi_c[, "theta", i],
          xi_s[, "theta", i],
          signed = FALSE
        )
      }
      if ("S_S_MAX2" %in% method) {
        stat[, "S_S_MAX2"] <- apply(S_S, 1, max, na.rm = TRUE)
        cp[, "S_S_MAX2"] <- apply(S_S, 1, which.max)
      }
      if (any(c("S_S_MAX1", "S_S_MIN1") %in% method)) {
        S_S[S_S < 0] <- 0
        S_S <- sign(xi_c[, "theta", ] - xi_s[, "theta", ]) * sqrt(S_S)
        if ("S_S_MAX1" %in% method) {
          stat[, "S_S_MAX1"] <- apply(S_S, 1, max, na.rm = TRUE)
          cp[, "S_S_MAX1"] <- apply(S_S, 1, which.max)
        }
        if ("S_S_MIN1" %in% method) {
          stat[, "S_S_MIN1"] <- apply(S_S, 1, min, na.rm = TRUE)
          cp[, "S_S_MIN1"] <- apply(S_S, 1, which.min)
        }
      }
    }
    if ("W_S" %in% extract(method, 1:2)) {
      W_S <- matrix(nrow = N, ncol = n - 1)
      for (i in cpi[1]:cpi[2]) {
        W_S[, i] <- compute_W_S(
          1:i,
          (i+1):n,
          p_0,
          p1_0,
          xi_c[, "theta", i],
          xi_s[, "theta", i],
          signed = FALSE
        )
      }
      if ("W_S_MAX2" %in% method) {
        stat[, "W_S_MAX2"] <- apply(W_S, 1, max, na.rm = TRUE)
        cp[, "W_S_MAX2"] <- apply(W_S, 1, which.max)
      }
      if (any(c("W_S_MAX1", "W_S_MIN1") %in% method)) {
        W_S[W_S < 0] <- 0
        W_S <- sign(xi_c[, "theta", ] - xi_s[, "theta", ]) * sqrt(W_S)
        if ("W_S_MAX1" %in% method) {
          stat[, "W_S_MAX1"] <- apply(W_S, 1, max, na.rm = TRUE)
          cp[, "W_S_MAX1"] <- apply(W_S, 1, which.max)
        }
        if ("W_S_MIN1" %in% method) {
          stat[, "W_S_MIN1"] <- apply(W_S, 1, min, na.rm = TRUE)
          cp[, "W_S_MIN1"] <- apply(W_S, 1, which.min)
        }
      }
    }
  }
  
  # Compute response time-based statistics
  if ("L_T" %in% extract(method, 1:2)) {
    L_T <- matrix(nrow = N, ncol = n - 1)
    for (i in cpi[1]:cpi[2]) {
      L_T[, i] <- compute_L_T(
        1:i,
        (i+1):n,
        y,
        psi,
        xi[, "tau"],
        xi_c[, "tau", i],
        xi_s[, "tau", i],
        signed = FALSE
      )
    }
    if ("L_T_MAX2" %in% method) {
      stat[, "L_T_MAX2"] <- apply(L_T, 1, max, na.rm = TRUE)
      cp[, "L_T_MAX2"] <- apply(L_T, 1, which.max)
    }
    if (any(c("L_T_MAX1", "L_T_MIN1") %in% method)) {
      L_T[L_T < 0] <- 0
      L_T <- sign(xi_c[, "tau", ] - xi_s[, "tau", ]) * sqrt(L_T)
      if ("L_T_MAX1" %in% method) {
        stat[, "L_T_MAX1"] <- apply(L_T, 1, max, na.rm = TRUE)
        cp[, "L_T_MAX1"] <- apply(L_T, 1, which.max)
      }
      if ("L_T_MIN1" %in% method) {
        stat[, "L_T_MIN1"] <- apply(L_T, 1, min, na.rm = TRUE)
        cp[, "L_T_MIN1"] <- apply(L_T, 1, which.min)
      }
    }
  }
  if ("W_T" %in% extract(method, 1:2)) {
    W_T <- matrix(nrow = N, ncol = n - 1)
    for (i in cpi[1]:cpi[2]) {
      W_T[, i] <- compute_W_T(
        1:i,
        (i+1):n,
        psi,
        xi_c[, "tau", i],
        xi_s[, "tau", i],
        signed = FALSE
      )
    }
    if ("W_T_MAX2" %in% method) {
      stat[, "W_T_MAX2"] <- apply(W_T, 1, max, na.rm = TRUE)
      cp[, "W_T_MAX2"] <- apply(W_T, 1, which.max)
    }
    if (any(c("W_T_MAX1", "W_T_MIN1") %in% method)) {
      W_T[W_T < 0] <- 0
      W_T <- sign(xi_c[, "tau", ] - xi_s[, "tau", ]) * sqrt(W_T)
      if ("W_T_MAX1" %in% method) {
        stat[, "W_T_MAX1"] <- apply(W_T, 1, max, na.rm = TRUE)
        cp[, "W_T_MAX1"] <- apply(W_T, 1, which.max)
      }
      if ("W_T_MIN1" %in% method) {
        stat[, "W_T_MIN1"] <- apply(W_T, 1, min, na.rm = TRUE)
        cp[, "W_T_MIN1"] <- apply(W_T, 1, which.min)
      }
    }
  }
  
  # Output
  list(stat = stat, cp = cp)
}
