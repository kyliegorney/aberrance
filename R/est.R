#' Estimate person parameters
#'
#' @noRd

est <- function(interval, psi, x = NULL, r = NULL, y = NULL) {
  N <- max(nrow(x), nrow(r), nrow(y))
  xi_names <- c("theta", "eta", "tau")[c(!is.null(x), !is.null(r), !is.null(y))]
  xi <- matrix(
    nrow = N,
    ncol = length(xi_names),
    dimnames = list(NULL, xi_names)
  )
  if (!is.null(x)) {
    m <- count(psi, ignore = "lambda1")
    xi[, "theta"] <- sapply(1:N, function(v)
      optimize(est_theta, interval, psi = psi, x = x[v, ], m = m,
               maximum = TRUE)$maximum)
  }
  if (!is.null(r)) {
    m <- count(psi)
    xi[, "eta"] <- sapply(1:N, function(v)
      optimize(est_eta, interval, psi = psi, r = r[v, ], m = m,
               maximum = TRUE)$maximum)
  }
  if (!is.null(y)) {
    xi[, "tau"] <- sapply(1:N, function(v)
      optimize(est_tau, interval, psi = psi, y = y[v, ],
               maximum = TRUE)$maximum)
  }
  xi
}

#' Estimate theta
#'
#' @noRd

est_theta <- function(theta, psi, x, m) {
  p <- est_p(m, psi, theta = theta, ignore = "lambda1")
  sum(outer(x, 0:(max(m)-1), "==") * log(p), na.rm = TRUE)
}

#' Estimate eta
#'
#' @noRd

est_eta <- function(eta, psi, r, m) {
  p <- est_p(m, psi, eta = eta, ignore = "b")
  sum(outer(r, 1:max(m), "==") * log(p), na.rm = TRUE)
}

#' Estimate tau
#'
#' @noRd

est_tau <- function(tau, psi, y) {
  sum(-(psi[, "alpha"] * (y - psi[, "beta"] + tau))^2, na.rm = TRUE)
}

#' Estimate IRT probabilities
#'
#' @noRd

est_p <- function(m, psi, theta = NULL, eta = NULL, ignore = NULL) {
  psi <- psi[, !colnames(psi) %in% ignore, drop = FALSE]
  n <- nrow(psi)
  p <- matrix(nrow = n, ncol = max(m))
  if ("b" %in% colnames(psi)) { # 3PL model
    p[, 2] <- psi[, "c"] + (1 - psi[, "c"]) /
      (1 + exp(psi[, "a"] * (psi[, "b"] - theta)))
    p[, 1] <- 1 - p[, 2]
  } else if ("b1" %in% colnames(psi)) { # graded response model
    for (i in 1:n) {
      tmp <- c(
        1,
        1 / (1 + exp(psi[i, "a"] * (psi[i, paste0("b", 1:(m[i]-1))] - theta))),
        0
      )
      p[i, 1:m[i]] <- tmp[1:m[i]] - tmp[2:(m[i]+1)]
    }
  } else if ("c0" %in% colnames(psi)) { # generalized partial credit model
    for (i in 1:n) {
      p[i, 1:m[i]] <- proportions(exp(cumsum(
        -psi[i, "a"] * (psi[i, paste0("c", 0:(m[i]-1))] - theta))))
    }
  } else { # nominal response model
    for (i in 1:n) {
      p[i, 1:m[i]] <- proportions(exp(
        psi[i, paste0("lambda", 1:m[i])] * eta +
          psi[i, paste0("zeta", 1:m[i])]))
    }
  }
  p
}
