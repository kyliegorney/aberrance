#' Compute the erasure detection index with continuity correction
#'
#' @noRd

compute_EDI_CO <- function(s, p, c) {
  (sum(rowSums(s - p, na.rm = TRUE)) + c) /
    sqrt(sum(rowSums(p * (1 - p), na.rm = TRUE)))
}

#' Compute the erasure detection index with Taylor series expansion
#'
#' @noRd

compute_EDI_TS <- function(s, p, p1, info) {
  sum(rowSums(s - p, na.rm = TRUE)) /
    sqrt(sum(rowSums(p * (1 - p), na.rm = TRUE)) +
           sum(rowSums(p1, na.rm = TRUE)^2 / info))
}

#' Compute the signed likelihood ratio test statistic
#'
#' @noRd

compute_L <- function(x, p_0, p_1, theta_c, theta_s) {
  l_0 <- irt_l(x, p_0)
  l_1 <- irt_l(x, p_1)
  L <- 2 * rowSums(l_1 - l_0)
  L[L < 0] <- 0
  sign(theta_c - theta_s) * sqrt(L)
}

#' Compute the omega statistic
#'
#' @noRd

compute_OMG <- function(s, p, c = 0) {
  (sum(s - p) + c) / sqrt(sum(p * (1 - p)))
}

#' Compute skewness corrections for standardized person-fit statistics
#'
#' @noRd

compute_SKEW <- function(method, correct, W, mu, sigma, gamma) {
  no <- (W - mu) / sigma
  if ((correct == "NO") || (correct == "TS")) {
    stat <- no
  } else if ((correct == "CF") || (correct == "TSCF")) {
    stat <- no - gamma * (no^2 - 1) / 12
  } else if ((correct == "CS") || (correct == "TSCS")) {
    nu <- 8 / gamma^2
    b <- sqrt(sigma^2 / (2 * nu))
    a <- b * nu
    if (method == "L") {
      stat <- pchisq(abs(W - mu - a) / b,
                     df = nu, lower.tail = FALSE, log.p = TRUE)
      stat <- qnorm(stat, lower.tail = TRUE, log.p = TRUE)
    } else {
      stat <- pchisq(abs(W - mu + a) / b,
                     df = nu, lower.tail = FALSE, log.p = TRUE)
      stat <- qnorm(stat, lower.tail = FALSE, log.p = TRUE)
    }
  } else {
    if (method == "L") {
      stat <- pnorm(no) - dnorm(no) * gamma * (no^2 - 1) / 6
      stat <- ifelse(stat <= 0 | stat >= 1,
                     pnorm(no, lower.tail = TRUE, log.p = TRUE), log(stat))
      stat <- qnorm(stat, lower.tail = TRUE, log.p = TRUE)
    } else {
      stat <- 1 - (pnorm(no) - dnorm(no) * gamma * (no^2 - 1) / 6)
      stat <- ifelse(stat <= 0 | stat >= 1,
                     pnorm(no, lower.tail = FALSE, log.p = TRUE), log(stat))
      stat <- qnorm(stat, lower.tail = FALSE, log.p = TRUE)
    }
  }
  stat
}

#' Compute standardized person-fit statistics for item scores
#'
#' @noRd

compute_SPF_S <- function(mdc, x, p, p1) {
  N <- dim(p)[1]
  n <- dim(p)[2]
  max_m <- dim(p)[3]
  stat <- matrix(nrow = N, ncol = length(mdc), dimnames = list(NULL, mdc))
  method <- unique(extract(mdc, 1))
  d <- unique(extract(mdc, 2))
  if (max_m == 2) {
    q <- p[, , 1]
    p <- p[, , 2]
    p1 <- p1[, , 2]
    for (m in method) {
      correct <- extract(mdc[extract(mdc, 1) == m], 3)
      if (m == "ECI2") {
        g_i <- colMeans(p)
        g <- mean(g_i)
        w <- matrix(g - g_i, nrow = N, ncol = n, byrow = TRUE)
      } else if (m == "ECI4") {
        w <- rowMeans(p) - p
      } else {
        w <- log(p / q)
      }
      W <- rowSums((x - p) * w)
      if (any(c("NO", "CF", "CS", "EW") %in% correct)) {
        mu <- 0
        sigma <- sqrt(rowSums(p * q * w^2))
        gamma <- rowSums(p * q * (q - p) * w^3) / sigma^3
        for (c in correct[correct %in% c("NO", "CF", "CS", "EW")]) {
          stat[, paste(m, d, c, sep = "_")] <-
            compute_SKEW(m, c, W, mu, sigma, gamma)
        }
      }
      if (any(c("TS", "TSCF", "TSCS", "TSEW") %in% correct)) {
        r <- p1 / (p * q)
        r0 <- 0
        c <- rowSums(p1 * w) / rowSums(p1 * r)
        w <- w - c * r
        mu <- -c * r0
        sigma <- sqrt(rowSums(p * q * w^2))
        gamma <- rowSums(p * q * (q - p) * w^3) / sigma^3
        for (c in correct[correct %in% c("TS", "TSCF", "TSCS", "TSEW")]) {
          stat[, paste(m, d, c, sep = "_")] <-
            compute_SKEW(m, c, W, mu, sigma, gamma)
        }
      }
    }
  } else {
    I <- outer(x, 0:(max_m-1), "==")
    I[is.na(p)] <- NA
    for (m in method) {
      correct <- extract(mdc[extract(mdc, 1) == m], 3)
      if (m == "ECI2") {
        g_ij <- apply(p, c(2, 3), mean)
        g_j <- colMeans(g_ij, na.rm = TRUE)
        w <- array(rep(t(g_j - t(g_ij)), each = N), dim = c(N, n, max_m))
      } else if (m == "ECI4") {
        w <- aperm(array(apply(p, c(1, 3), mean, na.rm = TRUE),
                         dim = c(N, max_m, n)), c(1, 3, 2)) - p
      } else {
        w <- log(p)
      }
      W <- apply((I - p) * w, 1, sum, na.rm = TRUE)
      if (any(c("NO", "CF", "CS", "EW") %in% correct)) {
        mu <- 0
        sigma <- gamma <- rep(NA, times = N)
        for (v in 1:N) {
          tmp2 <- sum(p[v, , ] * (1 - p[v, , ]) * w[v, , ]^2, na.rm = TRUE)
          tmp3 <- sum(p[v, , ] * (p[v, , ] - 1) * (2 * p[v, , ] - 1) *
                        w[v, , ]^3, na.rm = TRUE)
          for (j in 1:max_m) {
            for (k in (1:max_m)[-j]) {
              tmp2 <- tmp2 -
                sum(p[v, , j] * p[v, , k] * w[v, , j] * w[v, , k], na.rm = TRUE)
              tmp3 <- tmp3 +
                3 * sum(p[v, , j] * p[v, , k] * (2 * p[v, , j] - 1) *
                          w[v, , j]^2 * w[v, , k], na.rm = TRUE)
            }
          }
          for (j in 1:(max_m-2)) {
            for (k in (j+1):(max_m-1)) {
              for (l in (k+1):max_m) {
                tmp3 <- tmp3 +
                  6 * sum(2 * p[v, , j] * p[v, , k] * p[v, , l] *
                            w[v, , j] * w[v, , k] * w[v, , l], na.rm = TRUE)
              }
            }
          }
          sigma[v] <- sqrt(tmp2)
          gamma[v] <- tmp3 / sigma[v]^3
        }
        for (c in correct[correct %in% c("NO", "CF", "CS", "EW")]) {
          stat[, paste(m, d, c, sep = "_")] <-
            compute_SKEW(m, c, W, mu, sigma, gamma)
        }
      }
      if (any(c("TS", "TSCF", "TSCS", "TSEW") %in% correct)) {
        r <- p1 / p
        r0 <- 0
        c <- apply(p1 * w, 1, sum, na.rm = TRUE) /
          apply(p1 * r, 1, sum, na.rm = TRUE)
        w <- w - c * r
        mu <- -c * r0
        sigma <- gamma <- rep(NA, times = N)
        for (v in 1:N) {
          tmp2 <- sum(p[v, , ] * (1 - p[v, , ]) * w[v, , ]^2, na.rm = TRUE)
          tmp3 <- sum(p[v, , ] * (p[v, , ] - 1) * (2 * p[v, , ] - 1) *
                        w[v, , ]^3, na.rm = TRUE)
          for (j in 1:max_m) {
            for (k in (1:max_m)[-j]) {
              tmp2 <- tmp2 -
                sum(p[v, , j] * p[v, , k] * w[v, , j] * w[v, , k], na.rm = TRUE)
              tmp3 <- tmp3 +
                3 * sum(p[v, , j] * p[v, , k] * (2 * p[v, , j] - 1) *
                          w[v, , j]^2 * w[v, , k], na.rm = TRUE)
            }
          }
          for (j in 1:(max_m-2)) {
            for (k in (j+1):(max_m-1)) {
              for (l in (k+1):max_m) {
                tmp3 <- tmp3 +
                  6 * sum(2 * p[v, , j] * p[v, , k] * p[v, , l] *
                            w[v, , j] * w[v, , k] * w[v, , l], na.rm = TRUE)
              }
            }
          }
          sigma[v] <- sqrt(tmp2)
          gamma[v] <- tmp3 / sigma[v]^3
        }
        for (c in correct[correct %in% c("TS", "TSCF", "TSCS", "TSEW")]) {
          stat[, paste(m, d, c, sep = "_")] <-
            compute_SKEW(m, c, W, mu, sigma, gamma)
        }
      }
    }
  }
  stat
}

#' Compute standardized person-fit statistics for item scores and response times
#'
#' @noRd

compute_SPF_ST <- function(mdc, x, y, p, p1, m, s) {
  N <- dim(p)[1]
  n <- dim(p)[2]
  max_m <- dim(p)[3]
  stat <- matrix(nrow = N, ncol = length(mdc), dimnames = list(NULL, mdc))
  d <- unique(extract(mdc, 2))
  correct <- extract(mdc, 3)
  if (max_m == 2) {
    q <- p[, , 1]
    p <- p[, , 2]
    p1 <- p1[, , 2]
    w <- log(p / q)
    W <- rowSums((x - p) * w - 0.5 * ((y - m) / s)^2 + 0.5)
    if ("NO" %in% correct) {
      mu <- 0
      sigma <- sqrt(rowSums(p * q * w^2 + 0.5))
      stat[, paste("L", d, "NO", sep = "_")] <- (W - mu) / sigma
    }
    if ("TS" %in% correct) {
      r <- p1 / (p * q)
      r0 <- 0
      c <- rowSums(p1 * w) / rowSums(p1 * r)
      w <- w - c * r
      mu <- -c * r0
      sigma <- sqrt(rowSums(p * q * w^2 + 0.5))
      stat[, paste("L", d, "TS", sep = "_")] <- (W - mu) / sigma
    }
  } else {
    I <- outer(x, 0:(max_m-1), "==")
    I[is.na(p)] <- NA
    w <- log(p)
    W <- apply((I - p) * w, 1, sum, na.rm = TRUE) +
      rowSums(-0.5 * ((y - m) / s)^2 + 0.5)
    if ("NO" %in% correct) {
      mu <- 0
      sigma <- rep(NA, times = N)
      for (v in 1:N) {
        tmp2 <- sum(p[v, , ] * (1 - p[v, , ]) * w[v, , ]^2, na.rm = TRUE)
        for (j in 1:max_m) {
          for (k in (1:max_m)[-j]) {
            tmp2 <- tmp2 -
              sum(p[v, , j] * p[v, , k] * w[v, , j] * w[v, , k], na.rm = TRUE)
          }
        }
        sigma[v] <- sqrt(tmp2 + 0.5 * n)
      }
      stat[, paste("L", d, "NO", sep = "_")] <- (W - mu) / sigma
    }
    if ("TS" %in% correct) {
      r <- p1 / p
      r0 <- 0
      c <- apply(p1 * w, 1, sum, na.rm = TRUE) /
        apply(p1 * r, 1, sum, na.rm = TRUE)
      w <- w - c * r
      mu <- -c * r0
      sigma <- rep(NA, times = N)
      for (v in 1:N) {
        tmp2 <- sum(p[v, , ] * (1 - p[v, , ]) * w[v, , ]^2, na.rm = TRUE)
        for (j in 1:max_m) {
          for (k in (1:max_m)[-j]) {
            tmp2 <- tmp2 -
              sum(p[v, , j] * p[v, , k] * w[v, , j] * w[v, , k], na.rm = TRUE)
          }
        }
        sigma[v] <- sqrt(tmp2 + 0.5 * n)
      }
      stat[, paste("L", d, "TS", sep = "_")] <- (W - mu) / sigma
    }
  }
  stat
}

#' Compute the weighted omega statistic
#'
#' @noRd

compute_WOMG <- function(s, p) {
  w <- -log(p)
  sum((s - p) * w) / sqrt(sum(p * (1 - p) * w^2))
}
