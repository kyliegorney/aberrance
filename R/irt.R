#' Compute IRT probabilities
#'
#' @noRd

irt_p <- function(m, psi, xi, ignore = NULL) {
  psi <- psi[, !colnames(psi) %in% ignore, drop = FALSE]
  n <- nrow(psi)
  p <- array(dim = c(nrow(xi), n, max(m)))
  if ("b" %in% colnames(psi)) {
    for (i in 1:n) {
      p[, i, m[i]] <- psi[i, "c"] + (1 - psi[i, "c"]) /
        (1 + exp(psi[i, "a"] * (psi[i, "b"] - xi[, "theta"])))
    }
    if ("lambda1" %in% colnames(psi)) { # nested logit model
      for (i in 1:n) {
        p[, i, 1:(m[i]-1)] <- (1 - p[, i, m[i]]) * t(apply(
          t(exp(outer(psi[i, paste0("lambda", 1:(m[i]-1))], xi[, "eta"], "*") +
                  psi[i, paste0("zeta", 1:(m[i]-1))])), 1, proportions))
      }
    } else { # 3PL model
      p[, , 1] <- 1 - p[, , 2]
    }
  } else if ("b1" %in% colnames(psi)) { # graded response model
    for (i in 1:n) {
      tmp <- cbind(
        1,
        1 / t(1 + exp(psi[i, "a"] * outer(
          psi[i, paste0("b", 1:(m[i]-1))], xi[, "theta"], "-"))),
        0
      )
      p[, i, 1:m[i]] <- tmp[, 1:m[i]] - tmp[, 2:(m[i]+1)]
    }
  } else if ("c0" %in% colnames(psi)) { # generalized partial credit model
    for (i in 1:n) {
      p[, i, 1:m[i]] <- t(apply(
        -t(psi[i, "a"] * outer(
          psi[i, paste0("c", 0:(m[i]-1))], xi[, "theta"], "-")), 1, function(j)
            proportions(exp(cumsum(j)))))
    }
  } else { # nominal response model
    for (i in 1:n) {
      p[, i, 1:m[i]] <- t(apply(
        t(exp(outer(psi[i, paste0("lambda", 1:m[i])], xi[, "eta"], "*") +
                psi[i, paste0("zeta", 1:m[i])])), 1, proportions))
    }
  }
  p
}

#' Compute the first derivative of p
#'
#' @noRd

irt_p1 <- function(p, m, psi, xi, ignore = NULL) {
  psi <- psi[, !colnames(psi) %in% ignore, drop = FALSE]
  n <- nrow(psi)
  p1 <- array(dim = dim(p))
  if ("b" %in% colnames(psi)) { # 3PL model
    e <- exp(-psi[, "a"] * outer(psi[, "b"], xi[, "theta"], "-"))
    p1[, , 2] <- t((1 - psi[, "c"]) * psi[, "a"] * e / (1 + e)^2)
    p1[, , 1] <- -p1[, , 2]
  } else if ("b1" %in% colnames(psi)) { # graded response model
    for (i in 1:n) {
      tmp <- cbind(
        1,
        1 / t(1 + exp(psi[i, "a"] * outer(
          psi[i, paste0("b", 1:(m[i]-1))], xi[, "theta"], "-"))),
        0
      )
      p1[, i, 1:m[i]] <- psi[i, "a"] *
        (tmp[, 1:m[i]] * (1 - tmp[, 1:m[i]]) -
           tmp[, 2:(m[i]+1)] * (1 - tmp[, 2:(m[i]+1)]))
    }
  } else if ("c0" %in% colnames(psi)) { # generalized partial credit model
    for (i in 1:n) {
      inn <- rowSums(t(0:(m[i]-1) * t(p[, i, 1:m[i]])))
      p1[, i, 1:m[i]] <- psi[i, "a"] * p[, i, 1:m[i]] *
        t(outer(0:(m[i]-1), inn, "-"))
    }
  } else { # nominal response model
    for (i in 1:n) {
      num <- t(exp(outer(psi[i, paste0("lambda", 1:m[i])], xi[, "eta"], "*") +
                     psi[i, paste0("zeta", 1:m[i])]))
      den <- rowSums(num)
      for (j in 1:m[i]) {
        p1[, i, j] <- num[, j] *
          rowSums(t(t(num) * (psi[i, paste0("lambda", j)] -
                                psi[i, paste0("lambda", 1:m[i])]))) / den^2
      }
    }
  }
  p1
}

#' Compute the log-likelihood
#'
#' @noRd

irt_l <- function(x, p) {
  N <- dim(p)[1]
  n <- dim(p)[2]
  m <- dim(p)[3]
  L <- matrix(nrow = N, ncol = n)
  for (i in 1:n) {
    I <- outer(x[, i], 0:(m-1), "==")
    L[, i] <- rowSums(I * p[, i, ], na.rm = TRUE)
  }
  log(L)
}

#' Compute the first derivative of l
#'
#' @noRd

irt_l1 <- function(x, p, p1) {
  N <- dim(p)[1]
  n <- dim(p)[2]
  m <- dim(p)[3]
  l1 <- matrix(nrow = N, ncol = n)
  for (i in 1:n) {
    I <- outer(x[, i], 0:(m-1), "==")
    l1[, i] <- rowSums(I * p1[, i, ] / p[, i, ], na.rm = TRUE)
  }
  l1
}

#' Compute item and test information
#'
#' @noRd

irt_info <- function(p, p1, tif = TRUE) {
  info <- apply(p1^2 / p, 1:2, sum, na.rm = TRUE)
  if (tif) {
    rowSums(info)
  } else {
    info
  }
}
