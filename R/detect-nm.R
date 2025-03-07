#' Detect nonparametric misfit
#'
#' @description Detect nonparametric misfit using person-fit statistics.
#'
#' @param method The person-fit statistic(s) to compute. Options for score-based
#'   statistics are:
#'   - `"G_S"` for the number of Guttman errors (Guttman, 1944; see also
#'     Molenaar, 1991).
#'   - `"NC_S"` for the norm conformity index (Tatsuoka & Tatsuoka, 1983).
#'     *Note:* This statistic cannot be computed for polytomous item scores.
#'   - `"U1_S"` for the \eqn{U1} statistic, also known as the \eqn{G^*}
#'     statistic (van der Flier, 1977; see also Emons, 2008).
#'   - `"U3_S"` for the \eqn{U3} statistic (van der Flier, 1982; see also Emons,
#'     2008).
#'   - `"ZU3_S"` for the \eqn{ZU3} statistic (van der Flier, 1982).
#'     *Note:* This statistic cannot be computed for polytomous item scores.
#'   - `"A_S"` for the agreement index (Kane & Brennan, 1980). *Note:* This
#'     statistic cannot be computed for polytomous item scores.
#'   - `"D_S"` for the disagreement index (Kane & Brennan, 1980). *Note:* This
#'     statistic cannot be computed for polytomous item scores.
#'   - `"E_S"` for the dependability index (Kane & Brennan, 1980). *Note:* This
#'     statistic cannot be computed for polytomous item scores.
#'   - `"C_S"` for the caution index (Sato, 1975). *Note:* This statistic cannot
#'     be computed for polytomous item scores.
#'   - `"MC_S"` for the modified caution index, also known as the \eqn{C^*}
#'     statistic (Harnisch & Linn, 1981). *Note:* This statistic cannot be
#'     computed for polytomous item scores.
#'   - `"PC_S"` for the personal point-biserial correlation (Donlon & Fischer,
#'     1968). *Note:* This statistic cannot be computed for polytomous item
#'     scores.
#'   - `"HT_S` for the \eqn{H^T} statistic (Sijtsma, 1986). *Note:* This
#'     statistic cannot be computed for polytomous item scores.
#'
#'   Options for response time-based statistics are:
#'   - `"KL_T"` for the Kullback-Leibler divergence (Man et al., 2018).
#'
#' @param x,y Matrices of raw data. `x` is for the item scores and `y` the item
#'   log response times.
#'
#' @returns A list is returned with the following elements:
#' \item{stat}{A matrix of nonparametric person-fit statistics.}
#'
#' @references
#' Donlon, T. F., & Fischer, F. E. (1968). An index of an individual's agreement
#' with group-determined item difficulties. *Educational and Psychological
#' Measurement*, *28*(1), 105--113.
#'
#' Emons, W. H. M. (2008). Nonparametric person-fit analysis of polytomous item
#' scores. *Applied Psychological Measurement*, *32*(3), 224--247.
#'
#' Guttman, L. (1944). A basis for scaling qualitative data. *American
#' Sociological Review*, *9*(2), 139--150.
#'
#' Harnisch, D. L., & Linn, R. L. (1981). Analysis of item response patterns:
#' Questionable test data and dissimilar curriculum practices. *Journal of
#' Educational Measurement*, *18*(3), 133--146.
#'
#' Kane, M. T., & Brennan, R. L. (1980). Agreement coefficients as indices of
#' dependability for domain referenced tests. *Applied Psychological
#' Measurement*, *4*(1), 105--126.
#'
#' Man, K., Harring, J. R., Ouyang, Y., & Thomas, S. L. (2018). Response time
#' based nonparametric Kullback-Leibler divergence measure for detecting
#' aberrant test-taking behavior. *International Journal of Testing*, *18*(2),
#' 155--177.
#'
#' Molenaar, I. W. (1991). A weighted Loevinger H-coefficient extending Mokken
#' scaling to multicategory items. *Kwantitatieve Methoden*, *12*(37), 97--117.
#'
#' Sato, T. (1975). *The construction and interpretation of S-P tables*.
#'
#' Sijtsma, K. (1986). A coefficient of deviance of response patterns.
#' *Kwantitatieve Methoden*, *7*(22), 131--145.
#'
#' Tatsuoka, K. K., & Tatsuoka, M. M. (1983). Spotting erroneous rules of
#' operation by the individual consistency index. *Journal of Educational
#' Measurement*, *20*(3), 221--230.
#'
#' van der Flier, H. (1977) Environmental factors and deviant response patterns.
#' In Y. H. Poortinga (Ed.), *Basic problems in cross-cultural psychology*.
#' Swets & Zeitlinger Publishers.
#'
#' van der Flier, H. (1982). Deviant response patterns and comparability of test
#' scores. *Journal of Cross-Cultural Psychology*, *13*(3), 267--298.
#'
#' @seealso [detect_pm()] to detect parametric misfit.
#'
#' @examples
#' # Setup for Examples 1 to 3 -------------------------------------------------
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
#' # Example 1: Dichotomous Item Scores ----------------------------------------
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
#' # Modify contaminated data by changing the item scores
#' x[cv, ci] <- rbinom(length(cv) * length(ci), size = 1, prob = 0.90)
#'
#' # Detect nonparametric misfit
#' out <- detect_nm(
#'   method = c("G_S", "NC_S", "U1_S", "U3_S", "ZU3_S", "A_S", "D_S", "E_S",
#'              "C_S", "MC_S", "PC_S", "HT_S"),
#'   x = x
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
#' # Detect nonparametric misfit
#' out <- detect_nm(
#'   method = c("G_S", "U1_S", "U3_S"),
#'   x = x
#' )
#'
#' # Example 3: Item Response Times --------------------------------------------
#'
#' # Generate person parameters for the lognormal model
#' xi <- cbind(tau = rnorm(N, mean = 0.00, sd = sqrt(0.25)))
#'
#' # Generate item parameters for the lognormal model
#' psi <- cbind(
#'   alpha = runif(n, min = 1.50, max = 2.50),
#'   beta = rnorm(n, mean = 3.50, sd = sqrt(0.15))
#' )
#'
#' # Simulate uncontaminated data
#' y <- sim(psi, xi)$y
#'
#' # Modify contaminated data by reducing the log response times
#' y[cv, ci] <- y[cv, ci] * 0.75
#'
#' # Detect nonparametric misfit
#' out <- detect_nm(
#'   method = "KL_T",
#'   y = y
#' )
#' @export

detect_nm <- function(method, x = NULL, y = NULL) {

  # Checks
  if (any(c("G_S", "NC_S", "U1_S", "U3_S", "ZU3_S", "A_S", "D_S", "E_S",
            "C_S", "MC_S", "PC_S", "HT_S") %in% method)) {
    if (max(x) > 1 &&
        any(c("NC_S", "ZU3_S", "A_S", "D_S", "E_S",
              "C_S", "MC_S", "PC_S", "HT_S") %in% method)) {
      method <- setdiff(method, c("NC_S", "ZU3_S", "A_S", "D_S", "E_S",
                                  "C_S", "MC_S", "PC_S", "HT_S"))
      warning("The NC_S, ZU3_S, A_S, D_S, E_S, C_S, MC_S, PC_S, and HT_S ",
              "statistics cannot be computed for polytomous item scores.",
              call. = FALSE)
    }
  }
  method <- match.arg(
    arg = unique(method),
    choices = c("G_S",  "NC_S", "U1_S", "U3_S", "ZU3_S", "A_S", "D_S", "E_S",
                "C_S", "MC_S", "PC_S", "HT_S", "KL_T"),
    several.ok = TRUE
  )
  check_data(x, y)

  # Setup
  N <- max(nrow(x), nrow(y))
  n <- max(ncol(x), ncol(y))
  stat <- matrix(
    nrow = N, ncol = length(method),
    dimnames = list(
      person = 1:N,
      method = method
    )
  )

  # Compute score-based statistics
  if (any(c("G_S", "NC_S", "U1_S", "U3_S", "ZU3_S", "A_S", "D_S", "E_S",
            "C_S", "MC_S", "PC_S", "HT_S") %in% method)) {
    m <- count(x = x)
    max_m <- max(m)
    if (max_m == 2) {
      p <- colMeans(x)
      s <- rowSums(x)
      incl <- which((s > 0) & (s < n))
      excl <- setdiff(1:N, incl)
      gutt <- order(p, decreasing = TRUE)
      if (any(c("G_S", "NC_S", "U1_S") %in% method)) {
        G_S <- U1_S <- rep(0, times = length(incl))
        for (i in 1:(n-1)) {
          G_S <- G_S + rowSums(cbind(x[incl, gutt[i]] < x[incl, gutt[(i+1):n]]))
        }
        U1_S <- G_S / (s[incl] * (n - s[incl]))
        if ("G_S" %in% method) {
          stat[incl, "G_S"] <- G_S
        }
        if ("U1_S" %in% method) {
          stat[incl, "U1_S"] <- U1_S
        }
        if ("NC_S" %in% method) {
          stat[incl, "NC_S"] <- 1 - 2 * U1_S
        }
      }
      if (any(c("U3_S", "ZU3_S") %in% method)) {
        w <- log(p / (1 - p))
        w[is.infinite(w)] <- 0
        W <- rowSums(t(t(x[incl, ]) * w))
        W_max <- cumsum(w[gutt])[s[incl]]
        W_min <- cumsum(w[rev(gutt)])[s[incl]]
        U3_S <- (W_max - W) / (W_max - W_min)
        if ("U3_S" %in% method) {
          stat[incl, "U3_S"] <- U3_S
        }
        if ("ZU3_S" %in% method) {
          A <- sum(p * w) + sum(p * (1 - p) * w) * (s[incl] - sum(p)) /
            sum(p * (1 - p))
          B <- sum(p * (1 - p) * w^2) - sum(p * (1 - p) * w)^2 /
            sum(p * (1 - p))
          mu <- (W_max - A) / (W_max - W_min)
          sigma <- sqrt(B) / abs(W_max - W_min)
          stat[incl, "ZU3_S"] <- (U3_S - mu) / sigma
        }
      }
      if (any(c("A_S", "D_S", "E_S", "C_S", "MC_S") %in% method)) {
        P <- rowSums(t(t(x[incl, ]) * p))
        P_max <- cumsum(p[gutt])[s[incl]]
        P_min <- cumsum(p[rev(gutt)])[s[incl]]
        if ("A_S" %in% method) {
          stat[incl, "A_S"] <- P
        }
        if ("D_S" %in% method) {
          stat[incl, "D_S"] <- P_max - P
        }
        if ("E_S" %in% method) {
          stat[incl, "E_S"] <- P / P_max
        }
        if ("C_S" %in% method) {
          stat[incl, "C_S"] <-
            (n * P_max - n * P) / (n * P_max - s[incl] * sum(p))
        }
        if ("MC_S" %in% method) {
          stat[incl, "MC_S"] <- (P_max - P) / (P_max - P_min)
        }
      }
      if ("PC_S" %in% method) {
        stat[incl, "PC_S"] <- cor(t(x[incl, ]), cbind(p))
      }
      if ("HT_S" %in% method) {
        cov_mat <- cov(t(x))
        s_mat <- outer(s / n, 1 - s / n, "*")
        diag(cov_mat) <- diag(s_mat) <- NA
        stat[incl, "HT_S"] <-
          rowSums(cov_mat[incl, ] * (n - 1) / n, na.rm = TRUE) /
          rowSums(pmin(s_mat, t(s_mat)), na.rm = TRUE)[incl]
      }
    } else {
      n_step <- sum(m - 1)
      x_step <- matrix(nrow = N, ncol = n_step)
      ms <- cumsum(c(0, m - 1))
      for (i in 1:n) {
        for (j in 1:(m[i]-1)) {
          x_step[, ms[i]+j] <- ifelse(x[, i] >= j, 1, 0)
        }
      }
      p <- colMeans(x_step)
      s <- rowSums(x_step)
      incl <- which((s > 0) & (s < n_step))
      excl <- setdiff(1:N, incl)
      gutt <- order(p, decreasing = TRUE)
      if (any(c("G_S", "U1_S") %in% method)) {
        G_S <- rep(0, times = length(incl))
        for (i in 1:(n_step-1)) {
          G_S <- G_S + rowSums(cbind(
            x_step[incl, gutt[i]] < x_step[incl, gutt[(i+1):n_step]]))
        }
        if ("G_S" %in% method) {
          stat[incl, "G_S"] <- G_S
        }
        if ("U1_S" %in% method) {
          w <- rank(-p, ties.method = "first")
          R <- matrix(0, nrow = n, ncol = max_m)
          for (i in 1:n) {
            R[i, 2:m[i]] <- w[(ms[i]+1):ms[i+1]]
          }
          R <- t(apply(R, 1, cumsum))
          V <- outer(R[1, 1:m[1]], R[2, 1:m[2]], "+")
          T <- sapply(2:sum(dim(V)), function(j) max(V[row(V) + col(V) == j]))
          for (i in 2:(n-1)) {
            V <- outer(T, R[i+1, 1:m[i+1]], "+")
            T <- sapply(2:sum(dim(V)), function(j) max(V[row(V) + col(V) == j]))
          }
          stat[incl, "U1_S"] <- G_S /
            (T[s[incl]+1] - 0.5 * s[incl] * (s[incl] + 1))
        }
      }
      if ("U3_S" %in% method) {
        w <- log(p / (1 - p))
        w[is.infinite(w)] <- 0
        W <- rowSums(t(t(x_step[incl, ]) * w))
        W_max <- cumsum(w[gutt])[s[incl]]
        R <- matrix(0, nrow = n, ncol = max_m)
        for (i in 1:n) {
          R[i, 2:m[i]] <- w[(ms[i]+1):ms[i+1]]
        }
        R <- t(apply(R, 1, cumsum))
        V <- outer(R[1, 1:m[1]], R[2, 1:m[2]], "+")
        T <- sapply(2:sum(dim(V)), function(j) min(V[row(V) + col(V) == j]))
        for (i in 2:(n-1)) {
          V <- outer(T, R[i+1, 1:m[i+1]], "+")
          T <- sapply(2:sum(dim(V)), function(j) min(V[row(V) + col(V) == j]))
        }
        stat[incl, "U3_S"] <- (W_max - W) / (W_max - T[s[incl]+1])
      }
    }
  }

  # Compute response time-based statistics
  if ("KL_T" %in% method) {
    f <- proportions(colMeans(exp(y)))
    g <- apply(exp(y), 1, proportions)
    stat[, "KL_T"] <- rowSums(t(f * log(f / g)))
  }

  # Output
  list(stat = stat)
}
