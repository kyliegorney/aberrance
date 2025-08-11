#' Detect test tampering
#'
#' @description Detect test tampering at the person level or at the group level.
#'
#' @param method The test tampering statistic(s) to compute. Options for score
#'   and distractor-based statistics are:
#'   - `"EDI_SD_*"` for the erasure detection index (Wollack et al., 2015;
#'     Wollack & Eckerly, 2017).
#'   - `"GBT_SD"` for the generalized binomial test statistic (Sinharay &
#'     Johnson, 2017). *Note:* This statistic cannot be computed at the group
#'     level.
#'   - `"L_SD"` for the signed likelihood ratio test statistic (Sinharay et al.,
#'     2017). *Note:* This statistic cannot be computed at the group level.
#'
#'   Options for response-based statistics are:
#'   - `"EDI_R_*"` for the erasure detection index (Wollack et al., 2015;
#'     Wollack & Eckerly, 2017).
#'   - `"GBT_R"` for the generalized binomial test statistic (Sinharay &
#'     Johnson, 2017). *Note:* This statistic cannot be computed at the group
#'     level.
#'   - `"L_R"` for the signed likelihood ratio test statistic (Sinharay et al.,
#'     2017). *Note:* This statistic cannot be computed at the group level.
#'
#'   Statistics ending in `"*"` can be computed using various corrections.
#'   Options are:
#'   - `"*"` for all possible corrections.
#'   - `"NO"` for no correction (Sinharay, 2018; Wollack et al., 2015).
#'   - `"CO"` for the continuity correction (Wollack et al., 2015; Wollack &
#'     Eckerly, 2017). The value of the continuity correction can be specified
#'     using `c`.
#'   - `"TS"` for the Taylor series expansion (Sinharay, 2018).
#'
#' @param xi,xi_c,xi_s Matrices of person parameters. Rows correspond to persons
#'   and columns to parameters. `xi` is based on all items, `xi_c` is based on
#'   items with changed answers, and `xi_s` is based on items with the same
#'   answers. If `NULL` (default), person parameters are estimated using maximum
#'   likelihood estimation.
#'
#' @param x,d,r Matrices of final data. Rows correspond to persons and columns
#'   to items. `x` is for the item scores, `d` the item distractors, and `r` the
#'   item responses.
#'
#' @param x_0,d_0,r_0 Matrices of initial data. Rows correspond to persons and
#'   columns to items. `x_0` is for the item scores, `d_0` the item distractors,
#'   and `r_0` the item responses.
#'
#' @param group A vector indicating group membership. If `NULL` (default),
#'   statistics are computed at the person level.
#'
#' @param c Use with the erasure detection index. A value indicating the
#'   continuity correction. Default is `-0.5`.
#'
#' @inheritParams detect_pm
#'
#' @returns A list is returned with the following elements:
#' \item{stat}{A matrix of test tampering detection statistics.}
#' \item{pval}{A matrix of *p*-values.}
#' \item{flag}{An array of flagging results. The first dimension corresponds to
#'   persons/groups, the second dimension to methods, and the third dimension to
#'   significance levels.}
#'
#' @references
#' Sinharay, S., Duong, M. Q., & Wood, S. W. (2017). A new statistic for
#' detection of aberrant answer changes. *Journal of Educational Measurement*,
#' *54*(2), 200--217.
#'
#' Sinharay, S., & Johnson, M. S. (2017). Three new methods for analysis of
#' answer changes. *Educational and Psychological Measurement*, *77*(1), 54--81.
#'
#' Sinharay, S. (2018). Detecting fraudulent erasures at an aggregate level.
#' *Journal of Educational and Behavioral Statistics*, *43*(3), 286--315.
#'
#' Wollack, J. A., Cohen, A. S., & Eckerly, C. A. (2015). Detecting test
#' tampering using item response theory. *Educational and Psychological
#' Measurement*, *75*(6), 931--953.
#'
#' Wollack, J. A., & Eckerly, C. A. (2017). Detecting test tampering at the
#' group level. In G. J. Cizek & J. A. Wollack (Eds.), *Handbook of quantitative
#' methods for detecting cheating on tests* (pp. 214--231). Routledge.
#'
#' @examples
#' # Setup for Examples 1 and 2 ------------------------------------------------
#'
#' # Settings
#' set.seed(0)     # seed for reproducibility
#' N <- 500        # number of persons
#' n <- 40         # number of items
#' G <- 20         # number of groups
#'
#' # Create groups
#' group <- rep(1:G, each = N / G)
#'
#' # Randomly select 20% tampered groups with 20% tampered persons
#' cg <- sample(1:G, size = G * 0.20)
#' cv <- NULL
#' for (g in cg) {
#'   cv <- c(cv, sample(which(group == g), size = N / G * 0.20))
#' }
#'
#' # Create vectors of indicators (1 = tampered, 0 = non-tampered)
#' group_ind <- ifelse(1:G %in% cg, 1, 0)
#' person_ind <- ifelse(1:N %in% cv, 1, 0)
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
#' x_0 <- x <- dat$x
#' d_0 <- d <- dat$d
#'
#' # Simulate 5% random erasures for non-tampered persons
#' r_0 <- r <- ifelse(x == 1, 4, d)
#' for (v in setdiff(1:N, cv)) {
#'   ci <- sample(1:n, size = n * 0.05)
#'   for (i in ci) {
#'     r_0[v, i] <- sample((1:4)[-r[v, i]], size = 1)
#'   }
#'   x_0[v, ci] <- ifelse(r_0[v, ci] == 4, 1, 0)
#'   d_0[v, ci] <- ifelse(r_0[v, ci] == 4, NA, r_0[v, ci])
#' }
#' rm(r_0, r)
#'
#' # Modify contaminated data by tampering with 20% of the scores and
#' # distractors
#' for (v in cv) {
#'   ci <- sample(1:n, size = n * 0.20)
#'   x[v, ci] <- 1
#'   d[v, ci] <- NA
#' }
#'
#' # Example 1: Person-Level Statistics ----------------------------------------
#'
#' # Detect test tampering
#' out <- detect_tt(
#'   method = c("EDI_SD_*", "GBT_SD", "L_SD"),
#'   psi = psi,
#'   x = x,
#'   d = d,
#'   x_0 = x_0,
#'   d_0 = d_0
#' )
#'
#' # Example 2: Group-Level Statistics -----------------------------------------
#'
#' # Detect test tampering
#' out <- detect_tt(
#'   method = "EDI_SD_*",
#'   psi = psi,
#'   x = x,
#'   d = d,
#'   x_0 = x_0,
#'   d_0 = d_0,
#'   group = group
#' )
#' @export

detect_tt <- function(method,
                      psi,
                      xi = NULL,
                      xi_c = NULL,
                      xi_s = NULL,
                      x = NULL,
                      d = NULL,
                      r = NULL,
                      x_0 = NULL,
                      d_0 = NULL,
                      r_0 = NULL,
                      interval = c(-4, 4),
                      alpha = 0.05,
                      group = NULL,
                      c = -0.5) {

  # Checks
  if (any(extract(method, 3) == "*")) {
    tmp <- NULL
    for (m in method) {
      if (extract(m, 3) == "*") {
        if (is.null(group)) {
          tmp <- c(
            tmp,
            paste(
              extract(m, 1:2),
              c("NO", "CO"),
              sep = "_"
            )
          )
        } else {
          tmp <- c(
            tmp,
            paste(
              extract(m, 1:2),
              c("NO", "CO", "TS"),
              sep = "_"
            )
          )
        }
      } else {
        tmp <- c(tmp, m)
      }
    }
    method <- tmp
  }
  if (any("SD" %in% extract(method, 2)) && any("R" %in% extract(method, 2))) {
    stop(
      "`method` may contain either score and distractor-based statistics or ",
      "response-based statistics, but not both.",
      call. = FALSE
    )
  }
  if (any("SD" %in% extract(method, 2))) {
    check_par("x", psi)
    if (is.null(group)) {
      if ("EDI_SD_TS" %in% method) {
        method <- setdiff(method, "EDI_SD_TS")
        warning(
          "The EDI_SD_TS statistic cannot be computed at the person level.",
          call. = FALSE
        )
      }
    } else {
      if ("GBT_SD" %in% method) {
        method <- setdiff(method, "GBT_SD")
        warning(
          "The GBT_SD statistic cannot be computed at the group level.",
          call. = FALSE
        )
      }
      if ("L_SD" %in% method) {
        method <- setdiff(method, "L_SD")
        warning(
          "The L_SD statistic cannot be computed at the group level.",
          call. = FALSE
        )
      }
    }
  } else if (any("R" %in% extract(method, 2))) {
    check_par("r", psi)
    if (is.null(group)) {
      if ("EDI_R_TS" %in% method) {
        method <- setdiff(method, "EDI_R_TS")
        warning(
          "The EDI_R_TS statistic cannot be computed at the person level.",
          call. = FALSE
        )
      }
    } else {
      if ("GBT_R" %in% method) {
        method <- setdiff(method, "GBT_R")
        warning(
          "The GBT_R statistic cannot be computed at the group level.",
          call. = FALSE
        )
      }
      if ("L_R" %in% method) {
        method <- setdiff(method, "L_R")
        warning(
          "The L_R statistic cannot be computed at the group level.",
          call. = FALSE
        )
      }
    }
  }
  check_data(x, d, r)

  # Setup
  if (is.null(group)) {
    N <- max(nrow(x), nrow(d), nrow(r))
    method <- match.arg(
      arg = unique(method),
      choices = c(
        t(outer(
          c("EDI_SD", "EDI_R"),
          c("NO", "CO"),
          paste, sep = "_"
        )),
        "GBT_SD", "L_SD", "GBT_R", "L_R"
      ),
      several.ok = TRUE
    )
    stat <- pval <- matrix(
      nrow = N,
      ncol = length(method),
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
  } else {
    N <- length(table(group))
    method <- match.arg(
      arg = unique(method),
      choices = c(
        t(outer(
          c("EDI_SD", "EDI_R"),
          c("NO", "CO", "TS"),
          paste, sep = "_"
        ))
      ),
      several.ok = TRUE
    )
    stat <- pval <- matrix(
      nrow = N,
      ncol = length(method),
      dimnames = list(
        group = names(table(group)),
        method = method
      )
    )
    flag <- array(
      dim = c(N, length(method), length(alpha)),
      dimnames = list(
        group = row.names(stat),
        method = method,
        alpha = alpha
      )
    )
  }
  n <- max(ncol(x), ncol(d), ncol(r))
  if ("SD" %in% extract(method, 2)) {
    m <- count(psi, ignore = "lambda1")
    x_s <- ifelse(
      ((x == 1) & (x_0 == 1)) | ((x == 0) & (x_0 == 0) & (d == d_0)),
      x,
      NA
    )
    xi_s <- est(interval, psi, x = x_s)
    p_s <- irt_p(m, psi, xi_s, ignore = "lambda1")
  } else {
    m <- count(psi)
    r_s <- ifelse(r == r_0, r, NA)
    xi_s <- est(interval, psi, r = r_s)
    p_s <- irt_p(m, psi, xi_s)
  }

  # Compute person-level statistics
  if (is.null(group)) {

    # Compute score and distractor-based statistics
    if (any(c("EDI_SD_NO", "EDI_SD_CO", "GBT_SD") %in% method)) {
      for (v in which(rowSums(is.na(x_s)) > 0)) {
        ci <- which(is.na(x_s[v, ]))
        s <- x[v, ci]
        p <- p_s[v, ci, 2]
        if ("EDI_SD_NO" %in% method) {
          stat[v, "EDI_SD_NO"] <- compute_OMG(s, p, c = 0)
        }
        if ("EDI_SD_CO" %in% method) {
          stat[v, "EDI_SD_CO"] <- compute_OMG(s, p, c = c)
        }
        if ("GBT_SD" %in% method) {
          stat[v, "GBT_SD"] <- compute_GBT(s, p)
        }
      }
    }
    if ("L_SD" %in% method) {
      x_c <- ifelse(is.na(x_s), x, NA)
      xi_c <- est(interval, psi, x = x_c)
      p_c <- irt_p(m, psi, xi_c, ignore = "lambda1")
      xi <- est(interval, psi, x = x)
      p_0 <- p_1 <- irt_p(m, psi, xi, ignore = "lambda1")
      p_1[, , 1] <- ifelse(is.na(x_s), p_c[, , 1], p_s[, , 1])
      p_1[, , 2] <- ifelse(is.na(x_s), p_c[, , 2], p_s[, , 2])
      stat[, "L_SD"] <- compute_L_S(
        x,
        p_0,
        p_1,
        xi_c[, "theta"],
        xi_s[, "theta"],
        signed = TRUE
      )
      stat[which(rowSums(is.na(x_s)) == 0), "L_SD"] <- NA
    }

    # Compute response-based statistics
    if (any(c("EDI_R_NO", "EDI_R_CO", "GBT_R") %in% method)) {
      for (v in which(rowSums(is.na(r_s)) > 0)) {
        ci <- which(is.na(r_s[v, ]))
        s <- as.integer(r[v, ci] == m[ci])
        p <- rep(NA, times = length(ci))
        for (i in 1:length(ci)) {
          p[i] <- p_s[v, ci[i], m[ci[i]]]
        }
        if ("EDI_R_NO" %in% method) {
          stat[v, "EDI_R_NO"] <- compute_OMG(s, p, c = 0)
        }
        if ("EDI_R_CO" %in% method) {
          stat[v, "EDI_R_CO"] <- compute_OMG(s, p, c = c)
        }
        if ("GBT_R" %in% method) {
          stat[v, "GBT_R"] <- compute_GBT(s, p)
        }
      }
    }
    if ("L_R" %in% method) {
      r_c <- ifelse(is.na(r_s), r, NA)
      xi_c <- est(interval, psi, r = r_c)
      p_c <- irt_p(m, psi, xi_c)
      xi <- est(interval, psi, r = r)
      p_0 <- p_1 <- irt_p(m, psi, xi)
      for (j in 1:max(m)) {
        p_1[, , j] <- ifelse(is.na(r_s), p_c[, , j], p_s[, , j])
      }
      stat[, "L_R"] <- compute_L_S(
        r - 1,
        p_0,
        p_1,
        xi_c[, "eta"],
        xi_s[, "eta"],
        signed = TRUE
      )
      stat[which(rowSums(is.na(r_s)) == 0), "L_R"] <- NA
    }

  # Compute group-level statistics
  } else {

    # Compute score and distractor-based statistics
    if (any(c("EDI_SD_NO", "EDI_SD_CO", "EDI_SD_TS") %in% method)) {
      if ("EDI_SD_TS" %in% method) {
        p1_s <- irt_p1(p_s, m, psi, xi_s, ignore = "lambda1")
        info_s <- irt_info(p_s, p1_s, tif = FALSE)
      }
      for (v in 1:N) {
        g <- which(group == names(table(group))[v])
        s <- ifelse(is.na(x_s[g, ]), x[g, ], NA)
        if (sum(s, na.rm = TRUE) > 0) {
          p <- ifelse(is.na(x_s[g, ]), p_s[g, , 2], NA)
          if ("EDI_SD_NO" %in% method) {
            stat[v, "EDI_SD_NO"] <- compute_EDI_CO(s, p, c = 0)
          }
          if ("EDI_SD_CO" %in% method) {
            stat[v, "EDI_SD_CO"] <- compute_EDI_CO(s, p, c = c)
          }
          if ("EDI_SD_TS" %in% method) {
            p1 <- ifelse(is.na(x_s[g, ]), p1_s[g, , 2], NA)
            info <- rowSums(
              ifelse(!is.na(x_s[g, ]), info_s[g, ], NA), na.rm = TRUE)
            stat[v, "EDI_SD_TS"] <- compute_EDI_TS(s, p, p1, info)
          }
        }
      }
    }

    # Compute response-based statistics
    if (any(c("EDI_R_NO", "EDI_R_CO", "EDI_R_TS") %in% method)) {
      if ("EDI_R_TS" %in% method) {
        p1_s <- irt_p1(p_s, m, psi, xi_s)
        info_s <- irt_info(p_s, p1_s, tif = FALSE)
      }
      for (v in 1:N) {
        g <- which(group == names(table(group))[v])
        s <- ifelse(is.na(r_s[g, ]), t(t(r[g, ]) == m), NA)
        if (sum(s, na.rm = TRUE) > 0) {
          p <- matrix(nrow = length(g), ncol = n)
          for (i in 1:n) {
            p[, i] <- ifelse(is.na(r_s[g, i]), p_s[g, i, m[i]], NA)
          }
          if ("EDI_R_NO" %in% method) {
            stat[v, "EDI_R_NO"] <- compute_EDI_CO(s, p, c = 0)
          }
          if ("EDI_R_CO" %in% method) {
            stat[v, "EDI_R_CO"] <- compute_EDI_CO(s, p, c = c)
          }
          if ("EDI_R_TS" %in% method) {
            p1 <- matrix(nrow = length(g), ncol = n)
            for (i in 1:n) {
              p1[, i] <- ifelse(is.na(r_s[g, i]), p1_s[g, i, m[i]], NA)
            }
            info <- rowSums(
              ifelse(!is.na(r_s[g, ]), info_s[g, ], NA), na.rm = TRUE)
            stat[v, "EDI_R_TS"] <- compute_EDI_TS(s, p, p1, info)
          }
        }
      }
    }
  }

  # Compute p-values
  pval[, grep("EDI", colnames(pval))] <-
    pnorm(stat[, grep("EDI", colnames(stat))], lower.tail = FALSE)
  pval[, grep("GBT", colnames(pval))] <-
    stat[, grep("GBT", colnames(stat))]
  pval[, grep("L", colnames(pval))] <-
    pnorm(stat[, grep("L", colnames(stat))], lower.tail = FALSE)

  # Compute flagging rates
  for (a in 1:length(alpha)) {
    flag[, , a] <- pval <= alpha[a]
  }

  # Output
  list(stat = stat, pval = pval, flag = flag)
}
