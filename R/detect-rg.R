#' Detect rapid guessing
#'
#' @description Detect rapid guessing using item-level response time
#'   information.
#'
#' @param method The rapid guessing detection method to apply. Options for
#'   visual inspection methods are:
#'   - `"VI"` for the visual inspection method (Schnipke, 1995). Each plot
#'     contains a histogram of the item response times.
#'   - `"VITP"` for the visual inspection with proportion correct method (Lee &
#'     Jia, 2014; Ma et al., 2011). Each plot contains a histogram of the item
#'     response times, a dashed red line indicating the proportion correct, and
#'     a solid red line indicating the chance rate of success (see `chance`).
#'
#'   Options for threshold methods are:
#'   - `"CT"` for the custom threshold method (Wise et al., 2004; Wise & Kong,
#'     2005). The thresholds can be specified using `thr`.
#'   - `"NT"` for the normative threshold method (Martinez & Rios, 2023; Wise &
#'     Ma, 2012). The percentage(s) of the mean item response time can be
#'     specified using `nt`.
#'
#'   Options for visual inspection and threshold methods are:
#'   - `"CUMP"` for the cumulative proportion correct method (Guo et al., 2016).
#'     Each plot contains a histogram of the item response times, a dashed red
#'     line indicating the cumulative proportion correct, and a solid red line
#'     indicating the chance rate of success (see `chance`). *Note:* No
#'     thresholds are returned for items for which the cumulative proportion
#'     correct is consistently above or below `chance`.
#'
#' @param t,x Matrices of raw data. Rows correspond to persons and columns to
#'   items. `t` is for the item response times and `x` the item scores.
#'
#' @param outlier The percentile(s) above which to delete outliers in `t`.
#'   Length must be equal to 1 or equal to the total number of items. Default is
#'   `100`, for which no response times are identified as outliers to be
#'   deleted.
#'
#' @param chance Use with the visual inspection with proportion correct method
#'   and the cumulative proportion correct method. Value(s) indicating the
#'   chance rate(s) of success. Length must be equal to 1 or equal to the total
#'   number of items. Default is `0.25`.
#'
#' @param thr Use with the custom threshold method. Value(s) indicating the
#'   response time thresholds. Length must be equal to 1 or equal to the total
#'   number of items. Default is `3`.
#'
#' @param nt Use with the normative threshold method. Value(s) indicating the
#'   percentage(s) of the mean item response time to be used as thresholds. If
#'   length is equal to 1, one normative threshold is applied to all items (Wise
#'   et al., 2004). Else if length is greater than 1, multiple normative
#'   thresholds are applied to all items (Martinez & Rios, 2023). Default is
#'   `10`, for NT10.
#'
#' @param limits Use with threshold methods. A vector of length 2 indicating
#'   the minimum and maximum possible thresholds. Default is `c(0, Inf)`.
#'
#' @param min_item The minimum number of items used to identify unmotivated
#'   persons. Default is `1`.
#'
#' @returns A list is returned. If a visual inspection method is used, the list
#'   contains the following elements:
#' \item{plots}{A list containing one plot per item.}
#'
#' If a threshold method is used, the list contains the following elements:
#' \item{thr}{A vector or matrix of response time thresholds.}
#' \item{flag}{A matrix or array of flagging results.}
#' \item{rte}{A vector or matrix of response time effort, equal to 1 minus the
#'   proportion of flagged responses per person (Wise & Kong, 2005).}
#' \item{rtf}{A vector or matrix of response time fidelity, equal to 1 minus the
#'   proportion of flagged responses per item (Wise, 2006).}
#' \item{unmotivated}{The proportion of unmotivated persons.}
#'
#' @references
#' Guo, H., Rios, J. A., Haberman, S., Liu, O. L., Wang, J., & Paek, I. (2016).
#' A new procedure for detection of students' rapid guessing responses using
#' response time. *Applied Measurement in Education*, *29*(3), 173--183.
#'
#' Lee, Y.-H., & Jia, Y. (2014). Using response time to investigate students'
#' test-taking behaviors in a NAEP computer-based study. *Large-Scale
#' Assessments in Education*, *2*, Article 8.
#'
#' Ma, L., Wise, S. L., Thum, Y. M., & Kingsbury, G. (2011, April). *Detecting
#' response time threshold under the computer adaptive testing environment*
#' \[Paper presentation]. National Council of Measurement in Education, New
#' Orleans, LA, United States.
#'
#' Martinez, A. J., & Rios, J. A. (2023, April). *The impact of rapid guessing
#' on model fit and factor-analytic reliability* \[Paper presentation]. National
#' Council on Measurement in Education, Chicago, IL, United States.
#'
#' Schnipke, D. L. (1995, April). *Assessing speededness in computer-based tests
#' using item response times* \[Paper presentation]. National Council on
#' Measurement in Education, San Francisco, CA, United States.
#'
#' Wise, S. L. (2006). An investigation of the differential effort received by
#' items on a low-stakes computer-based test. *Applied Measurement in
#' Education*, *19*(2), 95--114.
#'
#' Wise, S. L., Kingsbury, G. G., Thomason, J., & Kong, X. (2004, April). *An
#' investigation of motivation filtering in a statewide achievement testing
#' program* \[Paper presentation]. National Council on Measurement in Education,
#' San Diego, CA, United States.
#'
#' Wise, S. L., & Kong, X. (2005). Response time effort: A new measure of
#' examinee motivation in computer-based tests. *Applied Measurement in
#' Education*, *18*(2), 163--183.
#'
#' Wise, S. L., & Ma, L. (2012, April). *Setting response time thresholds for a
#' CAT item pool: The normative threshold method* \[Paper presentation].
#' National Council on Measurement in Education, Vancouver, BC, Canada.
#'
#' @examples
#' # Setup for Examples 1 to 3 -------------------------------------------------
#'
#' # Settings
#' set.seed(0)     # seed for reproducibility
#' N <- 5000       # number of persons
#' n <- 40         # number of items
#'
#' # Randomly select 20% unmotivated persons
#' cv <- sample(1:N, size = N * 0.20)
#'
#' # Create vector of indicators (1 = unmotivated, 0 = motivated)
#' ind <- ifelse(1:N %in% cv, 1, 0)
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
#' t <- exp(dat$y)
#'
#' # Modify contaminated data by guessing on 20% of the items
#' for (v in cv) {
#'   ci <- sample(1:n, size = n * 0.20)
#'   x[v, ci] <- rbinom(length(ci), size = 1, prob = 0.25)
#'   t[v, ci] <- runif(length(ci), min = 1, max = 10)
#' }
#'
#' # Example 1: Visual Inspection Methods --------------------------------------
#'
#' # Detect rapid guessing using the visual inspection method
#' out <- detect_rg(
#'   method = "VI",
#'   t = t,
#'   outlier = 90
#' )
#'
#' # Detect rapid guessing using the visual inspection with proportion correct
#' # method
#' out <- detect_rg(
#'   method = "VITP",
#'   t = t,
#'   x = x,
#'   outlier = 90
#' )
#'
#' # Example 2: Threshold Methods ----------------------------------------------
#'
#' # Detect rapid guessing using the custom threshold method with a common
#' # three-second threshold
#' out <- detect_rg(
#'   method = "CT",
#'   t = t,
#'   thr = 3
#' )
#'
#' # Detect rapid guessing using the custom threshold method with 10% of the
#' # median item response time
#' out <- detect_rg(
#'   method = "CT",
#'   t = t,
#'   thr = apply(t, 2, function(i) 0.10 * median(i))
#' )
#'
#' # Detect rapid guessing using the normative threshold method with 10% of the
#' # mean item response time
#' out <- detect_rg(
#'   method = "NT",
#'   t = t,
#'   nt = 10
#' )
#'
#' # Detect rapid guessing using the normative threshold method with 5 to 35% of
#' # the mean item response time
#' out <- detect_rg(
#'   method = "NT",
#'   t = t,
#'   nt = seq(5, 35, by = 5)
#' )
#'
#' # Example 3: Visual Inspection and Threshold Methods ------------------------
#'
#' # Detect rapid guessing using the cumulative proportion correct method
#' out <- detect_rg(
#'   method = "CUMP",
#'   t = t,
#'   x = x,
#'   outlier = 90
#' )
#' @export

detect_rg <- function(method,
                      t,
                      x = NULL,
                      outlier = 100,
                      chance = 0.25,
                      thr = 3,
                      nt = 10,
                      limits = c(0, Inf),
                      min_item = 1) {

  # Checks
  method <- match.arg(
    arg = method,
    choices = c("VI", "VITP", "CT", "NT", "CUMP"),
    several.ok = FALSE
  )
  check_data(t, x)

  # Setup
  N <- nrow(t)
  n <- ncol(t)

  # Remove outliers
  if (length(outlier) == 1) {
    outlier <- rep(outlier, times = n)
  } else if (length(outlier) != n) {
    stop(
      "The length of `outlier` must be equal to 1 or equal to the total ",
      "number of items.",
      call. = FALSE
    )
  }
  if (any(outlier < 0) || any(outlier > 100)) {
    stop("`outlier` must be between 0 and 100.", call. = FALSE)
  }
  for (i in 1:n) {
    t[t[, i] > quantile(t[, i], probs = outlier[i] / 100), i] <- NA
  }

  # Apply visual inspection methods
  if (method %in% c("VI", "VITP", "CUMP")) {
    plots <- vector("list", n)
    names(plots) <- main <- paste("Item", 1:n)
    if (method == "VI") {
      for (i in 1:n) {
        t_i <- round(t[, i])
        u_i <- sort(unique(t_i))
        hist(
          t_i,
          breaks = max(u_i),
          main = main[i],
          xlab = "Response Time",
          ylab = "Frequency",
          xaxt = "n"
        )
        axis(1, at = min(u_i):max(u_i))
        plots[[i]] <- recordPlot()
      }
    } else {
      if (is.null(x)) {
        stop("`x` must be provided.", call. = FALSE)
      }
      if (length(chance) == 1) {
        chance <- rep(chance, times = n)
      } else if (length(chance) != n) {
        stop(
          "The length of `chance` must be equal to 1 or equal to the total ",
          "number of items.",
          call. = FALSE
        )
      }
      if (any(chance < 0) || any(chance > 1)) {
        stop("`chance` must be between 0 and 1.", call. = FALSE)
      }
      opar <- par(mar = par("mar"))
      par(mar = par("mar")[c(1, 2, 3, 2)])
      thr <- rep(NA, times = n)
      for (i in 1:n) {
        t_i <- round(t[, i])
        u_i <- sort(unique(t_i))
        tab <- table(t_i, x[, i])
        freq <- rowSums(tab)
        conv <- max(c(freq[1] + freq[2], freq))
        hist(
          t_i,
          breaks = max(u_i),
          main = main[i],
          xlab = "Response Time",
          ylab = "Frequency",
          xaxt = "n"
        )
        axis(1, at = min(u_i):max(u_i))
        axis(
          4,
          at = seq(0, conv, by = conv / 10),
          labels = seq(0, 1, by = 0.1),
          col = 2,
          col.axis = 2
        )
        if (method == "VITP") {
          p <- tab[, 2] / freq
          mtext("Proportion Correct", side = 4, line = 2.5, col = 2)
        } else {
          p <- cumsum(tab[, 2]) / cumsum(freq)
          mtext("Cumulative Proportion Correct", side = 4, line = 2.5, col = 2)
          if (any(p <= chance[i]) && any(p > chance[i])) {
            thr[i] <- u_i[max(which(p <= chance[i]))]
          }
        }
        lines(u_i, p * conv, col = 2, lty = 2)
        segments(
          x0 = min(u_i),
          y0 = chance[i] * conv,
          x1 = max(u_i),
          y1 = chance[i] * conv,
          col = 2
        )
        plots[[i]] <- recordPlot()
      }
      on.exit(par(opar), add = TRUE)
    }
    out <- list(plots = plots)
  } else {
    out <- NULL
  }

  # Apply threshold methods
  if (method %in% c("CT", "NT", "CUMP")) {
    if (method == "CT") {
      if (length(thr) == 1) {
        thr <- rep(thr, times = n)
      } else if (length(thr) != n) {
        stop(
          "The length of `thr` must be equal to 1 or equal to the total ",
          "number of items.",
          call. = FALSE
        )
      }
    } else if (method == "NT") {
      if (any(nt < 0) || any(nt > 100)) {
        stop("`nt` must be between 0 and 100.", call. = FALSE)
      }
      K <- length(nt)
      if (K == 1) {
        thr <- colMeans(t, na.rm = TRUE) * nt / 100
      } else {
        thr <- outer(colMeans(t, na.rm = TRUE), nt / 100, "*")
      }
    }
    if (length(limits) != 2) {
      stop("The length of `limits` must be equal to 2.", call. = FALSE)
    }
    thr[thr < limits[1]] <- limits[1]
    thr[thr > limits[2]] <- limits[2]
    if (is.matrix(thr)) {
      flag <- array(dim = c(N, n, K))
      for (k in 1:K) {
        flag[, , k] <- t(ifelse(t(t) < thr[, k], 1, 0))
      }
      rte <- 1 - apply(flag, c(1, 3), mean, na.rm = TRUE)
      rtf <- 1 - apply(flag, c(2, 3), mean, na.rm = TRUE)
      unmotivated <- apply(flag, 3, function(f)
        mean(rowSums(f, na.rm = TRUE) >= min_item))
      colnames(thr) <- dimnames(flag)[[3]] <- colnames(rte) <- colnames(rtf) <-
        names(unmotivated) <- paste0("NT", nt)
    } else {
      flag <- t(ifelse(t(t) < thr, 1, 0))
      rte <- 1 - rowMeans(flag, na.rm = TRUE)
      rtf <- 1 - colMeans(flag, na.rm = TRUE)
      unmotivated <- mean(rowSums(flag, na.rm = TRUE) >= min_item)
    }
    out <- c(
      out,
      list(
        thr = thr,
        flag = flag,
        rte = rte,
        rtf = rtf,
        unmotivated = unmotivated
      )
    )
  }

  # Output
  out
}
