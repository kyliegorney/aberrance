#' Check data
#'
#' @noRd

check_data <- function(...) {
  data <- list(...)
  data <- data[lengths(data) > 0]
  if (length(data) > 1) {
    if (any(sapply(data[-1], function(d) !identical(dim(d), dim(data[[1]]))))) {
      stop("All data must have the same dimensions.", call. = FALSE)
    }
  }
}

#' Check parameters
#'
#' @noRd

check_par <- function(data, psi, xi = NULL) {
  if ("x" %in% data) {
    if (sum(
      all(c("a", "b", "c") %in% colnames(psi)),
      all(c("a", "b1", "b2") %in% colnames(psi)),
      all(c("a", "c0", "c1") %in% colnames(psi))
    ) != 1) {
      stop("`psi` must contain columns: (a) a, b, c, (b) a, b1, b2, ..., or ",
           "(c) a, c0, c1, ...", call. = FALSE)
    } else if ("b" %in% colnames(psi)) {
      if (any(is.na(psi[, c("a", "b", "c")]))) {
        stop("`psi` must not contain missing values: a, b, c.", call. = FALSE)
      }
    } else if ("b1" %in% colnames(psi)) {
      tmp <- unique(c(
        "a",
        "b1",
        "b2",
        paste0("b", 1:sum(colnames(psi) %in% paste0("b", 1:9)))
      ))
      if (any(!tmp %in% colnames(psi))) {
        stop("`psi` must contain columns: ", paste(tmp, collapse = ", "), ".",
             call. = FALSE)
      } else if (any(is.na(psi[, c("a", "b1")]))) {
        stop("`psi` must not contain missing values: a, b1.", call. = FALSE)
      }
    } else {
      tmp <- unique(c(
        "a",
        "c0",
        "c1",
        paste0("c", 1:sum(colnames(psi) %in% paste0("c", 1:9)))
      ))
      if (any(!tmp %in% colnames(psi))) {
        stop("`psi` must contain columns: ", paste(tmp, collapse = ", "), ".",
             call. = FALSE)
      } else if (any(is.na(psi[, c("a", "c0", "c1")]))) {
        stop("`psi` must not contain missing values: a, c0, c1.", call. = FALSE)
      }
    }
  }
  if (("d" %in% data ) || ("r" %in% data)) {
    if (sum(grepl("lambda", colnames(psi))) !=
        sum(grepl("zeta", colnames(psi)))) {
      stop("`psi` must contain equal numbers of lambdas and zetas.",
           call. = FALSE)
    }
    tmp <- unique(c(
      "lambda1",
      "lambda2",
      paste0("lambda", 1:sum(colnames(psi) %in% paste0("lambda", 1:9))),
      "zeta1",
      "zeta2",
      paste0("zeta", 1:sum(colnames(psi) %in% paste0("lambda", 1:9)))
    ))
    if (any(!tmp %in% colnames(psi))) {
      stop("`psi` must contain columns: ", paste(tmp, collapse = ", "), ".",
           call. = FALSE)
    } else if (any(is.na(psi[, c("lambda1", "zeta1")]))) {
      stop("`psi` must not contain missing values: lambda1, zeta1.",
           call. = FALSE)
    }
  }
  if ("y" %in% data) {
    if (any(!c("alpha", "beta") %in% colnames(psi))) {
      stop("`psi` must contain columns: alpha, beta.", call. = FALSE)
    } else if (any(is.na(psi[, c("alpha", "beta")]))) {
      stop("`psi` must not contain missing values: alpha, beta.", call. = FALSE)
    }
  }
  if (!is.null(xi)) {
    tmp <- c(x = "theta", r = "eta", d = "eta", y = "tau")[data]
    if (any(!tmp %in% colnames(xi))) {
      stop("`xi` must contain columns: ", paste(tmp, collapse = ", "), ".",
           call. = FALSE)
    }
  }
}
