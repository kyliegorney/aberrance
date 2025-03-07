#' Count the number of categories
#'
#' @noRd

count <- function(psi = NULL, ignore = NULL, x = NULL, d = NULL, r = NULL) {
  if (is.null(psi)) {
    if (!is.null(x)) {
      if (!is.null(d)) { # item scores and distractors
        apply(d, 2, max, na.rm = TRUE) + 1
      } else { # item scores
        apply(x, 2, max, na.rm = TRUE) + 1
      }
    } else { # item responses
      apply(r, 2, max, na.rm = TRUE)
    }
  } else {
    psi <- psi[, !colnames(psi) %in% ignore, drop = FALSE]
    if ("b" %in% colnames(psi)) {
      if ("lambda1" %in% colnames(psi)) { # nested logit model
        rowSums(
          !is.na(psi[, grep("lambda[1-9]", colnames(psi)), drop = FALSE])) + 1
      } else { # 3PL model
        rep(2, times = nrow(psi))
      }
    } else if ("b1" %in% colnames(psi)) { # graded response model
      rowSums(!is.na(psi[, grep("b[1-9]", colnames(psi)), drop = FALSE])) + 1
    } else if ("c0" %in% colnames(psi)) { # generalized partial credit model
      rowSums(!is.na(psi[, grep("c[0-9]", colnames(psi)), drop = FALSE]))
    } else { # nominal response model
      rowSums(!is.na(psi[, grep("lambda[1-9]", colnames(psi)), drop = FALSE]))
    }
  }
}

#' Extract substring by position
#'
#' @noRd

extract <- function(x, pos) {
  sapply(x, function(string)
    paste(unlist(strsplit(string, split = "_"))[pos], collapse = "_"),
    USE.NAMES = FALSE)
}
