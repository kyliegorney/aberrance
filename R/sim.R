#' Simulate data
#'
#' @description Simulate data using item response theory (IRT) models.
#'
#' @param psi A matrix of item parameters.
#'
#' @param xi A matrix of person parameters.
#'
#' @details
#'
#' # Models for Item Scores
#' The **Rasch**, **2PL**, and **3PL models** (Birnbaum, 1968; Rasch, 1960) are
#' given by
#' \deqn{P(X_{vi} = 1 | \theta_v, a_i, b_i, c_i) =
#' c_i + \frac{1 - c_i}{1 + \exp\{-a_i(\theta_v - b_i)\}}.}
#' - `psi` must contain columns named `"a"`, `"b"`, and `"c"` for the item
#'   discrimination, difficulty, and pseudo-guessing parameters, respectively.
#' - `xi` must contain a column named `"theta"` for the person ability
#'   parameters.
#'
#' The **partial credit model** (PCM; Masters, 1982) and the **generalized
#' partial credit model** (GPCM; Muraki, 1992) are given by
#' \deqn{P(X_{vi} = j | \theta_v, a_i, \boldsymbol{c_i}) =
#' \frac{\exp\{\sum_{k=0}^j a_i(\theta_v - c_{ik})\}}
#' {\sum_{l=0}^{m_i} \exp\{\sum_{k=0}^l a_i(\theta_v - c_{ik})\}}.}
#' - `psi` must contain columns named `"a"` for the item discrimination
#'   parameter and `"c0"`, `"c1"`, `...`, for the item category parameters.
#' - `xi` must contain a column named `"theta"` for the person ability
#'   parameters.
#'
#' The **graded response model** (GRM; Samejima, 1969) is given by
#' \deqn{P(X_{vi} = j | \theta_v, a_i, \boldsymbol{b_i}) =
#' P(X_{vi} \ge j | \theta_v, a_i, \boldsymbol{b_i}) -
#' P(X_{vi} \ge j + 1 | \theta_v, a_i, \boldsymbol{b_i}),}
#' where
#' \deqn{P(X_{vi} \ge j | \theta_v, a_i, \boldsymbol{b_i}) = \begin{cases}
#' 1 &\text{if } j = 0, \\
#' \frac{1}{1 + \exp\{-a_i(\theta_v - b_{ij})\}}
#' &\text{if } 1 \le j \le m_i, \\
#' 0 &\text{if } j = m_i + 1.
#' \end{cases}}
#' - `psi` must contain columns named `"a"` for the item discrimination
#'   parameter and `"b1"`, `"b2"`, `...`, for the item location parameters
#'   listed in increasing order.
#' - `xi` must contain a column named `"theta"` for the person ability
#'   parameters.
#'
#' # Models for Item Distractors
#' The **nested logit model** (NLM; Bolt et al., 2012) is given by
#' \deqn{P(D_{vi} = j | \theta_v, \eta_v,
#' a_i, b_i, c_i, \boldsymbol{\lambda_i}, \boldsymbol{\zeta_i}) =
#' [1 - P(X_{vi} = 1 | \theta_v, a_i, b_i, c_i)] \times
#' P(D_{vi} = j | X_{vi} = 0, \eta_v,
#' \boldsymbol{\lambda_i}, \boldsymbol{\zeta_i}),}
#' where
#' \deqn{P(D_{vi} = j | X_{vi} = 0, \eta_v,
#' \boldsymbol{\lambda_i}, \boldsymbol{\zeta_i}) =
#' \frac{\exp(\lambda_{ij} \eta_v + \zeta_{ij})}
#' {\sum_{k=1}^{m_i-1} \exp(\lambda_{ik} \eta_v + \zeta_{ik})}.}
#' - `psi` must contain columns named `"a"`, `"b"`, and `"c"` for the item
#'   discrimination, difficulty, and pseudo-guessing parameters, respectively,
#'   `"lambda1"`, `"lambda2"`, `...`, for the item slope parameters, and
#'   `"zeta1"`, `"zeta2"`, `...`, for the item intercept parameters.
#' - `xi` must contain columns named `"theta"` and `"eta"` for the person
#'   parameters that govern response correctness and distractor selection,
#'   respectively.
#'
#' # Models for Item Responses
#' The **nominal response model** (NRM; Bock, 1972) is given by
#' \deqn{P(R_{vi} = j | \eta_v,
#' \boldsymbol{\lambda_i}, \boldsymbol{\zeta_i}) =
#' \frac{\exp(\lambda_{ij} \eta_v + \zeta_{ij})}
#' {\sum_{k=1}^{m_i} \exp(\lambda_{ik} \eta_v + \zeta_{ik})}.}
#' - `psi` must contain columns named `"lambda1"`, `"lambda2"`, `...`, for the
#'   item slope parameters and `"zeta1"`, `"zeta2"`, `...`, for the item
#'   intercept parameters. If there is a correct response category, its
#'   parameters should be listed last.
#' - `xi` must contain a column named `"eta"` for the person parameters that
#'   govern response selection.
#'
#' # Models for Item Log Response Times
#' The **lognormal model** (van der Linden, 2006) is given by
#' \deqn{f(Y_{vi} | \tau_v, \alpha_i, \beta_i) =
#' \frac{\alpha_i}{\sqrt{2 \pi}}
#' \exp\{-\frac{1}{2}[\alpha_i(Y_{vi} - (\beta_i - \tau_v))]^2\}.}
#' - `psi` must contain columns named `"alpha"` and `"beta"` for the item time
#'   discrimination and time intensity parameters, respectively.
#' - `xi` must contain a column named `"tau"` for the person speed parameters.
#'
#' @returns A list is returned. Possible elements include:
#' \item{x}{A matrix of item scores.}
#' \item{d}{A matrix of item distractors.}
#' \item{r}{A matrix of item responses.}
#' \item{y}{A matrix of item log response times.}
#'
#' @references
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an
#' examinee's ability. In F. M. Lord & M. R. Novick (Eds.), *Statistical
#' theories of mental test scores* (pp. 397--479). Addison-Wesley.
#'
#' Bock, R. D. (1972). Estimating item parameters and latent ability when
#' responses are scored in two or more nominal categories. *Psychometrika*,
#' *37*(1), 29--51.
#'
#' Bolt, D. M., Wollack, J. A., & Suh, Y. (2012). Application of a
#' multidimensional nested logit model to multiple-choice test items.
#' *Psychometrika*, *77*(2), 339--357.
#'
#' Masters, G. N. (1982). A Rasch model for partial credit scoring.
#' *Psychometrika*, *47*(2), 149--174.
#'
#' Muraki, E. (1992). A generalized partial credit model: Application of an EM
#' algorithm. *Applied Psychological Measurement*, *16*(2), 159--176.
#'
#' Rasch, G. (1960). *Probabilistic models for some intelligence and attainment
#' tests*. Danish Institute for Educational Research.
#'
#' Samejima, F. (1969). Estimation of latent ability using a response pattern of
#' graded scores. *Psychometrika*, *34*(S1), 1--97.
#'
#' van der Linden, W. J. (2006). A lognormal model for response times on test
#' items. *Journal of Educational and Behavioral Statistics*, *31*(2), 181--204.
#'
#' @examples
#' # Setup for Examples 1 to 5 -------------------------------------------------
#'
#' # Settings
#' set.seed(0)     # seed for reproducibility
#' N <- 500        # number of persons
#' n <- 40         # number of items
#'
#' # Example 1: 3PL Model and Lognormal Model ----------------------------------
#'
#' # Generate person parameters
#' xi <- MASS::mvrnorm(
#'   N,
#'   mu = c(theta = 0.00, tau = 0.00),
#'   Sigma = matrix(c(1.00, 0.25, 0.25, 0.25), ncol = 2)
#' )
#'
#' # Generate item parameters
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
#' # Simulate item scores and log response times
#' dat <- sim(psi, xi)
#' x <- dat$x
#' y <- dat$y
#'
#' # Example 2: Generalized Partial Credit Model -------------------------------
#'
#' # Generate person parameters
#' xi <- cbind(theta = rnorm(N, mean = 0.00, sd = 1.00))
#'
#' # Generate item parameters
#' psi <- cbind(
#'   a = rlnorm(n, meanlog = 0.00, sdlog = 0.25),
#'   c0 = 0,
#'   c1 = rnorm(n, mean = -1.00, sd = 0.50),
#'   c2 = rnorm(n, mean = 0.00, sd = 0.50),
#'   c3 = rnorm(n, mean = 1.00, sd = 0.50)
#' )
#'
#' # Simulate item scores
#' x <- sim(psi, xi)$x
#'
#' # Example 3: Graded Response Model ------------------------------------------
#'
#' # Generate person parameters
#' xi <- cbind(theta = rnorm(N, mean = 0.00, sd = 1.00))
#'
#' # Generate item parameters
#' psi <- cbind(
#'   a = rlnorm(n, meanlog = 0.00, sdlog = 0.25),
#'   b1 = rnorm(n, mean = -1.00, sd = 0.50),
#'   b2 = rnorm(n, mean = 0.00, sd = 0.50),
#'   b3 = rnorm(n, mean = 1.00, sd = 0.50)
#' )
#'
#' # Sort item location parameters in increasing order
#' psi[, paste0("b", 1:3)] <- t(apply(psi[, paste0("b", 1:3)], 1, sort))
#'
#' # Simulate item scores
#' x <- sim(psi, xi)$x
#'
#' # Example 4: Nested Logit Model ---------------------------------------------
#'
#' # Generate person parameters
#' xi <- MASS::mvrnorm(
#'   N,
#'   mu = c(theta = 0.00, eta = 0.00),
#'   Sigma = matrix(c(1.00, 0.80, 0.80, 1.00), ncol = 2)
#' )
#'
#' # Generate item parameters
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
#' # Simulate item scores and distractors
#' dat <- sim(psi, xi)
#' x <- dat$x
#' d <- dat$d
#'
#' # Example 5: Nominal Response Model -----------------------------------------
#'
#' # Generate person parameters
#' xi <- cbind(eta = rnorm(N, mean = 0.00, sd = 1.00))
#'
#' # Generate item parameters
#' psi <- cbind(
#'   lambda1 = rnorm(n, mean = -0.50, sd = 0.50),
#'   lambda2 = rnorm(n, mean = -0.50, sd = 0.50),
#'   lambda3 = rnorm(n, mean = -0.50, sd = 0.50),
#'   lambda4 = rnorm(n, mean = 1.50, sd = 0.50),
#'   zeta1 = rnorm(n, mean = -0.50, sd = 0.50),
#'   zeta2 = rnorm(n, mean = -0.50, sd = 0.50),
#'   zeta3 = rnorm(n, mean = -0.50, sd = 0.50),
#'   zeta4 = rnorm(n, mean = 1.50, sd = 0.50)
#' )
#'
#' # Simulate item responses
#' r <- sim(psi, xi)$r
#' @export

sim <- function(psi, xi) {

  # Setup
  N <- nrow(xi)
  n <- nrow(psi)

  # Simulate item scores, distractors, and responses
  if ("b" %in% colnames(psi)) {
    if ("lambda1" %in% colnames(psi)) {
      check_par(c("x", "d"), psi, xi)
      m <- count(psi)
      p <- irt_p(m, psi, xi)
      x <- d <- r <- matrix(nrow = N, ncol = n)
      for (v in 1:N) {
        for (i in 1:n) {
          r[v, i] <- sample(1:m[i], size = 1, prob = p[v, i, 1:m[i]])
        }
        x[v, ] <- ifelse(r[v, ] == m, 1, 0)
        d[v, ] <- ifelse(r[v, ] == m, NA, r[v, ])
      }
      r <- NULL
    } else {
      check_par("x", psi, xi)
      m <- count(psi)
      p <- irt_p(m, psi, xi)
      x <- matrix(rbinom(N * n, size = 1, prob = p[, , 2]), nrow = N, ncol = n)
      d <- r <- NULL
    }
  } else if (("b1" %in% colnames(psi)) || ("c0" %in% colnames(psi))) {
    check_par("x", psi, xi)
    m <- count(psi)
    p <- irt_p(m, psi, xi)
    x <- matrix(nrow = N, ncol = n)
    for (v in 1:N) {
      for (i in 1:n) {
        x[v, i] <- sample(0:(m[i]-1), size = 1, prob = p[v, i, 1:m[i]])
      }
    }
    d <- r <- NULL
  } else if ("lambda1" %in% colnames(psi)) {
    check_par("r", psi, xi)
    m <- count(psi)
    p <- irt_p(m, psi, xi)
    r <- matrix(nrow = N, ncol = n)
    for (v in 1:N) {
      for (i in 1:n) {
        r[v, i] <- sample(1:m[i], size = 1, prob = p[v, i, 1:m[i]])
      }
    }
    x <- d <- NULL
  } else {
    x <- d <- r <- NULL
  }

  # Simulate item log response times
  if ("beta" %in% colnames(psi)) {
    check_par("y", psi, xi)
    mu <- t(outer(psi[, "beta"], xi[, "tau"], "-"))
    sigma <- matrix(1 / psi[, "alpha"], nrow = N, ncol = n, byrow = TRUE)
    y <- matrix(rnorm(N * n, mean = mu, sd = sigma), nrow = N, ncol = n)
  } else {
    y <- NULL
  }

  # Output
  out <- list(x = x, d = d, r = r, y = y)
  out[lengths(out) > 0]
}
