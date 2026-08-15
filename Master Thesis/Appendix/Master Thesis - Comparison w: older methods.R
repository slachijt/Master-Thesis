################################################################################
##
## Comparison of four ways of testing the hot hand on a binary shooting
## sequence.  Revised version.
##
##   eps_cell   : conditional e-value of the thesis, conditions on (h, nH).
##   eps_orbit  : same model alternative, conditions on h only.
##   eps_koning : Koning (2026) log-optimal e-value for exchangeability with
##                p_hot = p_neutral^beta.
##   RR         : Ritzwoller and Romano (2022) one-sided permutation test on
##                Dhat_{n,k}, individually and pooled.
##
##
## NOTE ON A REMAINING ASYMMETRY.  eps_cell and eps_orbit condition on the
## post-prefix trials with the prefix held fixed; eps_koning and the RR
## permutation act on the full sequence of n shots.  This mirrors what each
## paper actually does and is left in place deliberately, but it should be
## stated whenever the numbers are compared.
################################################################################

rm(list = ls())
set.seed(2026)

logit     <- function(p) log(p / (1 - p))
logit_inv <- function(z) 1 / (1 + exp(-z))

## ----------------------------- GVT data -------------------------------------
gvt_data <- list(
  "101" = c(0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,1,1,1,0,1,
            0,1,1,1,0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,
            0,0,0,0,0,1,0,1,1,1,1,0,1,0,1,0,0,1,0,1,1,1,0,1,1,0,1,1,1,1,
            1,1,1,0,1,0,1,0,1,1),
  "102" = c(0,0,0,1,1,0,0,0,0,1,0,1,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,0,1,
            0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,
            0,1,0,1,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,0,
            1,0,0,0,0,0,0,0,0,1),
  "103" = c(1,1,1,1,0,0,1,0,1,1,1,0,1,0,0,0,1,1,0,1,1,1,1,1,1,0,1,0,0,1,
            0,0,0,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,
            1,1,1,0,0,1,1,1,1,1,0,0,1,1,0,0,1,1,1,0,0,0,0,1,1,1,1,1,1,0,
            0,0,0,1,1,1,1,1,0,0),
  "104" = c(0,0,1,0,1,0,0,0,1,0,0,1,1,0,0,0,0,0,1,1,0,1,1,0,1,0,1,0,1,1,
            0,0,0,0,0,1,0,1,1,1,0,1,1,0,1,1,0,0,0,0,1,0,1,0,1,1,0,0,0,0,
            1,1,1,1,0,0,1,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0),
  "105" = c(0,1,0,1,0,0,0,1,0,1,0,1,1,0,1,1,1,0,0,0,1,0,0,1,0,0,1,0,0,0,
            1,0,1,1,1,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,0,0,1,1,1,1,0,
            0,1,0,0,1,0,1,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,0,1,0,1,0,0,
            0,0,1,0,1,1,0,0,1,0),
  "106" = c(0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,1,1,1,1,1,0,1,0,0,0,1,1,1,1,1,
            1,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,
            1,1,0,0,1,1,0,1,1,1,0,0,1,0,0,0,0,0,1,1,1,1,0,1,0,0,1,1,0,1,
            0,0,1,1,1,1,1,0,0,0),
  "107" = c(0,1,1,1,1,1,1,1,0,0,0,1,1,0,0,1,1,0,1,1,1,1,1,1,0,1,0,0,1,0,
            0,1,1,1,1,1,0,1,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1,1,0,0,1,1,0,
            0,1,0,0,1,1,1,1,0,0,0,0,0,0,0),
  "108" = c(1,1,0,1,1,1,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,1,1,1,0,1,
            0,0,0,1,0,1,0,0,1,0,0,0,1,0,0,1,0,1,0,1),
  "109" = c(0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,
            0,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,
            1,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,0,1,1,0,0,0,1,0,1),
  "110" = c(0,0,0,0,0,1,0,1,0,1,1,1,1,1,1,0,1,1,0,1,0,0,0,1,1,0,1,0,1,1,
            0,0,1,1,1,1,0,1,1,1,0,0,0,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,0,1,
            0,1,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,1,0,1,0,1,1,1,1,
            1,1,1,0,1,1,1,1,1,0),
  "111" = c(1,1,0,0,1,1,0,0,1,0,0,1,1,0,0,0,1,0,1,1,0,0,1,1,0,0,0,0,1,1,
            1,0,0,1,0,1,1,1,1,0,1,0,1,0,1,1,1,1,1,1,0,0,1,0,1,0,1,0,0,1,
            1,1,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,0,
            0,0,1,0,0,0,0,0,1,0),
  "112" = c(0,1,1,1,1,0,0,1,1,0,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,
            0,1,0,0,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,1,0,0,1,1,1,0,0,
            1,1,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,
            0,1,1,1,1,1,0,1,0,1),
  "113" = c(0,0,0,0,1,0,0,1,0,1,0,1,1,0,1,1,0,1,0,1,0,1,0,1,1,1,0,0,1,1,
            1,0,0,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,1,1,1,1,1,0,1,1,1,1,0,1,
            0,1,1,1,1,1,1,0,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1,1,1,0,1,1,1,0,
            0,1,1,0,1,0,1,1,1,1),
  "114" = c(1,0,0,0,0,1,1,1,0,0,1,1,0,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,0,
            1,1,1,1,1,1,1,0,1,0,0,1,1,0,0,0,1,0,1,1,0,0,1,1,0,0,0,0,0,1,
            1,1,0,1,0,1,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,0,1,0,1,1,0,0,1,
            1,0,0,1,0,0,1,0,1,0),
  "201" = c(0,0,1,1,0,1,0,0,0,1,1,1,0,1,0,1,1,1,0,0,0,0,1,1,0,1,1,0,1,1,
            1,1,0,0,0,1,0,1,1,0,0,1,0,0,1,0,1,0,0,1,0,1,1,1,1,0,0,0,0,0,
            1,1,1,0,0,1,0,0,1,0,1,1,1,1,0,0,0,1,0,1,0,1,0,1,0,0,1,0,0,1,
            0,1,0,1,1,0,0,0,1,0),
  "202" = c(0,0,0,1,0,1,0,0,0,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,0,0,0,1,0,
            0,0,1,0,1,1,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,
            0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,1,0,0,1,1,1,0,0,0,0,1,1,1,0,
            0,0,0,0,0,1,0,1,0,0),
  "203" = c(0,0,0,1,0,0,0,1,0,0,0,1,0,1,1,1,0,1,1,0,0,1,1,1,0,0,1,1,1,1,
            0,0,0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,
            0,0,0,1,1,1,1,1,1,0,0,0,1,0,0,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,
            0,0,0,0,0,0,0,0,1,0),
  "204" = c(0,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,
            0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,1,0,0,1,
            0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,
            0,0,1,1,1,1,0,0,0,1),
  "205" = c(1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,1,0,0,
            1,0,1,0,0,1,1,1,0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,0,1,
            0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1,1,1,1,
            0,1,0,0,1,0,0,1,1,1),
  "206" = c(1,0,1,1,1,0,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,1,0,1,1,0,1,0,1,0,
            0,1,0,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,0,1,1,1,0,1,0,0,1,
            0,0,0,1,1,1,0,1,0,1,1,1,1,0,0,1,1,0,1,0,0,0,0,1,0,1,0,0,1,0,
            1,0,0,1,1,0,0,1,0,1),
  "207" = c(0,1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,1,1,1,
            1,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,0,1,1,1,0,1,0,0,0,1,0,1,0,0,
            0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,1,0,1,0,
            1,1,1,1,0,0,0,0,0,0),
  "208" = c(0,0,0,1,1,0,1,0,1,0,1,0,1,1,0,1,0,0,1,1,1,1,1,0,0,1,0,1,1,0,
            0,0,1,0,1,0,0,1,0,1,1,0,1,0,1,0,1,0,0,0,1,0,0,0,1,1,0,1,0,0,
            0,0,0,0,1,1,1,1,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1,1,1,0,0,1,1,0,
            0,1,0,1,1,1,1,1,0,1),
  "209" = c(1,1,0,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0,1,0,1,1,0,0,1,1,0,
            0,1,0,1,0,0,1,1,1,1,1,0,1,0,0,1,0,0,1,0,0,0,0,1,1,1,0,0,1,0,
            0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,1,0,0,0,1,0,1,0,1,
            0,1,0,1,1,0,1,0,1,1),
  "210" = c(0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,1,0,0,1,1,1,1,1,1,
            1,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,1,0,
            1,0,1,1,1,0,0,1,1,0,1,0,1,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,1,1,
            1,1,1,1,1,1,0,1,0,1),
  "211" = c(0,1,0,1,1,0,1,1,0,0,0,0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,0,0,1,
            1,1,0,0,0,0,0,1,0,0,1,0,1,0,1,1,1,0,1,1,1,1,0,1,1,1,0,0,1,0,
            0,0,1,0,0,0,0,1,1,0,1,1,1,0,1,0,0,1,1,1,0,1,0,0,1,1,0,0,1,1,
            0,1,0,0,1,1,1,1,1,1),
  "212" = c(0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,
            1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,
            0,0,0,0,0,1,0,0,0,1,1,0,1,1,0,0,0,1,1,0,0,1,0,0,0,0,0,0,1,1,
            0,1,0,0,0,1,0,0,1,0)
)

################################################################################
## 1. EXACT DYNAMIC PROGRAMME FOR THE CELL PARTITION FUNCTION
##    Z_delta(h, nH) = sum over x in X(h, nH) of exp(delta * hH(x)).
##    Returns logZ[nH+1, h+1], with -Inf on empty cells.  Unchanged: verified
##    correct against a direct enumeration for n <= 16.
################################################################################

dp_logZ <- function(m, k, r0, delta) {
  ed  <- exp(delta)
  W   <- array(0, dim = c(k + 1L, m + 1L, m + 1L))
  W[r0 + 1L, 1L, 1L] <- 1
  off <- 0
  for (step in seq_len(m)) {
    Wn <- array(0, dim = c(k + 1L, m + 1L, m + 1L))
    for (r in 0:k) {
      A <- W[r + 1L, , , drop = TRUE]
      if (max(A) == 0) next
      hot <- (r == k)
      rp  <- min(r + 1L, k)
      if (hot) {
        Wn[1L, 2:(m + 1L), ] <- Wn[1L, 2:(m + 1L), ] + A[1:m, , drop = FALSE]
        Wn[rp + 1L, 2:(m + 1L), 2:(m + 1L)] <-
          Wn[rp + 1L, 2:(m + 1L), 2:(m + 1L)] + ed * A[1:m, 1:m, drop = FALSE]
      } else {
        Wn[1L, , ] <- Wn[1L, , ] + A
        Wn[rp + 1L, , 2:(m + 1L)] <- Wn[rp + 1L, , 2:(m + 1L)] +
          A[, 1:m, drop = FALSE]
      }
    }
    mx <- max(Wn)
    if (mx > 0) { Wn <- Wn / mx; off <- off + log(mx) }
    W <- Wn
  }
  Zmat <- apply(W, c(2L, 3L), sum)
  ifelse(Zmat > 0, log(Zmat) + off, -Inf)
}

prefix_run <- function(prefix, k) {
  r <- 0L
  for (v in rev(prefix)) { if (v == 1L) r <- r + 1L else break }
  min(r, k)
}

.dp_cache <- new.env(parent = emptyenv())
get_logZ <- function(m, k, r0, delta) {
  key <- sprintf("%d_%d_%d_%.10f", m, k, r0, delta)
  if (!exists(key, envir = .dp_cache)) {
    assign(key, dp_logZ(m, k, r0, delta), envir = .dp_cache)
  }
  get(key, envir = .dp_cache)
}

################################################################################
## 2. COUNT STATISTICS AND THE THREE E-VALUES
################################################################################

counts_one <- function(x, k) {
  n    <- length(x)
  hot  <- vapply((k + 1L):n, function(t) all(x[(t - k):(t - 1L)] == 1L),
                 logical(1))
  post <- x[(k + 1L):n]
  list(m = n - k, h = sum(post), nH = sum(hot), hH = sum(post[hot]),
       nN = sum(!hot), hN = sum(post[!hot]))
}

log_eps_cell <- function(x, k, delta1, delta0 = 0) {
  s  <- counts_one(x, k)
  if (s$nH == 0L) return(0)
  r0 <- prefix_run(x[1:k], k)
  L1 <- get_logZ(s$m, k, r0, delta1)[s$nH + 1L, s$h + 1L]
  L0 <- get_logZ(s$m, k, r0, delta0)[s$nH + 1L, s$h + 1L]
  (delta1 - delta0) * s$hH + L0 - L1
}

## eta = NULL plugs in the null MLE logit(h/m), which is orbit-measurable and
## therefore costs nothing in validity.  Pass a number to use a fixed eta.
log_eps_orbit <- function(x, k, delta1, eta = NULL) {
  s  <- counts_one(x, k)
  if (s$h == 0L || s$h == s$m) return(0)
  if (is.null(eta)) eta <- logit(s$h / s$m)
  r0 <- prefix_run(x[1:k], k)
  col <- get_logZ(s$m, k, r0, delta1)[, s$h + 1L]
  c1  <- log1p(exp(eta + delta1)) - log1p(exp(eta))
  nHv <- 0:s$m
  ok  <- is.finite(col)
  a   <- col[ok] - c1 * nHv[ok]
  M   <- max(a)
  logS1 <- M + log(sum(exp(a - M)))
  delta1 * s$hH - c1 * s$nH - logS1 + lchoose(s$m, s$h)
}

## Koning (2026).  Verified against his worked example: the sequence with
## counts (4 hits, 2 misses) written 111010, k = 2, beta = 0.9, returns
## q = 0.0673029 and eps = 1.0095446, matching his 0.0673 and 1.0095.
log_eps_koning <- function(x, k, beta) {
  n <- length(x); nh <- sum(x)
  if (nh == 0L || nh == n) return(0)
  rem_h <- nh; rem_s <- n; lp <- 0
  for (t in seq_len(n)) {
    p_n <- rem_h / rem_s
    hot <- (t > k) && all(x[(t - k):(t - 1L)] == 1L)
    p_t <- if (hot && p_n > 0 && p_n < 1) p_n^beta else p_n
    lp  <- lp + if (x[t] == 1L) log(p_t) else log(1 - p_t)
    rem_s <- rem_s - 1L; rem_h <- rem_h - x[t]
  }
  lp + lchoose(n, nh)
}

## Matching rule for beta.  Kept for reference, but see the grid search below:
## it is badly misspecified for large delta because Koning's p_neutral is the
## running remaining-hit rate, and because his alternative carries no
## -c(eta, delta) * nH offset.
beta_from_delta <- function(theta_N, delta) {
  if (delta == 0) return(1)
  theta_H <- logit_inv(logit(theta_N) + delta)
  log(theta_H) / log(theta_N)
}

################################################################################
## 3. RITZWOLLER AND ROMANO: Dhat_{n,k}, INDIVIDUAL AND POOLED
################################################################################

## Vectorised Dhat over a matrix of sequences, one sequence per row.
D_mat <- function(X, k) {
  n  <- ncol(X)
  Cs <- cbind(0, t(apply(X, 1L, cumsum)))
  S  <- Cs[, (k + 1L):n, drop = FALSE] - Cs[, 1:(n - k), drop = FALSE]
  hot  <- S == k
  cold <- S == 0
  Y  <- X[, (k + 1L):n, drop = FALSE]
  nh <- rowSums(hot); nc <- rowSums(cold)
  out <- numeric(nrow(X))
  ok  <- nh > 0 & nc > 0
  out[ok] <- rowSums(Y * hot)[ok] / nh[ok] - rowSums(Y * cold)[ok] / nc[ok]
  out
}

D_stat <- function(x, k) D_mat(matrix(x, nrow = 1L), k)

## M uniform reorderings of one sequence, as an M x n matrix.
perm_mat <- function(x, M) {
  n   <- length(x)
  idx <- t(replicate(M, sample.int(n)))
  matrix(x[as.vector(idx)], nrow = M)
}

## Individual one-sided permutation p-value with the (M+1) correction, and the
## full permutation draw, which the pooled test reuses.
rr_one <- function(x, k, M) {
  obs  <- D_stat(x, k)
  perm <- D_mat(perm_mat(x, M), k)
  list(p = (1 + sum(perm >= obs)) / (M + 1), obs = obs, perm = perm)
}

## Pooled tests of the joint null over a block of shooters.  Xlist is a list of
## sequences.  Returns RR's stratified permutation p-value for Dbar_k (their
## Section 5.3) and Fisher's combination of the individual p-values.  The two
## differ in what they need: Dbar needs the raw sequences re-permuted jointly,
## Fisher needs only the 26 marginal p-values, which is the fair structural
## counterpart of multiplying 26 e-values.
rr_pooled <- function(Xlist, k, M) {
  fits <- lapply(Xlist, rr_one, k = k, M = M)
  pind <- vapply(fits, function(f) f$p, numeric(1))
  obs  <- mean(vapply(fits, function(f) f$obs, numeric(1)))
  PD   <- do.call(rbind, lapply(fits, function(f) f$perm))   # (#shooters) x M
  Dbar <- colMeans(PD)
  list(p_ind    = pind,
       p_dbar   = (1 + sum(Dbar >= obs)) / (M + 1),
       p_fisher = pchisq(-2 * sum(log(pind)), df = 2 * length(pind),
                         lower.tail = FALSE))
}

################################################################################
## 4. p-VALUE -> e-VALUE CALIBRATION
##    f_kappa(p) = kappa * p^(kappa-1), kappa in (0,1), is admissible (Vovk
##    2021).  Every calibration is lossy, so a calibrated p-value understates
##    what the test can do.  kappa must be fixed in advance; report a grid.
################################################################################

log_calib <- function(p, kappa = 0.4) log(kappa) + (kappa - 1) * log(p)

################################################################################
## 5. EXACT QUANTITIES UNDER THE MODEL (Part A)
################################################################################

exact_joint <- function(m, k, r0, eta, delta) {
  L  <- get_logZ(m, k, r0, delta)
  nH <- matrix(0:m, nrow = m + 1L, ncol = m + 1L)
  h  <- matrix(0:m, nrow = m + 1L, ncol = m + 1L, byrow = TRUE)
  lg <- L + eta * h - nH * log1p(exp(eta + delta)) - (m - nH) * log1p(exp(eta))
  lg[!is.finite(L)] <- -Inf
  Q <- exp(lg - max(lg[is.finite(lg)]))
  Q / sum(Q)
}

exact_EhH <- function(m, k, r0, delta, eps = 1e-3) {
  Lp <- get_logZ(m, k, r0, delta + eps)
  Lm <- get_logZ(m, k, r0, delta - eps)
  out <- (Lp - Lm) / (2 * eps)
  out[!is.finite(Lp) | !is.finite(Lm)] <- 0
  out
}

## eta_orbit = "true" reproduces the identity of Remark [exact price] exactly,
## because that remark takes the alternative to be the true conditional law.
## eta_orbit = "mle" is the feasible version used in the simulation, and then
## gap_check is strictly below divergence.
exact_price <- function(m, k, r0, eta, delta, eta_orbit = c("true", "mle")) {
  eta_orbit <- match.arg(eta_orbit)
  Q    <- exact_joint(m, k, r0, eta, delta)
  L1   <- get_logZ(m, k, r0, delta)
  L0   <- get_logZ(m, k, r0, 0)
  EhH  <- exact_EhH(m, k, r0, delta)
  nHv  <- 0:m
  
  term_cell <- delta * EhH + L0 - L1
  term_cell[!is.finite(L1) | !is.finite(L0)] <- 0
  E_cell <- sum(Q * term_cell)
  
  E_orb <- 0; div <- 0
  for (hi in 0:m) {
    col <- L1[, hi + 1L]; ok <- is.finite(col)
    if (!any(ok)) next
    eta_h <- if (eta_orbit == "true") eta else
      if (hi == 0 || hi == m) eta else logit(hi / m)
    c1 <- log1p(exp(eta_h + delta)) - log1p(exp(eta_h))
    a  <- col[ok] - c1 * nHv[ok]; M <- max(a)
    logS1 <- M + log(sum(exp(a - M)))
    t_orb <- delta * EhH[, hi + 1L] - c1 * nHv - logS1 + lchoose(m, hi)
    t_orb[!ok] <- 0
    E_orb <- E_orb + sum(Q[, hi + 1L] * t_orb)
    
    qh <- sum(Q[, hi + 1L])
    if (qh > 0) {
      qn <- Q[, hi + 1L] / qh
      p0 <- exp(L0[, hi + 1L] - lchoose(m, hi))
      p0[!is.finite(L0[, hi + 1L])] <- 0
      pos <- qn > 0
      div <- div + qh * sum(qn[pos] * log(qn[pos] / p0[pos]))
    }
  }
  list(E_logeps_cell = E_cell, E_logeps_orbit = E_orb,
       divergence = div, gap_check = E_orb - E_cell)
}

################################################################################
## 6. DATA GENERATING PROCESS
################################################################################

simulate_seq <- function(n, k, theta_N, delta) {
  eta <- logit(theta_N)
  x <- integer(n)
  x[1:k] <- as.integer(runif(k) < theta_N)
  for (t in (k + 1L):n) {
    hot  <- all(x[(t - k):(t - 1L)] == 1L)
    x[t] <- as.integer(runif(1) < logit_inv(eta + delta * hot))
  }
  x
}

################################################################################
## 7. PART A: exact price of the finer conditioning
##
##    With eta_orbit = "true" the gap equals the divergence to machine
##    precision, which is the identity of Remark [exact price].  This is the
##    version the thesis table reports, so the table note should say that the
##    orbit arm uses the true eta, not the null MLE.
################################################################################

cat("\n=== Part A: exact price of conditioning on nH as well as h ===\n")
cat(sprintf("%4s %4s %7s %7s %11s %11s %11s %11s\n",
            "n", "k", "thetaN", "delta", "E[log e_h]", "E[log e_c]",
            "gap", "divergence"))
gridA <- expand.grid(n = c(50, 100), k = 1:3, theta_N = c(0.4, 0.5),
                     delta = c(0.15, 0.30, 0.60))
for (i in seq_len(nrow(gridA))) {
  g  <- gridA[i, ]
  e  <- exact_price(g$n - g$k, g$k, g$k, logit(g$theta_N), g$delta,
                    eta_orbit = "true")
  cat(sprintf("%4d %4d %7.2f %7.2f %11.5f %11.5f %11.5f %11.5f\n",
              g$n, g$k, g$theta_N, g$delta,
              e$E_logeps_orbit, e$E_logeps_cell, e$gap_check, e$divergence))
}

################################################################################
## 8. PART B: comparison of the four methods, by streak length
##
##    Design.  For each delta we draw B independent blocks of npool shooters.
##    Individual quantities are computed on all B*npool sequences; pooled
##    quantities on the B blocks.  The alternative is delta_1 = delta for
##    delta > 0 and delta_1 = delta1_null in the null row, so that the null row
##    tests a genuine e-value rather than the constant one.
################################################################################

run_partB <- function(k, delta_grid, B = 200L, npool = 26L, M_perm = 199L,
                      n_sim = 100L, theta_N = 0.50, alpha = 0.05,
                      delta1_null = 0.40, kappa_cal = 0.40,
                      beta_grid = seq(0.50, 1.00, by = 0.05)) {
  
  thr <- log(1 / alpha)
  out <- data.frame(k = k, delta = delta_grid, delta1 = NA, beta_m = NA,
                    beta_o = NA,
                    elog_cell = NA, elog_orbit = NA, elog_kon = NA,
                    elog_konO = NA, elog_RRcal = NA,
                    mean_cell = NA, mean_orbit = NA, mean_kon = NA,
                    rej_cell = NA, rej_orbit = NA, rej_kon = NA,
                    rej_konO = NA, rej_RR = NA,
                    pool_cell = NA, pool_orbit = NA, pool_kon = NA,
                    pool_konO = NA, pool_RRdbar = NA, pool_RRfisher = NA)
  
  for (di in seq_along(delta_grid)) {
    d  <- delta_grid[di]
    d1 <- if (d > 0) d else delta1_null
    bm <- beta_from_delta(theta_N, d1)
    
    R  <- B * npool
    X  <- vector("list", R)
    for (i in seq_len(R)) X[[i]] <- simulate_seq(n_sim, k, theta_N, d)
    
    lc <- vapply(X, log_eps_cell,  numeric(1), k = k, delta1 = d1)
    lo <- vapply(X, log_eps_orbit, numeric(1), k = k, delta1 = d1)
    lk <- vapply(X, log_eps_koning, numeric(1), k = k, beta = bm)
    
    ## give Koning his best beta on the same draws
    bmean <- vapply(beta_grid, function(b)
      mean(vapply(X, log_eps_koning, numeric(1), k = k, beta = b)), numeric(1))
    bo  <- beta_grid[which.max(bmean)]
    lkb <- vapply(X, log_eps_koning, numeric(1), k = k, beta = bo)
    
    ## RR, block by block
    p_ind <- numeric(R); p_db <- numeric(B); p_fi <- numeric(B)
    for (b in seq_len(B)) {
      idx <- ((b - 1L) * npool + 1L):(b * npool)
      rr  <- rr_pooled(X[idx], k, M_perm)
      p_ind[idx] <- rr$p_ind
      p_db[b]    <- rr$p_dbar
      p_fi[b]    <- rr$p_fisher
    }
    
    pool <- function(v) mean(rowSums(matrix(v, nrow = B, byrow = TRUE)) >= thr)
    
    out$delta1[di]    <- d1
    out$beta_m[di]    <- bm
    out$beta_o[di]    <- bo
    out$elog_cell[di] <- mean(lc);  out$elog_orbit[di] <- mean(lo)
    out$elog_kon[di]  <- mean(lk);  out$elog_konO[di]  <- mean(lkb)
    out$elog_RRcal[di] <- mean(log_calib(p_ind, kappa_cal))
    ## validity: the ARITHMETIC mean must be <= 1, not the geometric mean
    out$mean_cell[di]  <- mean(exp(lc))
    out$mean_orbit[di] <- mean(exp(lo))
    out$mean_kon[di]   <- mean(exp(lk))
    out$rej_cell[di]  <- mean(lc  >= thr)
    out$rej_orbit[di] <- mean(lo  >= thr)
    out$rej_kon[di]   <- mean(lk  >= thr)
    out$rej_konO[di]  <- mean(lkb >= thr)
    out$rej_RR[di]    <- mean(p_ind <= alpha)          ## the fix
    out$pool_cell[di]  <- pool(lc)
    out$pool_orbit[di] <- pool(lo)
    out$pool_kon[di]   <- pool(lk)
    out$pool_konO[di]  <- pool(lkb)
    out$pool_RRdbar[di]   <- mean(p_db <= alpha)
    out$pool_RRfisher[di] <- mean(p_fi <= alpha)
    
    cat(sprintf("  k=%d delta=%.2f (d1=%.2f, beta_m=%.3f, beta_o=%.2f)\n",
                k, d, d1, bm, bo))
  }
  out
}

cat("\n=== Part B: simulation comparison ===\n")
delta_grid <- seq(0, 0.8, by = 0.20)
resB <- do.call(rbind, lapply(c(1L, 2L), run_partB, delta_grid = delta_grid,
                              B = 150L, M_perm = 199L))

cat("\n--- validity at delta = 0 (mean e-value should be <= 1) ---\n")
z <- resB[resB$delta == 0, ]
print(round(z[, c("k", "delta1", "mean_cell", "mean_orbit", "mean_kon",
                  "rej_cell", "rej_orbit", "rej_kon", "rej_RR",
                  "pool_RRdbar", "pool_RRfisher")], 4), row.names = FALSE)

cat("\n--- e-power ---\n")
print(round(resB[, c("k", "delta", "elog_cell", "elog_orbit", "elog_kon",
                     "elog_konO", "elog_RRcal")], 4), row.names = FALSE)
cat("\n--- rejection rate at alpha = 0.05, single sequence ---\n")
print(round(resB[, c("k", "delta", "rej_cell", "rej_orbit", "rej_kon",
                     "rej_konO", "rej_RR")], 4), row.names = FALSE)
cat("\n--- pooled rejection rate at alpha = 0.05, 26 shooters ---\n")
print(round(resB[, c("k", "delta", "pool_cell", "pool_orbit", "pool_kon",
                     "pool_konO", "pool_RRdbar", "pool_RRfisher")], 4),
      row.names = FALSE)

################################################################################
## 9. PART C: the GVT data under the four methods
##    gvt_data.rds should hold a named list of 26 integer vectors.
################################################################################

if (file.exists("gvt_data.rds")) {
  gvt_data <- readRDS("gvt_data.rds")
} else {
  gvt_data <- NULL
  cat("\n[Part C skipped: put the 26 GVT sequences in gvt_data.rds]\n")
}

if (!is.null(gvt_data)) {
  cat("\n=== Part C: GVT data ===\n")
  ids    <- names(gvt_data)
  ks     <- c(1L, 2L, 3L)
  d1     <- 0.22
  bt     <- beta_from_delta(0.50, d1)
  M_gvt  <- 9999L
  kappas <- c(0.2, 0.4, 0.6)
  cat(sprintf("delta_1 = %.2f, matched beta = %.4f, M = %d\n", d1, bt, M_gvt))
  
  for (k in ks) {
    ec <- eo <- ek <- pr <- numeric(length(ids))
    Xk <- lapply(gvt_data, as.integer)
    for (j in seq_along(ids)) {
      x <- Xk[[j]]
      ec[j] <- exp(log_eps_cell(x, k, d1))
      eo[j] <- exp(log_eps_orbit(x, k, d1))
      ek[j] <- exp(log_eps_koning(x, k, bt))
    }
    rr <- rr_pooled(Xk, k, M_gvt)
    pr <- rr$p_ind
    
    tab <- data.frame(shooter = ids, cell = ec, orbit = eo, koning = ek,
                      p_RR = pr)
    for (kap in kappas)
      tab[[sprintf("RRcal_%.1f", kap)]] <- exp(log_calib(pr, kap))
    cat(sprintf("\n-- k = %d --\n", k))
    print(round(tab[, -1], 4), row.names = FALSE)
    
    cat(sprintf("  products: cell = %.6f  orbit = %.6f  koning = %.6f\n",
                prod(ec), prod(eo), prod(ek)))
    for (kap in kappas)
      cat(sprintf("  RR calibrated product (kappa = %.1f) = %.6f\n",
                  kap, prod(exp(log_calib(pr, kap)))))
    cat(sprintf("  RR Fisher p = %.5f    RR Dbar_%d stratified permutation p = %.5f\n",
                rr$p_fisher, k, rr$p_dbar))
    cat(sprintf("  RR individual rejections at 0.05: %d of %d; Holm: %d\n",
                sum(pr <= 0.05), length(pr), sum(p.adjust(pr, "holm") <= 0.05)))
  }
}
