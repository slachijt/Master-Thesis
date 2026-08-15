rm(list = ls())
cat("\014")

################################################################################
## Conditional e-values on the GVT (Gilovich, Vallone, Tversky 1985) data --
## Every Z_delta(h,n_H) is computed by
## the recursion of the appendix, one cell row extracted per shooter.
##
## Caching: the exact table only depends on (n, k, r0), where r0 is the
## trailing run length in the shooter's own prefix x_{1:k}, capped at k -- NOT
## on the full prefix identity and NOT on the shooter's own (h, n_H). So all
## shooters sharing (n, k, r0) reuse one built table, cutting 26 x 3 = 78
## potential builds down to a handful (18 for the three k in this script).
##
## Two blocks, unchanged in spirit from the MC version:
##  (A) "Hot hand" mixture e-values against H0: delta <= delta0, using a
##      truncated-normal mixture over delta1 > delta0.
##  (B) "Sceptic" point e-values against H0: delta >= delta0, using the FIXED
##      alternative delta1 = 0, for delta0 in {0.6, 0.8, 1.0}.
################################################################################

## ----------------------------- CONFIG (edit me) -----------------------------
k_vec <- c(1, 2, 3)

## --- Block A: mixture e-values against delta <= delta0 -----------------------
delta0_vec_A <- c(0, 0.1, 0.2)
prior_sd     <- 0.1
prior_mu_vec <- c(0.1, 0.2, 0.3)
prior_hi     <- 3.0
prior_ngrid  <- 200

## --- Block B: sceptic's point e-values against delta >= delta0 ---------------
delta0_vec_B <- c(0.6, 0.8, 1.0)
delta1_B     <- 0

out_tex <- "evalue_table_exact.tex"


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

## ----------------------------------------------------------------------------
trail_run <- function(prefix, k) {
  r <- 0L
  for (j in k:1) if (prefix[j] == 1L) r <- r + 1L else break
  r
}

## N_t(r, h, nH, hH): sequences of the first t trials with trailing hit run
## min(r,k) and post-prefix counts (h, nH, hH). Exact for every k (verified
## against brute-force enumeration up to k = 4).
count_table <- function(n, k, prefix) {
  m  <- n - k
  r0 <- trail_run(prefix, k)
  N <- array(0, dim = c(k + 1L, m + 1L, m + 1L, m + 1L))
  N[r0 + 1L, 1L, 1L, 1L] <- 1
  for (t in (k + 1L):n) {
    u <- t - k
    B <- array(0, dim = dim(N))
    for (r in 0:k) {
      src <- N[r + 1L, 1:u, 1:u, 1:u, drop = FALSE]
      if (all(src == 0)) next
      src <- array(src, dim = c(u, u, u))
      rp  <- min(r + 1L, k)
      if (r == k) {                                   # trial t is Hot
        B[k + 1L, 2:(u + 1L), 2:(u + 1L), 2:(u + 1L)] <-
          B[k + 1L, 2:(u + 1L), 2:(u + 1L), 2:(u + 1L)] + src   # hit
        B[1L, 1:u, 2:(u + 1L), 1:u] <-
          B[1L, 1:u, 2:(u + 1L), 1:u] + src                     # miss
      } else {                                        # trial t is Neutral
        B[rp + 1L, 2:(u + 1L), 1:u, 1:u] <-
          B[rp + 1L, 2:(u + 1L), 1:u, 1:u] + src                # hit
        B[1L, 1:u, 1:u, 1:u] <-
          B[1L, 1:u, 1:u, 1:u] + src                            # miss
      }
    }
    N <- B
  }
  apply(N, c(2, 3, 4), sum)          # N[h+1, nH+1, hH+1]
}

## Cell table: one row per non-empty (h, nH) cell, holding log N(h,nH,hH)
## across hH = 0,...,m. idx[h+1,nH+1] gives the row for that cell.
build_cells <- function(n, k, prefix) {
  m <- n - k
  N <- count_table(n, k, prefix)
  mass <- apply(N, c(1, 2), sum)
  sel  <- which(mass > 0, arr.ind = TRUE)
  W <- matrix(0, nrow(sel), m + 1L)
  for (i in seq_len(nrow(sel))) W[i, ] <- N[sel[i, 1], sel[i, 2], ]
  logW <- log(W)
  list(n = n, k = k, m = m,
       h = as.numeric(sel[, 1] - 1L), nH = as.numeric(sel[, 2] - 1L),
       v = 0:m, logW = logW, ncell = nrow(sel),
       idx = { M <- matrix(NA_integer_, m + 1L, m + 1L)
       M[sel] <- seq_len(nrow(sel)); M })
}

## Numerically stable log( sum_v exp(logW_row + delta*v) )
logsumexp_row <- function(logW_row, v, delta) {
  a  <- logW_row + delta * v
  mx <- max(a[is.finite(a)])
  mx + log(sum(exp(a - mx)))
}

## Exact log e-value against the point null delta = delta0, for ONE cell row.
log_evalue_point_exact <- function(d0, d1, logW_row, v, hH_obs) {
  logsumexp_row(logW_row, v, d0) - logsumexp_row(logW_row, v, d1) + (d1 - d0) * hH_obs
}

## Exact log mixture e-value: log( sum_j weights[j] * eps_{d0, grid[j]}(x) )
log_evalue_mixture_exact <- function(d0, grid, weights, logW_row, v, hH_obs) {
  log_eps <- vapply(grid, function(d1)
    log_evalue_point_exact(d0, d1, logW_row, v, hH_obs), numeric(1))
  M <- max(log_eps)
  M + log(sum(weights * exp(log_eps - M)))
}
## ----------------------------- helper functions ------------------------------
counts_one <- function(x, k) {
  n    <- length(x)
  hot  <- vapply((k + 1):n, function(t) all(x[(t - k):(t - 1)] == 1), logical(1))
  post <- x[(k + 1):n]
  list(h = sum(post), nH = sum(hot), hH = sum(post[hot]))
}

tn_grid_weights <- function(mu, sd, lo, hi, ngrid) {
  grid    <- seq(lo, hi, length.out = ngrid)
  weights <- dnorm(grid, mu, sd)
  weights <- weights / sum(weights)
  list(grid = grid, weights = weights)
}

shooter_ids <- names(gvt_data)
J           <- length(shooter_ids)

## ----------------------------- shared cell computation ------------------------
## For each k, build one exact table per DISTINCT (n, r0) and reuse it across
## every shooter that shares it. This is the exact replacement for the old
## "cells_by_k" MC cache -- same shape, no sampling.
cat("=== building exact recursion tables (cached by n, k, r0) ===\n")
cells_by_k <- vector("list", length(k_vec))
names(cells_by_k) <- paste0("k", k_vec)

for (k in k_vec) {
  ## observed counts per shooter, and the cache key (n, r0) for each
  obs_list <- lapply(gvt_data, counts_one, k = k)
  r0_of <- function(x) trail_run(x[1:k], k)
  keys  <- sapply(gvt_data, function(x) paste0(length(x), "_", r0_of(x)))
  
  table_cache <- new.env(hash = TRUE)
  cells <- vector("list", J)
  names(cells) <- shooter_ids
  
  for (j in seq_len(J)) {
    x    <- gvt_data[[shooter_ids[j]]]
    n_j  <- length(x)
    key  <- keys[j]
    if (!exists(key, envir = table_cache, inherits = FALSE)) {
      t0 <- proc.time()[3]
      CL <- build_cells(n_j, k, prefix = x[1:k])
      assign(key, CL, envir = table_cache)
      cat(sprintf("  k=%d built table for key=%s (n=%d): %.2fs, %d cells\n",
                  k, key, n_j, proc.time()[3] - t0, CL$ncell))
    }
    CL  <- get(key, envir = table_cache, inherits = FALSE)
    obs <- obs_list[[j]]
    row <- CL$idx[obs$h + 1L, obs$nH + 1L]
    if (is.na(row)) stop(sprintf("shooter %s, k=%d: observed cell (h=%d,nH=%d) not found in table",
                                 shooter_ids[j], k, obs$h, obs$nH))
    cells[[j]] <- list(hH_obs = obs$hH, logW_row = CL$logW[row, ], v = CL$v)
  }
  cells_by_k[[paste0("k", k)]] <- cells
  cat(sprintf("k=%d done, %d distinct tables built for %d shooters\n\n",
              k, length(ls(table_cache)), J))
}

## ======================= BLOCK A: hot-hand mixture e-values ===================
## H0: delta <= delta0, tested with delta1 ~ TN(mu, sd) on (delta0, prior_hi).
log_ev_grid_A <- list()

for (k in k_vec) {
  cells <- cells_by_k[[paste0("k", k)]]
  for (idx in seq_along(delta0_vec_A)) {
    delta0   <- delta0_vec_A[idx]
    prior_mu <- prior_mu_vec[idx]
    tn <- tn_grid_weights(prior_mu, prior_sd, lo = delta0, hi = prior_hi, ngrid = prior_ngrid)
    
    log_ev <- numeric(J)
    for (j in seq_len(J)) {
      cl <- cells[[j]]
      log_ev[j] <- log_evalue_mixture_exact(delta0, tn$grid, tn$weights,
                                            cl$logW_row, cl$v, cl$hH_obs)
    }
    
    key <- paste0("k", k, "_d0_", delta0)
    log_ev_grid_A[[key]] <- log_ev
    
    cat(sprintf("[A] done: k=%d, delta0=%.2f (prior TN(%.1f,%.1f) on (%.2f,%.1f))\n",
                k, delta0, prior_mu, prior_sd, delta0, prior_hi))
  }
}

## ======================= BLOCK B: sceptic's point e-values ====================
log_ev_grid_B <- list()

for (k in k_vec) {
  cells <- cells_by_k[[paste0("k", k)]]
  for (delta0 in delta0_vec_B) {
    stopifnot(delta1_B < delta0)
    
    log_ev <- numeric(J)
    for (j in seq_len(J)) {
      cl <- cells[[j]]
      log_ev[j] <- log_evalue_point_exact(delta0, delta1_B, cl$logW_row, cl$v, cl$hH_obs)
    }
    
    key <- paste0("k", k, "_d0_", delta0)
    log_ev_grid_B[[key]] <- log_ev
    
    cat(sprintf("[B] done: k=%d, delta0=%.2f, delta1=%.2f (sceptic's test, H0: delta>=%.2f)\n",
                k, delta0, delta1_B, delta0))
  }
}

## ----------------------------- assemble panels --------------------------------
fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

build_panels <- function(log_ev_grid, delta0_vec, k_vec) {
  panels <- vector("list", length(delta0_vec))
  for (idx in seq_along(delta0_vec)) {
    keys   <- sapply(k_vec, function(k) paste0("k", k, "_d0_", delta0_vec[idx]))
    ev_mat <- sapply(keys, function(key) exp(log_ev_grid[[key]]))
    pooled <- sapply(keys, function(key) exp(sum(log_ev_grid[[key]])))
    panels[[idx]] <- list(ev = ev_mat, pooled = pooled)
  }
  panels
}

panels_A <- build_panels(log_ev_grid_A, delta0_vec_A, k_vec)
panels_B <- build_panels(log_ev_grid_B, delta0_vec_B, k_vec)

## ----------------------------- scientific notation formatter ------------------
fmt_sci <- function(x, digits = 2, sci_threshold = 1e4) {
  vapply(x, function(v) {
    if (is.na(v)) return("NA")
    if (abs(v) < sci_threshold) return(formatC(v, format = "f", digits = 3))
    e <- floor(log10(abs(v)))
    m <- v / 10^e
    sprintf("$%.*f\\times 10^{%d}$", digits, m, e)
  }, character(1))
}

## ----------------------------- generic LaTeX table builder --------------------
build_latex_table <- function(panels, panel_labels, k_vec, shooter_ids,
                              caption, label, notes = NULL,
                              pooled_sci = TRUE) {
  n_k <- length(k_vec)
  lines <- character(0)
  lines <- c(lines, "\\begin{table}[H]")
  lines <- c(lines, "\\centering")
  lines <- c(lines, sprintf("\\begin{tabular}{l*{%d}{c}*{%d}{c}*{%d}{c}}", n_k, n_k, n_k))
  lines <- c(lines, "\\toprule")
  
  header1 <- paste0("Shooter & ",
                    paste(sprintf("\\multicolumn{%d}{c}{%s}", n_k, panel_labels), collapse = " & "),
                    " \\\\")
  lines <- c(lines, header1)
  
  cmid <- character(0)
  start <- 2
  for (p in seq_along(panel_labels)) {
    end <- start + n_k - 1
    cmid <- c(cmid, sprintf("\\cmidrule(lr){%d-%d}", start, end))
    start <- end + 1
  }
  lines <- c(lines, paste(cmid, collapse = ""))
  
  kheader <- paste0(" & ", paste(rep(sprintf("$k=%d$", k_vec), times = length(panel_labels)), collapse = " & "), " \\\\")
  lines <- c(lines, kheader)
  lines <- c(lines, "\\midrule")
  
  J <- length(shooter_ids)
  for (j in seq_len(J)) {
    vals <- character(0)
    for (idx in seq_along(panel_labels)) vals <- c(vals, fmt_num(panels[[idx]]$ev[j, ]))
    lines <- c(lines, paste0(shooter_ids[j], " & ", paste(vals, collapse = " & "), " \\\\"))
  }
  
  lines <- c(lines, "\\midrule")
  pooled_vals <- character(0)
  for (idx in seq_along(panel_labels)) {
    p <- panels[[idx]]$pooled
    pooled_vals <- c(pooled_vals, if (pooled_sci) fmt_sci(p) else fmt_num(p))
  }
  lines <- c(lines, paste0("Pooled & ", paste(pooled_vals, collapse = " & "), " \\\\"))
  
  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular}\\\\")
  if (!is.null(notes)) {
    lines <- c(lines, sprintf("\\footnotesize \\textit{Notes.} %s", notes))
  }
  lines <- c(lines, sprintf("\\caption{%s}", caption))
  lines <- c(lines, sprintf("\\label{%s}", label))
  lines <- c(lines, "\\end{table}")
  
  paste(lines, collapse = "\n")
}

## ----------------------------- Table A: hot-hand mixture -----------------------
panel_labels_A <- sprintf("$\\delta_0=%.1f$, $\\delta_1\\sim TN(%.1f,%.1f)$",
                          delta0_vec_A, prior_mu_vec, prior_sd)

tex_A <- build_latex_table(
  panels_A, panel_labels_A, k_vec, shooter_ids,
  caption = "Mixture e-values $\varepsilon_{\delta_0,\pi}(x)$ against $H_0: \delta \leq \delta_0$ for each shooter, by streak length $k$.",
  label   = "tab:gvt_evalues_hothand"
)

clipr::write_clip(tex_A)

## ----------------------------- Table B: sceptic's point test -------------------
panel_labels_B <- sprintf("$\\delta_0=%.1f$, $\\delta_1=%.1f$ (sceptic)",
                          delta0_vec_B, delta1_B)

notes_B <- paste0(
  "The sceptic's point e-value $\\varepsilon_{\\delta_0,0}(x)$ of Equation \\eqref{eq: point e-value} is computed exactly for each shooter via the recursion of Appendix \\ref{app: recursion}. ",
  "The fixed alternative is $\\delta_1=0$ throughout, so by Theorem \\ref{thm: one-sided e-value} each entry is an e-value against $H_0:\\delta\\ge\\delta_0$. ",
  "Pooled e-values are the product of the individual shooters' e-values."
)

tex_B <- build_latex_table(
  panels_B, panel_labels_B, k_vec, shooter_ids,
  caption = "Sceptic's point e-values $\epsilon_{\delta_0,0}(x)$ against $H_0:\delta\ge\delta_0$ for each shooter, by streak length $k$.",
  label   = "tab:gvt_evalues_sceptic",
  notes   = notes_B,
  pooled_sci = TRUE
)

clipr::write_clip(tex_B)

tex_out <- paste(tex_A, "\n\n", tex_B, collapse = "\n")
writeLines(tex_out, out_tex)
cat(sprintf("\nLaTeX tables written to %s\n", out_tex))
