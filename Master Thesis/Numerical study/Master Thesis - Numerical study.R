
################################################################################
## Simulation subsection: exact evaluation of the conditional hot-hand e-value.
##
##   eps_{d0,d1}(x) = ( Z_{d0}(h,nH) / Z_{d1}(h,nH) ) * exp( (d1 - d0) * hH(x) )
##
##   point null   D = {d0}           any d1 != d0
##   lower null   D = (-Inf, d0]     valid iff d1 > d0   (test FOR the hot hand)
##   upper null   D = [d0, +Inf)     valid iff d1 < d0   (test AGAINST it)
##
## Everything here is EXACT. The count table N(h, nH, hH) is built by the
## recursion of the appendix, so Z_delta, Lambda', the cell probabilities and
## hence E[eps] and E[log eps] are finite sums with no Monte Carlo error. The
## only randomness in the whole script is in Figure 5, where the DATA GENERATING
## law is outside the model and has to be simulated; the e-value itself is still
## read off the exact table there.
##
## Blocks:
##   0. recursion, cell table, exact functionals
##   1. validation: recursion against brute-force enumeration at small n
##   F1 validity and exactness in both directions, and the role of eta
##   F2 e-power against the alternative magnitude d1
##   F3 point bet against mixture bets, swept over the centre of the bet
################################################################################

RUN_VALIDATION <- TRUE

n_main    <- 100
k_values  <- c(1, 2, 3)
theta_main <- 0.50

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Install ggplot2 first:  install.packages('ggplot2')")
library(ggplot2)

k_colours <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73")
k_labels  <- c("1" = "k = 1", "2" = "k = 2", "3" = "k = 3")
th_colours <- c("0.35" = "#D55E00", "0.5" = "#0072B2", "0.65" = "#009E73")
th_labels  <- c("0.35" = "theta_N = 0.35", "0.5" = "theta_N = 0.50",
                "0.65" = "theta_N = 0.65")

## ============================ 0. PRIMITIVES ===================================
logit_inv <- function(z) 1 / (1 + exp(-z))

## trailing run of hits in the prefix, capped at k: the only thing about the
## prefix that the hotness of any later trial depends on
trail_run <- function(prefix, k) {
  r <- 0L
  for (j in k:1) if (prefix[j] == 1L) r <- r + 1L else break
  r
}

## Recursion of the appendix. N_t(r, h, nH, hH) counts sequences of the first t
## trials with trailing hit run min(r,k) and post-prefix counts (h, nH, hH).
## Before trial t all three counts are at most t-1-k, which is why every slice
## below is restricted to 1:u with u = t-k: that restriction is what keeps the
## cost at O(n^4) element operations rather than O(n^4) with a full array each
## step.
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

## The cell table. One row per non-empty cell (h, nH), holding the log counts of
## each value of hH. Everything downstream is vectorised over these rows, which
## is what makes a grid of a few hundred deltas cheap.
build_cells <- function(n, k, prefix = rep(1L, k)) {
  m <- n - k
  N <- count_table(n, k, prefix)
  mass <- apply(N, c(1, 2), sum)
  sel  <- which(mass > 0, arr.ind = TRUE)
  W <- matrix(0, nrow(sel), m + 1L)
  for (i in seq_len(nrow(sel))) W[i, ] <- N[sel[i, 1], sel[i, 2], ]
  logW <- log(W)                       # -Inf where a value of hH is impossible
  list(n = n, k = k, m = m,
       h = as.numeric(sel[, 1] - 1L), nH = as.numeric(sel[, 2] - 1L),
       v = 0:m, logW = logW, ncell = nrow(sel),
       idx = { M <- matrix(NA_integer_, m + 1L, m + 1L)
       M[sel] <- seq_len(nrow(sel)); M })
}

## log Z_delta and Lambda'(delta) for every cell at once, cached on delta.
lam_env <- new.env(hash = TRUE, parent = emptyenv())
lam <- function(CL, d) {
  key <- sprintf("%d|%d|%.10f", CL$n, CL$k, d)
  if (exists(key, envir = lam_env, inherits = FALSE)) return(get(key, envir = lam_env))
  A  <- sweep(CL$logW, 2, d * CL$v, "+")
  mx <- apply(A, 1, max)
  E  <- exp(A - mx)
  s  <- rowSums(E)
  out <- list(L = mx + log(s), Lp = as.numeric(E %*% CL$v) / s)
  assign(key, out, envir = lam_env)
  out
}
## conditional pmf of hH within each cell under delta: rows sum to one
cond_hH <- function(CL, d) {
  A  <- sweep(CL$logW, 2, d * CL$v, "+")
  mx <- apply(A, 1, max)
  E  <- exp(A - mx)
  E / rowSums(E)
}

## log P_{eta,delta}( X in cell ), for every cell
cell_logmass <- function(CL, eta, delta) {
  eta * CL$h -
    CL$nH * log((1 + exp(eta + delta)) / (1 + exp(eta))) -
    CL$m * log(1 + exp(eta)) +
    lam(CL, delta)$L
}

## E_{eta,delta}[ eps_{d0,d1} ]. Equals 1 at delta = d0 for every eta.
mean_eps <- function(CL, eta, delta, d0, d1) {
  a  <- d1 - d0
  lp <- cell_logmass(CL, eta, delta)
  ce <- lam(CL, d0)$L - lam(CL, d1)$L + lam(CL, delta + a)$L - lam(CL, delta)$L
  sum(exp(lp + ce))
}
## E_{eta,delta}[ log eps_{d0,d1} ], the e-power
epower <- function(CL, eta, delta, d0, d1) {
  a  <- d1 - d0
  lp <- cell_logmass(CL, eta, delta)
  sum(exp(lp) * (lam(CL, d0)$L - lam(CL, d1)$L + a * lam(CL, delta)$Lp))
}
## the same two quantities for a mixture over d1 with weights w on the grid dj
mean_eps_mix <- function(CL, eta, delta, d0, dj, w)
  sum(w * vapply(dj, function(d1) mean_eps(CL, eta, delta, d0, d1), numeric(1)))

epower_mix <- function(CL, eta, delta, d0, dj, w) {
  lp <- cell_logmass(CL, eta, delta)
  P  <- cond_hH(CL, delta)
  L0 <- lam(CL, d0)$L
  ## log eps at every (cell, hH), assembled by a two pass log-sum-exp over dj
  mx <- matrix(-Inf, CL$ncell, CL$m + 1L)
  for (j in seq_along(dj)) {
    Mj <- (L0 - lam(CL, dj[j])$L) + outer(rep(1, CL$ncell), (dj[j] - d0) * CL$v)
    mx <- pmax(mx, Mj)
  }
  S <- matrix(0, CL$ncell, CL$m + 1L)
  for (j in seq_along(dj)) {
    Mj <- (L0 - lam(CL, dj[j])$L) + outer(rep(1, CL$ncell), (dj[j] - d0) * CL$v)
    S  <- S + w[j] * exp(Mj - mx)
  }
  sum(exp(lp) * rowSums(P * (mx + log(S))))
}

## truncated normal weights on a grid, restricted to the admissible side of d0
mix_weights <- function(mu, sd, d0, side = c("lower", "upper"),
                        grid = seq(-3, 3, by = 0.02)) {
  side <- match.arg(side)
  g <- if (side == "lower") grid[grid > d0] else grid[grid < d0]
  w <- dnorm(g, mu, sd)
  list(d = g, w = w / sum(w))
}

CELLS <- list()
for (k in k_values) {
  cat(sprintf("building exact cell table: n = %d, k = %d ... ", n_main, k))
  t0 <- proc.time()[3]
  CELLS[[as.character(k)]] <- build_cells(n_main, k)
  cat(sprintf("%d cells, %.1fs\n", CELLS[[as.character(k)]]$ncell, proc.time()[3] - t0))
}
eta_of <- function(th) log(th / (1 - th))

## ======================= 1. VALIDATION OF THE RECURSION =======================
## Brute-force enumeration is used here and nowhere else: at n = 14 the whole
## sample space fits in memory, so the recursion can be checked exactly rather
## than plotted.
if (RUN_VALIDATION) {
  cat("\n---- validation ----\n")
  enum_table <- function(n, k, prefix) {
    m <- n - k
    G <- as.matrix(expand.grid(rep(list(c(0L, 1L)), m)))
    X <- cbind(matrix(as.integer(prefix), nrow(G), k, byrow = TRUE), G)
    A <- X[, 1:m, drop = FALSE]
    if (k >= 2) for (j in 2:k) A <- A * X[, j:(m + j - 1), drop = FALSE]
    Hot  <- A == 1
    post <- X[, (k + 1):n, drop = FALSE]
    N <- array(0, dim = c(m + 1L, m + 1L, m + 1L))
    h  <- rowSums(post); nH <- rowSums(Hot); hH <- rowSums(Hot * post)
    for (i in seq_len(nrow(X))) N[h[i] + 1, nH[i] + 1, hH[i] + 1] <-
      N[h[i] + 1, nH[i] + 1, hH[i] + 1] + 1
    N
  }
  for (cfg in list(list(14, 1, c(1L)), list(14, 2, c(1L, 1L)),
                   list(14, 3, c(1L, 1L, 1L)), list(13, 2, c(0L, 1L)))) {
    n <- cfg[[1]]; k <- cfg[[2]]; pf <- cfg[[3]]
    d <- max(abs(count_table(n, k, pf) - enum_table(n, k, pf)))
    cat(sprintf("  n=%d k=%d prefix=(%s): max |recursion - enumeration| = %g\n",
                n, k, paste(pf, collapse = ""), d))
  }
  cat("\n  exactness at the boundary, n = 100 (target: mass 1, E[eps] 1)\n")
  for (k in k_values) for (th in c(0.35, 0.50, 0.65)) {
    CL <- CELLS[[as.character(k)]]; e <- eta_of(th)
    cat(sprintf("    k=%d theta_N=%.2f  mass=%.12f  E[eps | delta=d0=0]=%.12f\n",
                k, th, sum(exp(cell_logmass(CL, e, 0))),
                mean_eps(CL, e, 0, 0, 0.6)))
  }
}

## ================== F1: VALIDITY IN BOTH DIRECTIONS ===========================
## Rows: the two one-sided nulls. Columns: k. Colour: theta_N. The three colours
## must meet at delta = d0 and separate everywhere else, which is the picture of
## "exactness is free of the baseline skill, power is not".
dir_spec <- list(
  lower = list(d0 = 0.0, d1 = 0.6,
               lab = "atop(H[0]*': '*delta<=0, delta[1]==0.6)"),
  upper = list(d0 = 0.6, d1 = 0.0,
               lab = "atop(H[0]*': '*delta>=0.6, delta[1]==0)")
)
delta_grid_F1 <- seq(-0.4, 1.2, by = 0.02)
theta_grid_F1 <- c(0.35, 0.50, 0.65)

F1 <- do.call(rbind, lapply(names(dir_spec), function(sd_) {
  sp <- dir_spec[[sd_]]
  do.call(rbind, lapply(k_values, function(k) {
    CL <- CELLS[[as.character(k)]]
    do.call(rbind, lapply(theta_grid_F1, function(th) {
      e <- eta_of(th)
      data.frame(delta = delta_grid_F1,
                 mean_e = vapply(delta_grid_F1, function(dt)
                   mean_eps(CL, e, dt, sp$d0, sp$d1), numeric(1)),
                 k = factor(k), theta = factor(th), side = sd_,
                 d0 = sp$d0, lab = sp$lab)
    }))
  }))
}))

rect_F1 <- data.frame(
  side = c("lower", "upper"),
  lab  = vapply(dir_spec, function(z) z$lab, character(1)),
  xmin = c(min(delta_grid_F1), 0.6),
  xmax = c(0.0, max(delta_grid_F1)),
  d0   = c(0.0, 0.6))
rect_F1 <- merge(rect_F1, expand.grid(k = factor(k_values)), by = NULL)

pF1 <- ggplot(F1, aes(delta, mean_e, colour = theta)) +
  geom_rect(data = rect_F1, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
            fill = "grey88") +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey30") +
  geom_vline(aes(xintercept = d0), linetype = "dashed", colour = "grey30") +
  geom_line(linewidth = 0.7) +
  facet_grid(lab ~ k, labeller = labeller(lab = label_parsed,
                                          k = as_labeller(k_labels))) +
  scale_colour_manual(values = th_colours, labels = th_labels, name = NULL) +
  scale_y_log10() +
  labs(x = expression(delta~"(true)"),
       y = expression(E[eta*","*delta]*"["*epsilon*"]")) +
  #subtitle = sprintf("exact, n = %d, all-hits prefix; shaded region is the null, the three curves meet only at the boundary", n_main)) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
print(pF1)

## ================== F2: THE COST OF THE ALTERNATIVE ===========================
## The peak sits at d1 = delta_true, because on each cell the e-power is
## Lam(d0) - Lam(d1) + (d1 - d0) Lam'(delta_true), whose derivative in d1 is
## Lam'(delta_true) - Lam'(d1), and Lam'' > 0 makes the root unique.
##
## The d1 grids below are deliberately extended PAST d0 into the inadmissible
## side (d1 <= d0 for the lower null, d1 >= d0 for the upper null): a one-sided
## test that "bets" on the wrong side of d0 is not a valid test, and shading
## that zone (as in F1) makes the boundary visible instead of just plotting
## only the admissible side.
## Both grids now share the SAME absolute range, [-1.00, 1.50], covering the
## union of what each row needs. d1 is the same kind of quantity (a shift on
## the delta scale) in both the lower and upper null tests, so putting both
## rows on one shared axis makes slopes/steepness directly comparable between
## rows -- not just within a row. The cost: each row's peak gets a bit less
## zoom than it had on its own tighter range.
F2_x_range <- c(-1.00, 1.50)
F2_spec <- list(
  list(side = "lower", d0 = 0.0, truths = c(0.8, 1.0, 1.2),
       d1 = seq(F2_x_range[1], F2_x_range[2], by = 0.025)),
  list(side = "upper", d0 = 1.0, truths = c(0.0, 0.2, 0.4),
       d1 = seq(F2_x_range[1], F2_x_range[2], by = 0.025))
)
F2 <- do.call(rbind, lapply(F2_spec, function(sp)
  do.call(rbind, lapply(k_values, function(k) {
    CL <- CELLS[[as.character(k)]]; e <- eta_of(theta_main)
    do.call(rbind, lapply(sp$truths, function(dt)
      data.frame(d1 = sp$d1,
                 epow = vapply(sp$d1, function(d1)
                   epower(CL, e, dt, sp$d0, d1), numeric(1)),
                 delta_true = dt, k = factor(k), d0 = sp$d0,
                 panel = sprintf("%s null, d0 = %.1f, true delta = %.1f",
                                 sp$side, sp$d0, dt))))
  }))))

cat("\n---- F2: where the e-power peaks (should equal the true delta) ----\n")
for (p in unique(F2$panel)) for (kk in levels(F2$k)) {
  s <- subset(F2, panel == p & k == kk)
  cat(sprintf("  %-44s k=%s  peak at d1 = %+.3f (true %+.2f), value %.4f\n",
              p, kk, s$d1[which.max(s$epow)], s$delta_true[1], max(s$epow)))
}

## Lower-null shading now extends to the full grid boundary (-1.00), since the
## grid itself was widened to match the upper row; the null side is still
## everything with d1 <= d0, just plotted over a longer stretch now.
rect_F2 <- data.frame(
  panel = c("lower null, d0 = 0.0, true delta = 0.8",
            "lower null, d0 = 0.0, true delta = 1.0",
            "lower null, d0 = 0.0, true delta = 1.2",
            "upper null, d0 = 1.0, true delta = 0.0",
            "upper null, d0 = 1.0, true delta = 0.2",
            "upper null, d0 = 1.0, true delta = 0.4"),
  xmin = c(F2_x_range[1], F2_x_range[1], F2_x_range[1], 1.00, 1.00, 1.00),
  xmax = c(0.00, 0.00, 0.00, F2_x_range[2], F2_x_range[2], F2_x_range[2]),
  d0   = c(0.0, 0.0, 0.0, 1.0, 1.0, 1.0))

## No manual scale_y_continuous(limits=...) here: those numbers were tuned for
## the previous truths/d0 and would silently clip the new panels (truths now
## run up to 1.2 against d0 = 1.0, a different e-power range). facet_wrap with
## no "scales" argument gives every panel the same, data-driven x- AND y-axis
## now -- fully comparable across all 6 panels, top row and bottom row alike.
pF2 <- ggplot(F2, aes(d1, epow, colour = k)) +
  geom_rect(data = rect_F2, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey88") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_vline(aes(xintercept = d0), linetype = "dashed", colour = "grey30") +
  geom_vline(aes(xintercept = delta_true), linetype = "dotted", colour = "grey55") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ panel, nrow = 2) +
  scale_colour_manual(values = k_colours, labels = k_labels, name = NULL) +
  labs(x = expression(delta[1]~"(alternative)"),
       y = expression(E*"["*log~epsilon*"]"~"(nats)"))+
       #subtitle = "exact; shaded region is inadmissible (wrong side of d0), peak sits at the true delta") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
print(pF2)

## ============ F3: POINT BET AGAINST MIXTURE BETS, SWEPT OVER THE CENTRE =======
## Sweeping the truth confounds two things. Sweeping the centre mu of the bet at
## a fixed truth is the comparison that isolates what the spread buys.
##
## mu is extended down to -0.50 (past d0 = 0) so the null side (mu <= d0, i.e.
## centring the bet on or below the null) is visible and can be shaded, exactly
## as in F1. The mixture weights themselves stay truncated to the admissible
## side (d1 > d0) via mix_weights(), so a mixture centred at a negative mu just
## becomes a weak, near-degenerate bet piled up near d1 = 0 -- which is itself
## informative to see.
F3_mu     <- seq(0.05, 1.50, by = 0.05)
F3_sd     <- c(0, 0.3, 0.5)
F3_truths <- c(0.8, 1.2)
F3_k      <- 1
CL3 <- CELLS[[as.character(F3_k)]]; eta3 <- eta_of(theta_main)

F3 <- do.call(rbind, lapply(F3_truths, function(dt)
  do.call(rbind, lapply(F3_sd, function(s)
    data.frame(mu = F3_mu, sd = factor(s), delta_true = dt,
               epow = vapply(F3_mu, function(mu) {
                 if (s == 0) return(epower(CL3, eta3, dt, 0, mu))
                 mw <- mix_weights(mu, s, 0, "lower")
                 epower_mix(CL3, eta3, dt, 0, mw$d, mw$w)
               }, numeric(1)))))))

cat("\n---- F3: e-power at the correct centre, and the hedge premium ----\n")
for (dt in F3_truths) {
  s0 <- subset(F3, delta_true == dt & sd == "0"   & abs(mu - dt) < 1e-9)$epow
  for (s in F3_sd[-1]) {
    ss <- subset(F3, delta_true == dt & sd == as.character(s) & abs(mu - dt) < 1e-9)$epow
    cat(sprintf("  true delta = %.1f: point %.4f, sd = %.1f %.4f (%.1f%% of the point bet)\n",
                dt, s0, s, ss, 100 * ss / s0))
  }
}

## Null region shaded as mu <= d0 = 0 (betting centred on or past the null),
## and the two facets now share a fixed y-scale by default (facet_wrap with no
## "scales" argument), so e-power is directly comparable between the two true-
## delta panels, exactly as requested for F1.
rect_F3 <- data.frame(delta_true = F3_truths, xmin = -0.2, xmax = 0)

pF3 <- ggplot(F3, aes(mu, epow, colour = sd)) +
  geom_rect(data = rect_F3, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey88") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  geom_vline(aes(xintercept = delta_true), linetype = "dotted", colour = "grey55") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ delta_true,
             labeller = as_labeller(setNames(
               paste0("true delta = ", F3_truths), F3_truths))) +
  scale_colour_manual(values = c("0" = "#D55E00", "0.3" = "#0072B2", "0.5" = "#009E73"),
                      labels = c("0" = "point", "0.3" = "sd = 0.3", "0.5" = "sd = 0.5"),
                      name = "prior over the alternative") +
  labs(x = expression(mu~"(centre of the bet)"),
       y = expression(E*"["*log~epsilon*"]"~"(nats)"))+
       #subtitle = "exact, n = 100, k = 1, d0 = 0; shaded region is the null side of the bet (mu <= d0)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
print(pF3)