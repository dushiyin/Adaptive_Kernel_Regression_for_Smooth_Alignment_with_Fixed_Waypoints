set.seed(2025)
n      <- 1200
m_ctrl <- 3
ngrid  <- 500

h_grid      <- seq(0.04, 0.12, length.out = 15)
lambda_grid <- seq(0, 100, 5)

alpha_loss <- 1
Kfold <- 5

tau_c <- 0.10
tau_g <- 0.05
w_c   <- 0.50
w_g   <- 0.30
w_s   <- 0.20

if (!dir.exists("data")) dir.create("data")
if (!dir.exists("figs")) dir.create("figs")

m_true <- function(s) 0.5 * sin(8*s)

s <- sort(runif(n, 0, 1))

x <- 100*s + 5*sin(6*s)
y <- 100*(s^1.2) + 3*cos(4*s)

z <- m_true(s) + rnorm(n, mean = 0, sd = 0.08)

sg <- seq(0, 1, length.out = ngrid)
truth_on_grid <- m_true(sg)
pick_constraints_2d <- function(s, z, m) {
  idx <- sort(sample(seq_along(s), m))
  sC  <- s[idx]
  zC  <- z[idx]
  data.frame(sC = sC, zC = zC, idx = idx)
}
C2 <- pick_constraints_2d(s, z, m_ctrl)

K <- function(u) dnorm(u)

nw_2d <- function(s0, s, z, h) {
  w <- K((s0 - s) / h)
  sw <- sum(w)
  if (sw == 0) return(0)
  sum(w * z) / sw
}

anw_2d <- function(s0, s, z, sC, zC, h, lambda) {
  w  <- K((s0 - s) / h)
  wC <- lambda * K((s0 - sC) / h)
  num <- sum(w * z) + sum(wC * zC)
  den <- sum(w) + sum(wC)
  if (den == 0) return(0)
  num / den
}

pred_on_grid <- function(sg, fun) vapply(sg, fun, numeric(1))

rmse <- function(truth, pred) sqrt(mean((truth - pred)^2))

constraint_rmse <- function(sC, zC, pred_fun_point) {
  zhatC <- vapply(sC, pred_fun_point, numeric(1))
  sqrt(mean((zC - zhatC)^2))
}

smooth_index <- function(z_hat_grid) {
  d2 <- diff(z_hat_grid, differences = 2)
  mean(d2^2)
}

phi <- function(z) pmax(0, z - 1)

z_mean <- mean(z)
z_sd   <- sd(z)
z_std  <- (z - z_mean) / z_sd
zC_std <- (C2$zC - z_mean) / z_sd

set.seed(2025)
fold_id <- sample(rep(1:Kfold, length.out = n))

anw_pred_vec_2d <- function(s_train, z_train, s_eval, sC, zC, h, lambda) {
  vapply(
    s_eval,
    function(s0) anw_2d(s0, s_train, z_train, sC, zC, h, lambda),
    numeric(1)
  )
}

ds_anw_pred_vec_2d <- function(s_train, z0_std, s_eval, sC, zC_std, h, lambda, M = 1) {
  if (M == 0) {
    return(anw_pred_vec_2d(s_train, z0_std, s_eval, sC, zC_std, h, lambda))
  }
  zm <- z0_std
  for (m in 1:M) {
    fit_on_train <- anw_pred_vec_2d(s_train, zm, s_train, sC, zC_std, h, lambda)
    zm <- z0_std + (zm - fit_on_train)  
  }
  anw_pred_vec_2d(s_train, zm, s_eval, sC, zC_std, h, lambda)
}

cv_nw_std <- function(h) {
  mse_vec <- numeric(Kfold)
  for (k in 1:Kfold) {
    tr <- which(fold_id != k)
    te <- which(fold_id == k)
    
    pred_te <- vapply(
      s[te],
      function(s0) nw_2d(s0, s[tr], z_std[tr], h),
      numeric(1)
    )
    mse_vec[k] <- mean((z_std[te] - pred_te)^2)
  }
  mean(mse_vec)
}

loss_anw_std <- function(h, lambda, alpha = 1) {
  if (h <= 0 || lambda < 0) return(list(cv = Inf, constr = Inf, loss = Inf))
  
  mse_vec <- numeric(Kfold)
  for (k in 1:Kfold) {
    tr <- which(fold_id != k)
    te <- which(fold_id == k)
    
    pred_te <- anw_pred_vec_2d(s[tr], z_std[tr], s[te], C2$sC, zC_std, h, lambda)
    mse_vec[k] <- mean((z_std[te] - pred_te)^2)
  }
  cv_part <- mean(mse_vec)
  
  predC_std <- anw_pred_vec_2d(s, z_std, C2$sC, C2$sC, zC_std, h, lambda)
  constr_part <- mean((predC_std - zC_std)^2)
  
  list(cv = cv_part, constr = constr_part, loss = cv_part + alpha * constr_part)
}

loss_ds_anw_std <- function(h, lambda, M, alpha = 1) {
  if (h <= 0 || lambda < 0) return(list(cv = Inf, constr = Inf, loss = Inf))
  
  mse_vec <- numeric(Kfold)
  for (k in 1:Kfold) {
    tr <- which(fold_id != k)
    te <- which(fold_id == k)
    
    s_tr <- s[tr]
    z_tr <- z_std[tr]
    
    C_tr_idx <- which(tr %in% C2$idx)
    sC_tr    <- s_tr[C_tr_idx]
    zC_tr    <- z_tr[C_tr_idx]
    
    if (length(C_tr_idx) == 0) {
      sC_tr <- numeric(0)
      zC_tr <- numeric(0)
    }
    
    pred_te <- ds_anw_pred_vec_2d(
      s_train = s_tr, z0_std = z_tr, s_eval = s[te],
      sC = sC_tr, zC_std = zC_tr,
      h = h, lambda = lambda, M = M
    )
    
    mse_vec[k] <- mean((z_std[te] - pred_te)^2)
  }
  cv_part <- mean(mse_vec)
  
  predC_std <- ds_anw_pred_vec_2d(
    s_train = s, z0_std = z_std, s_eval = C2$sC,
    sC = C2$sC, zC_std = zC_std,
    h = h, lambda = lambda, M = M
  )
  constr_part <- mean((predC_std - zC_std)^2)
  
  list(cv = cv_part, constr = constr_part, loss = cv_part + alpha * constr_part)
}

cv_vals_nw <- sapply(h_grid, cv_nw_std)
h_nw_opt   <- h_grid[which.min(cv_vals_nw)]
cat("Optimal h (NW) =", h_nw_opt, "\n")
best_loss_anw <- Inf
best_anw <- list(h = NA_real_, L = NA_real_, cv = NA_real_, constr = NA_real_)

loss_records <- list()
ctr <- 0

for (h0 in h_grid) {
  for (L0 in lambda_grid) {
    out <- loss_anw_std(h0, L0, alpha = alpha_loss)
    ctr <- ctr + 1
    loss_records[[ctr]] <- data.frame(
      Method = "ANW", h = h0, lambda = L0,
      cv = out$cv, constr = out$constr, loss = out$loss
    )
    if (out$loss < best_loss_anw) {
      best_loss_anw <- out$loss
      best_anw <- list(h = h0, L = L0, cv = out$cv, constr = out$constr)
    }
  }
}

cat("Optimal (h, lambda) for ANW:\n")
cat("  h_opt      =", best_anw$h, "\n")
cat("  lambda_opt =", best_anw$L, "\n")
cat("  Loss       =", best_loss_anw, "\n")

best_loss_ds1 <- Inf
best_ds1 <- list(h = NA_real_, L = NA_real_, cv = NA_real_, constr = NA_real_)

for (h0 in h_grid) {
  for (L0 in lambda_grid) {
    out <- loss_ds_anw_std(h0, L0, M = 1, alpha = alpha_loss)
    ctr <- ctr + 1
    loss_records[[ctr]] <- data.frame(
      Method = "DS-ANW (M=1)", h = h0, lambda = L0,
      cv = out$cv, constr = out$constr, loss = out$loss
    )
    if (out$loss < best_loss_ds1) {
      best_loss_ds1 <- out$loss
      best_ds1 <- list(h = h0, L = L0, cv = out$cv, constr = out$constr)
    }
  }
}

cat("Optimal (h, lambda) for DS-ANW (M=1):\n")
cat("  h_opt      =", best_ds1$h, "\n")
cat("  lambda_opt =", best_ds1$L, "\n")
cat("  Loss       =", best_loss_ds1, "\n")

best_loss_ds2 <- Inf
best_ds2 <- list(h = NA_real_, L = NA_real_, cv = NA_real_, constr = NA_real_)

for (h0 in h_grid) {
  for (L0 in lambda_grid) {
    out <- loss_ds_anw_std(h0, L0, M = 2, alpha = alpha_loss)
    ctr <- ctr + 1
    loss_records[[ctr]] <- data.frame(
      Method = "DS-ANW (M=2)", h = h0, lambda = L0,
      cv = out$cv, constr = out$constr, loss = out$loss
    )
    if (out$loss < best_loss_ds2) {
      best_loss_ds2 <- out$loss
      best_ds2 <- list(h = h0, L = L0, cv = out$cv, constr = out$constr)
    }
  }
}

cat("Optimal (h, lambda) for DS-ANW (M=2):\n")
cat("  h_opt      =", best_ds2$h, "\n")
cat("  lambda_opt =", best_ds2$L, "\n")
cat("  Loss       =", best_loss_ds2, "\n")

loss_df <- do.call(rbind, loss_records)
write.csv(loss_df, "data/exp2d_loss_surface.csv", row.names = FALSE)

pred_nw <- pred_on_grid(
  sg,
  function(s0) nw_2d(s0, s, z, h_nw_opt)
)

pred_anw <- pred_on_grid(
  sg,
  function(s0) anw_2d(s0, s, z, C2$sC, C2$zC, best_anw$h, best_anw$L)
)

pred_ds1_std <- ds_anw_pred_vec_2d(
  s_train = s, z0_std = z_std, s_eval = sg,
  sC = C2$sC, zC_std = zC_std,
  h = best_ds1$h, lambda = best_ds1$L, M = 1
)
pred_ds1 <- z_mean + z_sd * pred_ds1_std

pred_ds2_std <- ds_anw_pred_vec_2d(
  s_train = s, z0_std = z_std, s_eval = sg,
  sC = C2$sC, zC_std = zC_std,
  h = best_ds2$h, lambda = best_ds2$L, M = 2
)
pred_ds2 <- z_mean + z_sd * pred_ds2_std

tab2 <- data.frame(
  Model = c("NW", "ANW", "DS-ANW (M=1)", "DS-ANW (M=2)"),
  RMSE  = c(rmse(truth_on_grid, pred_nw),
            rmse(truth_on_grid, pred_anw),
            rmse(truth_on_grid, pred_ds1),
            rmse(truth_on_grid, pred_ds2)),
  CErr  = c(
    constraint_rmse(C2$sC, C2$zC, function(s0) nw_2d(s0, s, z, h_nw_opt)),
    constraint_rmse(C2$sC, C2$zC, function(s0) anw_2d(s0, s, z, C2$sC, C2$zC, best_anw$h, best_anw$L)),
    constraint_rmse(C2$sC, C2$zC, function(s0){
      zhat_std <- ds_anw_pred_vec_2d(
        s_train = s, z0_std = z_std, s_eval = s0,
        sC = C2$sC, zC_std = zC_std,
        h = best_ds1$h, lambda = best_ds1$L, M = 1
      )
      z_mean + z_sd * zhat_std
    }),
    constraint_rmse(C2$sC, C2$zC, function(s0){
      zhat_std <- ds_anw_pred_vec_2d(
        s_train = s, z0_std = z_std, s_eval = s0,
        sC = C2$sC, zC_std = zC_std,
        h = best_ds2$h, lambda = best_ds2$L, M = 2
      )
      z_mean + z_sd * zhat_std
    })
  ),
  Smoothness = c(
    smooth_index(pred_nw),
    smooth_index(pred_anw),
    smooth_index(pred_ds1),
    smooth_index(pred_ds2)
  )
)

tau_s <- median(tab2$Smoothness)

tab2$CSS <- with(tab2,
                 w_c * phi(CErr / tau_c) +
                   w_g * phi(RMSE / tau_g) +
                   w_s * phi(Smoothness / tau_s))

write.csv(tab2, "data/exp2d_summary.csv", row.names = FALSE)

pdf("figs/fig_2d_main.pdf", width = 7, height = 5)
par(mar = c(4, 4, 2, 1))

plot(sg, truth_on_grid, type = "l", lwd = 3, col = "black",
     xlab = "s (track parameter)", ylab = "z(s)",
     main = "2D (parametric) reconstruction: NW vs ANW vs DS-ANW")

lines(sg, pred_nw,  col = "dodgerblue", lwd = 2, lty = 1)
lines(sg, pred_anw, col = "firebrick",  lwd = 2, lty = 2)
lines(sg, pred_ds1, col = "darkgreen",  lwd = 2, lty = 3)
lines(sg, pred_ds2, col = "orange",  lwd = 2, lty = 4)

points(s, z, pch = 16, cex = 0.5, col = rgb(0,0,0,0.15))
points(C2$sC, C2$zC, pch = 16, col = "red", cex = 1.4)

anw_lab <- bquote(ANW~"(h="*.(sprintf("%.3f", best_anw$h))*", "*lambda*"="*.(sprintf("%g", best_anw$L))*")")
ds1_lab <- bquote(DS-ANW~"(M=1,"~h==" "*.(sprintf("%.3f", best_ds1$h))*","~lambda==" "*.(sprintf("%g", best_ds1$L))*")")
ds2_lab <- bquote(DS-ANW~"(M=2,"~h==" "*.(sprintf("%.3f", best_ds2$h))*","~lambda==" "*.(sprintf("%g", best_ds2$L))*")")

legend("bottomleft",
       legend = c("Truth", "NW", anw_lab, ds1_lab, ds2_lab, "Control points", "Samples"),
       col    = c("black", "dodgerblue", "firebrick", "darkgreen", "orange", "red", rgb(0,0,0,0.15)),
       lty    = c(1, 1, 2, 3, 4, NA, NA),
       pch    = c(NA, NA, NA, NA, NA, 16, 16),
       lwd    = c(3, 2, 2, 2, 2, NA, NA),
       bty    = "n", cex = 0.85)

grid()
dev.off()

pdf("figs/fig_2d_track.pdf", width = 6, height = 6)
par(mar = c(4, 4, 2, 1))

plot(x, y, type = "l", lwd = 2, col = "gray40",
     xlab = "x(s)", ylab = "y(s)",
     main = "Track geometry with control points")

points(x[C2$idx], y[C2$idx], pch = 16, col = "red", cex = 1.4)

legend("bottomright",
       legend = c("track", "control points"),
       col = c("gray40", "red"),
       lty = c(1, NA), pch = c(NA, 16),
       bty = "n")
grid()
dev.off()

cat("\n=== 2D  ===\n")
print(tab2)
cat("\n print：data/exp2d_summary.csv\n")
cat("print：figs/fig_2d_main.pdf, figs/fig_2d_track.pdf\n")
cat("print：data/exp2d_loss_surface.csv (Loss surface; ANW/DS1/DS2)\n")
