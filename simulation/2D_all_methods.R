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

suppressPackageStartupMessages({
  library(splines)
  library(mgcv)
  library(GauPro)
})

m_true <- function(s) 0.5 * sin(8*s)

s <- sort(runif(n, 0, 1))

x <- 100*s + 5*sin(6*s)
y <- 100*(s^1.2) + 3*cos(4*s)

z <- m_true(s) + rnorm(n, mean = 0, sd = 0.08)

sg <- seq(0, 1, length.out = ngrid)
truth_on_grid <- m_true(sg)

pick_constraints_2d <- function(s, m) {
  idx <- sort(sample(seq_along(s), m))
  sC  <- s[idx]
  zC  <- z[idx]     
  data.frame(sC = sC, zC = zC, idx = idx)
}
C2 <- pick_constraints_2d(s, m_ctrl)

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
inv_std <- function(z_star) z_mean + z_sd * z_star

set.seed(2025)
fold_id <- sample(rep(1:Kfold, length.out = n))

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

cv_anw_loss_std <- function(h, lambda, alpha = 1) {
  if (h <= 0 || lambda < 0) {
    return(list(cv = Inf, constr = Inf, loss = Inf))
  }
  
  mse_vec <- numeric(Kfold)
  for (k in 1:Kfold) {
    tr <- which(fold_id != k)
    te <- which(fold_id == k)
    
    pred_te <- vapply(
      s[te],
      function(s0) anw_2d(s0, s[tr], z_std[tr], C2$sC, zC_std, h, lambda),
      numeric(1)
    )
    mse_vec[k] <- mean((z_std[te] - pred_te)^2)
  }
  cv_part <- mean(mse_vec)
  
  predC_std <- vapply(
    C2$sC,
    function(s0) anw_2d(s0, s, z_std, C2$sC, zC_std, h, lambda),
    numeric(1)
  )
  constr_part <- mean((predC_std - zC_std)^2)
  
  list(cv = cv_part, constr = constr_part, loss = cv_part + alpha * constr_part)
}

DS_anw_eval <- function(s_eval, s_tr, y_tr, sC, zC, h, lambda) {
  vapply(
    s_eval,
    function(s0) anw_2d(s0, s_tr, y_tr, sC, zC, h, lambda),
    numeric(1)
  )
}

DS_sharpen_train <- function(s_tr, y_base_tr, sC, zC, h, lambda, M) {
  y_m <- y_base_tr
  if (M <= 0) return(y_m)
  
  for (iter in 1:M) {
    fit_m <- DS_anw_eval(s_tr, s_tr, y_m, sC, zC, h, lambda)
    y_m <- y_base_tr + (y_m - fit_m)
  }
  y_m
}

DS_cv_part_std <- function(h, lambda, M) {
  if (h <= 0 || lambda < 0) return(Inf)
  
  mse_vec <- numeric(Kfold)
  for (k in 1:Kfold) {
    tr <- which(fold_id != k)
    te <- which(fold_id == k)
    
    y_base_tr <- z_std[tr]
    
    y_sh_tr <- DS_sharpen_train(
      s_tr = s[tr],
      y_base_tr = y_base_tr,
      sC = C2$sC,
      zC = zC_std,
      h = h, lambda = lambda, M = M
    )
    
    pred_te <- DS_anw_eval(
      s_eval = s[te],
      s_tr   = s[tr],
      y_tr   = y_sh_tr,
      sC     = C2$sC,
      zC     = zC_std,
      h      = h,
      lambda = lambda
    )
    
    mse_vec[k] <- mean((z_std[te] - pred_te)^2)
  }
  mean(mse_vec)
}

DS_constr_part_std <- function(h, lambda, M) {
  if (h <= 0 || lambda < 0) return(Inf)
  
  y_sh_full <- DS_sharpen_train(
    s_tr = s,
    y_base_tr = z_std,
    sC = C2$sC,
    zC = zC_std,
    h = h, lambda = lambda, M = M
  )
  
  predC_std <- DS_anw_eval(
    s_eval = C2$sC,
    s_tr   = s,
    y_tr   = y_sh_full,
    sC     = C2$sC,
    zC     = zC_std,
    h      = h,
    lambda = lambda
  )
  
  mean((predC_std - zC_std)^2)
}

DS_loss_std <- function(h, lambda, M, alpha = 1) {
  cv_part <- DS_cv_part_std(h, lambda, M)
  constr_part <- DS_constr_part_std(h, lambda, M)
  list(cv = cv_part, constr = constr_part, loss = cv_part + alpha * constr_part)
}

grid_search_loss <- function(filename, loss_fun) {
  loss_records <- list()
  ctr <- 0
  
  best_loss <- Inf
  h_opt  <- NA_real_
  L_opt  <- NA_real_
  best_cv <- NA_real_
  best_constr <- NA_real_
  
  for (ih in seq_along(h_grid)) {
    for (il in seq_along(lambda_grid)) {
      h0 <- h_grid[ih]
      L0 <- lambda_grid[il]
      
      out <- loss_fun(h0, L0)
      
      ctr <- ctr + 1
      loss_records[[ctr]] <- data.frame(
        h = h0,
        lambda = L0,
        cv = out$cv,
        constr = out$constr,
        loss = out$loss
      )
      
      if (out$loss < best_loss) {
        best_loss <- out$loss
        h_opt <- h0
        L_opt <- L0
        best_cv <- out$cv
        best_constr <- out$constr
      }
    }
  }
  
  loss_df <- do.call(rbind, loss_records)
  write.csv(loss_df, filename, row.names = FALSE)
  
  list(h = h_opt, lambda = L_opt,
       best_loss = best_loss, best_cv = best_cv, best_constr = best_constr,
       surface = loss_df)
}

cv_vals_nw <- sapply(h_grid, cv_nw_std)
h_nw_opt   <- h_grid[which.min(cv_vals_nw)]
cat("Optimal h (NW) =", h_nw_opt, "\n")

anw_fit <- grid_search_loss(
  filename = "data/exp2d_loss_surface.csv",
  loss_fun = function(h0, L0) cv_anw_loss_std(h0, L0, alpha = alpha_loss)
)

h_opt <- anw_fit$h
L_opt <- anw_fit$lambda
cat("Optimal (h, lambda) (ANW by Loss):\n")
cat("  h_opt      =", h_opt, "\n")
cat("  lambda_opt =", L_opt, "\n")
cat("  Loss       =", anw_fit$best_loss, "\n")
cat("  CV part    =", anw_fit$best_cv, "\n")
cat("  Constr part=", anw_fit$best_constr, "\n")

ds1_fit <- grid_search_loss(
  filename = "data/exp2d_loss_surface_ds1.csv",
  loss_fun = function(h0, L0) DS_loss_std(h0, L0, M = 1, alpha = alpha_loss)
)

h_ds1 <- ds1_fit$h
L_ds1 <- ds1_fit$lambda
cat("Optimal (h, lambda) (DS-ANW M=1 by its own Loss):\n")
cat("  h_ds1      =", h_ds1, "\n")
cat("  lambda_ds1 =", L_ds1, "\n")
cat("  Loss       =", ds1_fit$best_loss, "\n")
cat("  CV part    =", ds1_fit$best_cv, "\n")
cat("  Constr part=", ds1_fit$best_constr, "\n")

ds2_fit <- grid_search_loss(
  filename = "data/exp2d_loss_surface_ds2.csv",
  loss_fun = function(h0, L0) DS_loss_std(h0, L0, M = 2, alpha = alpha_loss)
)

h_ds2 <- ds2_fit$h
L_ds2 <- ds2_fit$lambda
cat("Optimal (h, lambda) (DS-ANW M=2 by its own Loss):\n")
cat("  h_ds2      =", h_ds2, "\n")
cat("  lambda_ds2 =", L_ds2, "\n")
cat("  Loss       =", ds2_fit$best_loss, "\n")
cat("  CV part    =", ds2_fit$best_cv, "\n")
cat("  Constr part=", ds2_fit$best_constr, "\n")

pred_nw <- pred_on_grid(
  sg,
  function(s0) nw_2d(s0, s, z, h_nw_opt)
)

pred_anw <- pred_on_grid(
  sg,
  function(s0) anw_2d(s0, s, z, C2$sC, C2$zC, h_opt, L_opt)
)

DS_pred_grid_original <- function(M, h, lambda) {
  y_sh_full <- DS_sharpen_train(
    s_tr = s,
    y_base_tr = z_std,
    sC = C2$sC,
    zC = zC_std,
    h = h, lambda = lambda, M = M
  )
  pred_std_grid <- DS_anw_eval(
    s_eval = sg,
    s_tr   = s,
    y_tr   = y_sh_full,
    sC     = C2$sC,
    zC     = zC_std,
    h      = h,
    lambda = lambda
  )
  inv_std(pred_std_grid)
}

pred_ds1 <- DS_pred_grid_original(M = 1, h = h_ds1, lambda = L_ds1)
pred_ds2 <- DS_pred_grid_original(M = 2, h = h_ds2, lambda = L_ds2)

spline_predict_vec <- function(x_train, y_train, x_out) {
  fit <- smooth.spline(x_train, y_train)
  as.numeric(predict(fit, x_out)$y)
}
spline_predict_point <- function(x_train, y_train, x0) {
  as.numeric(predict(smooth.spline(x_train, y_train), x0)$y)
}

bspline_predict_vec <- function(X_data, Y_data, x_out) {
  B <- bs(X_data)
  fit <- lm(Y_data ~ B)
  B_new <- bs(x_out)
  as.numeric(predict(fit, newdata = list(B = B_new)))
}
bspline_predict_point <- function(X_data, Y_data, x0) {
  B <- bs(X_data)
  fit <- lm(Y_data ~ B)
  as.numeric(predict(fit, newdata = list(B = bs(x0))))
}

pspline_predict_vec <- function(X_data, Y_data, x_out) {
  fit <- gam(Y_data ~ s(X_data, bs="ps", k=20, m=c(2,2)), method="REML")
  as.numeric(predict(fit, newdata = data.frame(X_data=x_out)))
}
pspline_predict_point <- function(X_data, Y_data, x0) {
  fit <- gam(Y_data ~ s(X_data, bs="ps", k=20, m=c(2,2)), method="REML")
  as.numeric(predict(fit, newdata = data.frame(X_data=x0)))
}

LOESS_predict_vec <- function(X_data, Y_data, x_out) {
  fit <- loess(Y_data ~ X_data, span=0.7, degree=2)
  as.numeric(predict(fit, data.frame(X_data=x_out)))
}
LOESS_predict_point <- function(X_data, Y_data, x0) {
  fit <- loess(Y_data ~ X_data, span=0.7, degree=2)
  as.numeric(predict(fit, data.frame(X_data=x0)))
}

GPR_predict_vec <- function(X_data, Y_data, x_out) {
  fit <- gpkm(X_data, Y_data)
  as.numeric(fit$predict(x_out))
}
GPR_predict_point <- function(X_data, Y_data, x0) {
  fit <- gpkm(X_data, Y_data)
  as.numeric(fit$predict(x0))
}

ATPS_predict_vec <- function(X_data, Y_data, x_out) {
  fit <- gam(Y_data ~ s(X_data, bs="ts", k=20) + s(X_data, bs="ts", k=20, m=2), method="REML")
  as.numeric(predict(fit, data.frame(X_data=x_out)))
}
ATPS_predict_point <- function(X_data, Y_data, x0) {
  fit <- gam(Y_data ~ s(X_data, bs="ts", k=20) + s(X_data, bs="ts", k=20, m=2), method="REML")
  as.numeric(predict(fit, data.frame(X_data=x0)))
}

pred_sp    <- spline_predict_vec(s, z, sg)
pred_bs    <- bspline_predict_vec(s, z, sg)
pred_ps    <- pspline_predict_vec(s, z, sg)
pred_LOESS <- LOESS_predict_vec(s, z, sg)
pred_gpr   <- GPR_predict_vec(s, z, sg)
pred_atps  <- ATPS_predict_vec(s, z, sg)

tab2 <- data.frame(
  Model = c("NW", "ANW", "DS-ANW (M=1)", "DS-ANW (M=2)",
            "Spline","BS","PS","LOESS","GPR","ATPS"),
  RMSE  = c(
    rmse(truth_on_grid, pred_nw),
    rmse(truth_on_grid, pred_anw),
    rmse(truth_on_grid, pred_ds1),
    rmse(truth_on_grid, pred_ds2),
    rmse(truth_on_grid, pred_sp),
    rmse(truth_on_grid, pred_bs),
    rmse(truth_on_grid, pred_ps),
    rmse(truth_on_grid[-c(1,length(truth_on_grid))], pred_LOESS[-c(1,length(truth_on_grid))]),
    rmse(truth_on_grid, pred_gpr),
    rmse(truth_on_grid, pred_atps)
  ),
  CErr  = c(
    constraint_rmse(C2$sC, C2$zC, function(s0) nw_2d(s0, s, z, h_nw_opt)),
    constraint_rmse(C2$sC, C2$zC, function(s0) anw_2d(s0, s, z, C2$sC, C2$zC, h_opt, L_opt)),
    constraint_rmse(C2$sC, C2$zC, function(s0) {
      y_sh_full <- DS_sharpen_train(s, z_std, C2$sC, zC_std, h_ds1, L_ds1, M = 1)
      inv_std(anw_2d(s0, s, y_sh_full, C2$sC, zC_std, h_ds1, L_ds1))
    }),
    constraint_rmse(C2$sC, C2$zC, function(s0) {
      y_sh_full <- DS_sharpen_train(s, z_std, C2$sC, zC_std, h_ds2, L_ds2, M = 2)
      inv_std(anw_2d(s0, s, y_sh_full, C2$sC, zC_std, h_ds2, L_ds2))
    }),
    constraint_rmse(C2$sC, C2$zC, function(xx) spline_predict_point(s, z, xx)),
    constraint_rmse(C2$sC, C2$zC, function(xx) bspline_predict_point(s, z, xx)),
    constraint_rmse(C2$sC, C2$zC, function(xx) pspline_predict_point(s, z, xx)),
    constraint_rmse(C2$sC, C2$zC, function(xx) LOESS_predict_point(s, z, xx)),
    constraint_rmse(C2$sC, C2$zC, function(xx) GPR_predict_point(s, z, xx)),
    constraint_rmse(C2$sC, C2$zC, function(xx) ATPS_predict_point(s, z, xx))
  ),
  Smoothness = c(
    smooth_index(pred_nw),
    smooth_index(pred_anw),
    smooth_index(pred_ds1),
    smooth_index(pred_ds2),
    smooth_index(pred_sp),
    smooth_index(pred_bs),
    smooth_index(pred_ps),
    smooth_index(pred_LOESS[-c(1,length(truth_on_grid))]),
    smooth_index(pred_gpr),
    smooth_index(pred_atps)
  )
)

tau_s <- median(tab2$Smoothness)

tab2$CSS <- with(tab2,
                 w_c * phi(CErr / tau_c) +
                   w_g * phi(RMSE / tau_g) +
                   w_s * phi(Smoothness / tau_s))

write.csv(tab2, "data/exp2d_summary.csv", row.names = FALSE)

pdf("figs/fig_2d_main_all.pdf", width = 7, height = 5)
par(mar = c(4, 4, 2, 1))

plot(sg, truth_on_grid, type = "l", lwd = 3, col = "black",
     xlab = "s (track parameter)", ylab = "z(s)",
     main = "2D (parametric) reconstruction: methods comparison")

lines(sg, pred_nw,   col = "dodgerblue",  lwd = 2, lty = 1)

lines(sg, pred_anw,  col = "firebrick",   lwd = 2, lty = 2)
lines(sg, pred_ds1,  col = "darkorange3", lwd = 2, lty = 3)
lines(sg, pred_ds2,  col = "purple3",     lwd = 2, lty = 4)

lines(sg, pred_sp,   col = "seagreen4", lwd = 2, lty = 5)
lines(sg, pred_bs,   col = "purple",    lwd = 2, lty = 6)
lines(sg, pred_ps,   col = "orange",    lwd = 2, lty = 7)
lines(sg, pred_LOESS,col = "blue",      lwd = 2, lty = 8)
lines(sg, pred_gpr,  col = "gold",      lwd = 2, lty = 9)
lines(sg, pred_atps, col = "skyblue",   lwd = 2, lty = 10)

points(s, z, pch = 16, cex = 0.5, col = rgb(0,0,0,0.15))
points(C2$sC, C2$zC, pch = 16, col = "red", cex = 1.4)

anw_lab <- bquote(ANW~"(h="*.(sprintf("%.3f", h_opt))*", "*lambda*"="*.(sprintf("%g", L_opt))*")")
ds1_lab <- bquote(DS-ANW~"(M=1, h="*.(sprintf("%.3f", h_ds1))*", "*lambda*"="*.(sprintf("%g", L_ds1))*")")
ds2_lab <- bquote(DS-ANW~"(M=2, h="*.(sprintf("%.3f", h_ds2))*", "*lambda*"="*.(sprintf("%g", L_ds2))*")")

legend("bottomleft",
       legend = c("Truth", "NW", anw_lab, ds1_lab, ds2_lab,
                  "Spline","BS","PS","LOESS","GPR","ATPS",
                  "Control points", "Samples"),
       col    = c("black","dodgerblue","firebrick","darkorange3","purple3",
                  "seagreen4","purple","orange","blue","gold","skyblue",
                  "red", rgb(0,0,0,0.15)),
       lty    = c(1,1,2,3,4,5,6,7,8,9,10, NA, NA),
       pch    = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA, 16, 16),
       lwd    = c(3,2,2,2,2,2,2,2,2,2,2, NA, NA),
       bty    = "n", cex = 0.7)

grid()
dev.off()

pdf("figs/fig_2d_track_all.pdf", width = 6, height = 6)
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

cat("\n=== 2D DS-ANW ===\n")
print(tab2)
cat("\n print：data/exp2d_summary.csv\n")
cat("print：figs/fig_2d_main_all.pdf, figs/fig_2d_track_all.pdf\n")
cat("print：\n")
cat("  data/exp2d_loss_surface.csv (ANW)\n")
cat("  data/exp2d_loss_surface_ds1.csv (DS-ANW M=1)\n")
cat("  data/exp2d_loss_surface_ds2.csv (DS-ANW M=2)\n")
