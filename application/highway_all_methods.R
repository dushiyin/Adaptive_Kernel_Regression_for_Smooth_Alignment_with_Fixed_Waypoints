suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(geosphere)
  library(ggplot2)
  library(ggrepel)
  library(splines)
  library(mgcv)
  library(GauPro)
})

shp_path <- "D:/Judy File/UBCO/research/Paper/data/Saskatchewan/HIGHWAY_OFFICIAL.shp"
roads_sf <- st_read(shp_path, quiet = TRUE)

hwy11_sf <- roads_sf %>%
  filter(RTNUMBER1 == 11 | RTNUMBER2 == 11 | RTNUMBER3 == 11 | RTNUMBER4 == 11)

hwy11_ll <- st_transform(hwy11_sf, 4326)  

coords_all <- st_coordinates(hwy11_ll) %>%
  as.data.frame() %>%
  dplyr::rename(lon = X, lat = Y)

coords_all <- coords_all %>% arrange(lat)

dists <- distHaversine(
  coords_all[, c("lon", "lat")][-nrow(coords_all), ],
  coords_all[, c("lon", "lat")][-1, ]
)
coords_all$s <- c(0, cumsum(dists))

total_length_km <- max(coords_all$s) / 1000
cat("Approx Highway 11 length (from shapefile) =",
    round(total_length_km, 1), "km\n")

km_list <- c(
  0.0, 3.6, 5.1, 5.6, 7.6, 9.6, 12.7,
  18.1, 22.4, 28.8, 32.4, 32.9,
  38.7, 46.8, 58.8, 59.8, 77.9,
  89.4, 116.5, 117.4, 119.7,
  131.9, 143.8, 147.0,
  179.8, 196.7, 215.5, 218.7, 223.6,
  246.9, 249.7, 250.9, 252.1, 253.8, 254.9,
  256.1, 258.1, 258.8, 258.9,
  259.6, 260.1, 262.8, 264.1, 265.3, 266.9,
  267.9, 272.3, 283.2, 285.6,
  288.8, 290.4, 309.1, 326.4,
  344.7, 393.3
)

max_s  <- max(coords_all$s)
max_km <- max(km_list)

get_location_from_km <- function(km_target, coords_df, max_km, max_s) {
  target_s <- (km_target / max_km) * max_s
  idx <- which.min(abs(coords_df$s - target_s))
  data.frame(
    km_raw = km_target,
    s_m    = coords_df$s[idx],
    lon    = coords_df$lon[idx],
    lat    = coords_df$lat[idx]
  )
}

results <- do.call(
  rbind,
  lapply(km_list, get_location_from_km,
         coords_df = coords_all,
         max_km = max_km,
         max_s  = max_s)
)
results$label <- paste0("P", seq_len(nrow(results)))
cat("Prepared intersection-like sampled points from shapefile.\n")

cons_hwy <- data.frame(
  name = c("Hawarden", "Martensville"),
  lat  = c(51.284722, 52.2890),
  lon  = c(-106.195833, -106.6670)
)

X_base <- results$lat
Y_base <- results$lon

Xc <- cons_hwy$lat
Yc <- cons_hwy$lon

Ybar <- mean(Y_base)
sY   <- sd(Y_base)
if (!is.finite(sY) || sY <= 0) stop("sd(Y_base) is not positive; cannot standardize.")

Y_base_star <- (Y_base - Ybar) / sY
Yc_star     <- (Yc     - Ybar) / sY

inv_std <- function(y_star) Ybar + sY * y_star

K1d <- function(u) dnorm(u)

nw_1d_point <- function(x0, X, Y, h) {
  w   <- K1d((x0 - X) / h)
  den <- sum(w)
  if (den == 0) return(NA_real_)
  sum(w * Y) / den
}

ANW_1d_point <- function(x0, X, Y, Xc, Yc, h, Lambda) {
  w_data <- K1d((x0 - X) / h)
  w_cons <- Lambda * K1d((x0 - Xc) / h)
  den <- sum(w_data) + sum(w_cons)
  if (den == 0) return(NA_real_)
  (sum(w_data * Y) + sum(w_cons * Yc)) / den
}

pred_on_grid_1d <- function(xg, FUN) vapply(xg, FUN, numeric(1))

loocv_mse_nw_1d <- function(h, X, Y) {
  if (h <= 0) return(Inf)
  n   <- length(X)
  se  <- numeric(n)
  for (i in 1:n) {
    y_hat <- nw_1d_point(X[i], X[-i], Y[-i], h)
    se[i] <- (Y[i] - y_hat)^2
  }
  mean(se, na.rm = TRUE)
}

Lx        <- diff(range(X_base))
h_grid_1d <- seq(0.1 * Lx, 0.5 * Lx, length.out = 25)

cv_nw_1d <- sapply(h_grid_1d, function(h) loocv_mse_nw_1d(h, X = X_base, Y = Y_base_star))
h_nw_opt_1d <- h_grid_1d[which.min(cv_nw_1d)]
cat(sprintf("Highway 11 NW optimal h (lat): %.6f degrees [standardized CV]\n", h_nw_opt_1d))

loocv_mse_ANW_1d_star <- function(Lambda, h, X, Ystar, Xc, Ycstar) {
  if (Lambda < 0 || h <= 0) return(Inf)
  n  <- length(X)
  se <- numeric(n)
  for (i in 1:n) {
    y_hat <- ANW_1d_point(
      x0 = X[i],
      X  = X[-i],
      Y  = Ystar[-i],
      Xc = Xc,
      Yc = Ycstar,
      h  = h,
      Lambda = Lambda
    )
    se[i] <- (Ystar[i] - y_hat)^2
  }
  mean(se, na.rm = TRUE)
}

constraint_sq_error_star <- function(Lambda, h, X, Ystar, Xc, Ycstar) {
  y_hatC <- vapply(
    Xc,
    function(x0) ANW_1d_point(x0, X = X, Y = Ystar, Xc = Xc, Yc = Ycstar, h = h, Lambda = Lambda),
    numeric(1)
  )
  sum((y_hatC - Ycstar)^2)
}

alpha <- 1  

loss_ANW_star <- function(h, Lambda, X, Ystar, Xc, Ycstar, alpha) {
  mse_loo <- loocv_mse_ANW_1d_star(Lambda, h, X, Ystar, Xc, Ycstar)
  c_sq    <- constraint_sq_error_star(Lambda, h, X, Ystar, Xc, Ycstar)
  mse_loo + alpha * c_sq
}

Lambda_grid_1d <- seq(0, 120, by = 0.2) 
h_grid_ANW     <- h_grid_1d

loss_mat <- matrix(NA_real_, nrow = length(h_grid_ANW), ncol = length(Lambda_grid_1d))
for (i in seq_along(h_grid_ANW)) {
  for (j in seq_along(Lambda_grid_1d)) {
    loss_mat[i, j] <- loss_ANW_star(
      h      = h_grid_ANW[i],
      Lambda = Lambda_grid_1d[j],
      X      = X_base,
      Ystar  = Y_base_star,
      Xc     = Xc,
      Ycstar = Yc_star,
      alpha  = alpha
    )
  }
}

idx_min <- which(loss_mat == min(loss_mat, na.rm = TRUE), arr.ind = TRUE)
h_ANW_opt_1d  <- h_grid_ANW[idx_min[1, 1]]
Lambda_opt_1d <- Lambda_grid_1d[idx_min[1, 2]]

cat(sprintf("Highway 11 ANW optimal h: %.6f, Lambda: %.2f (alpha=%.2f) [standardized loss]\n",
            h_ANW_opt_1d, Lambda_opt_1d, alpha))

anw_fit_at_X_star <- function(Ym_star) {
  vapply(
    X_base,
    function(x0) ANW_1d_point(x0, X = X_base, Y = Ym_star, Xc = Xc, Yc = Yc_star,
                              h = h_ANW_opt_1d, Lambda = Lambda_opt_1d),
    numeric(1)
  )
}

ds_anw_predict_star <- function(M = 1) {
  if (M <= 0) {
    return(vapply(
      xg,
      function(x0) ANW_1d_point(x0, X = X_base, Y = Y_base_star, Xc = Xc, Yc = Yc_star,
                                h = h_ANW_opt_1d, Lambda = Lambda_opt_1d),
      numeric(1)
    ))
  }
  Ym <- Y_base_star
  for (m in 1:M) {
    fit_m <- anw_fit_at_X_star(Ym)
    Ym <- Y_base_star + (Ym - fit_m)
  }
  vapply(
    xg,
    function(x0) ANW_1d_point(x0, X = X_base, Y = Ym, Xc = Xc, Yc = Yc_star,
                              h = h_ANW_opt_1d, Lambda = Lambda_opt_1d),
    numeric(1)
  )
}

y_hat_DS1_star <- ds_anw_predict_star(M = 1)
y_hat_DS2_star <- ds_anw_predict_star(M = 2)

y_hat_DS1 <- inv_std(y_hat_DS1_star)
y_hat_DS2 <- inv_std(y_hat_DS2_star)

xg <- seq(min(X_base), max(X_base), length.out = 1500)

y_hat_nw_star <- pred_on_grid_1d(xg, function(xx) {
  nw_1d_point(xx, X = X_base, Y = Y_base_star, h = h_nw_opt_1d)
})
y_hat_ANW_star <- pred_on_grid_1d(xg, function(xx) {
  ANW_1d_point(xx, X = X_base, Y = Y_base_star, Xc = Xc, Yc = Yc_star,
               h = h_ANW_opt_1d, Lambda = Lambda_opt_1d)
})

y_hat_nw  <- inv_std(y_hat_nw_star)
y_hat_ANW <- inv_std(y_hat_ANW_star)

baseline_rot <- coords_all %>%
  arrange(lat) %>%
  transmute(lat = lat, lon = lon, method = "Baseline")

get_truth_lon <- function(lat0) {
  idx <- which.min(abs(baseline_rot$lat - lat0))
  baseline_rot$lon[idx]
}
truth_on_grid <- vapply(xg, get_truth_lon, numeric(1))

rmse <- function(truth, pred) sqrt(mean((truth - pred)^2, na.rm = TRUE))

constraint_rmse <- function(Xc, Yc, pred_fun) {
  y_hat_c <- vapply(Xc, pred_fun, numeric(1))
  sqrt(mean((Yc - y_hat_c)^2, na.rm = TRUE))
}

smooth_index <- function(y) {
  d2 <- diff(y, differences = 2)
  mean(d2^2, na.rm = TRUE)
}

phi <- function(z) {
  z <- pmax(z, 0)
  z / (1 + z)
}

rmse_nw   <- rmse(truth_on_grid, y_hat_nw)
rmse_ANW  <- rmse(truth_on_grid, y_hat_ANW)
rmse_DS1  <- rmse(truth_on_grid, y_hat_DS1)
rmse_DS2  <- rmse(truth_on_grid, y_hat_DS2)

pred_nw_lon <- function(x0) inv_std(nw_1d_point(x0, X_base, Y_base_star, h_nw_opt_1d))
pred_anw_lon <- function(x0) inv_std(ANW_1d_point(x0, X_base, Y_base_star, Xc, Yc_star, h_ANW_opt_1d, Lambda_opt_1d))

ds_anw_pred_at_points_lon <- function(X_eval, M = 1) {
  Ym <- Y_base_star
  if (M > 0) {
    for (m in 1:M) {
      fit_m <- anw_fit_at_X_star(Ym)
      Ym <- Y_base_star + (Ym - fit_m)
    }
  }
  inv_std(vapply(
    X_eval,
    function(x0) ANW_1d_point(x0, X = X_base, Y = Ym, Xc = Xc, Yc = Yc_star,
                              h = h_ANW_opt_1d, Lambda = Lambda_opt_1d),
    numeric(1)
  ))
}

ce_nw  <- constraint_rmse(Xc, Yc, pred_nw_lon)
ce_ANW <- constraint_rmse(Xc, Yc, pred_anw_lon)

pred_DS1_at_Xc <- ds_anw_pred_at_points_lon(Xc, M = 1)
pred_DS2_at_Xc <- ds_anw_pred_at_points_lon(Xc, M = 2)
ce_DS1 <- sqrt(mean((Yc - pred_DS1_at_Xc)^2, na.rm = TRUE))
ce_DS2 <- sqrt(mean((Yc - pred_DS2_at_Xc)^2, na.rm = TRUE))

sm_nw   <- smooth_index(y_hat_nw)
sm_ANW  <- smooth_index(y_hat_ANW)
sm_DS1  <- smooth_index(y_hat_DS1)
sm_DS2  <- smooth_index(y_hat_DS2)

X_aug <- c(X_base, Xc)
Y_aug <- c(Y_base, Yc)

B <- bs(X_aug)
fit_bs <- lm(Y_aug ~ B)
bs_y  <- predict(fit_bs, newdata = list(B = bs(xg)))
bs_Yc <- predict(fit_bs, newdata = list(B = bs(Xc)))

fit_ps <- gam(Y_aug ~ s(X_aug, bs = "ps", k = 20, m = c(2, 2)), method = "REML")
ps_y  <- predict(fit_ps, newdata = data.frame(X_aug = xg))
ps_Yc <- predict(fit_ps, newdata = data.frame(X_aug = Xc))

fit_lo <- loess(Y_aug ~ X_aug)
lo_y  <- predict(fit_lo, newdata = data.frame(X_aug = xg))
lo_Yc <- predict(fit_lo, newdata = data.frame(X_aug = Xc))

fit_ss <- smooth.spline(X_aug, Y_aug, spar = 0.7)
ss_y  <- predict(fit_ss, xg)$y
ss_Yc <- predict(fit_ss, Xc)$y

fit_gpr <- gpkm(X_aug, Y_aug)
gpr_y  <- fit_gpr$predict(xg)
gpr_Yc <- fit_gpr$predict(Xc)

fit_atps <- gam(
  Y_aug ~ s(X_aug, bs = "ts", k = 20) + s(X_aug, bs = "ts", k = 20, m = 2),
  method = "REML"
)
atps_y  <- predict(fit_atps, newdata = data.frame(X_aug = xg))
atps_Yc <- predict(fit_atps, newdata = data.frame(X_aug = Xc))

tau_c <- 0.1
tau_g <- 0.05

sm_all <- c(sm_nw, sm_ANW, sm_DS1, sm_DS2,
            smooth_index(bs_y), smooth_index(ps_y), smooth_index(lo_y),
            smooth_index(ss_y), smooth_index(gpr_y), smooth_index(atps_y))
tau_s <- median(sm_all)

w_c <- 0.5; w_g <- 0.3; w_s <- 0.2
css_fun <- function(ce, g_rmse, sm) w_c * phi(ce / tau_c) + w_g * phi(g_rmse / tau_g) + w_s * phi(sm / tau_s)

tab <- data.frame(
  Model = c("NW", "ANW", "DS-ANW (M=1)", "DS-ANW (M=2)", "BS", "PS", "LOESS", "SS", "GPR", "ATPS"),
  RMSE = c(
    rmse_nw, rmse_ANW, rmse_DS1, rmse_DS2,
    rmse(truth_on_grid, bs_y),
    rmse(truth_on_grid, ps_y),
    rmse(truth_on_grid, lo_y),
    rmse(truth_on_grid, ss_y),
    rmse(truth_on_grid, gpr_y),
    rmse(truth_on_grid, atps_y)
  ),
  ConstraintErr = c(
    ce_nw, ce_ANW, ce_DS1, ce_DS2,
    sqrt(mean((Yc - bs_Yc)^2)),
    sqrt(mean((Yc - ps_Yc)^2)),
    sqrt(mean((Yc - lo_Yc)^2)),
    sqrt(mean((Yc - ss_Yc)^2)),
    sqrt(mean((Yc - gpr_Yc)^2)),
    sqrt(mean((Yc - atps_Yc)^2))
  ),
  Smoothness = c(
    sm_nw, sm_ANW, sm_DS1, sm_DS2,
    smooth_index(bs_y),
    smooth_index(ps_y),
    smooth_index(lo_y),
    smooth_index(ss_y),
    smooth_index(gpr_y),
    smooth_index(atps_y)
  )
)
tab$CSS <- mapply(css_fun, tab$ConstraintErr, tab$RMSE, tab$Smoothness)
print(tab)

curves_rot <- bind_rows(
  baseline_rot,
  data.frame(lat = xg, lon = y_hat_nw,  method = "NW"),
  data.frame(lat = xg, lon = y_hat_ANW, method = "ANW"),
  data.frame(lat = xg, lon = y_hat_DS1, method = "DS-ANW (M=1)"),
  data.frame(lat = xg, lon = y_hat_DS2, method = "DS-ANW (M=2)"),
  data.frame(lat = xg, lon = bs_y,      method = "BS"),
  data.frame(lat = xg, lon = ps_y,      method = "PS"),
  data.frame(lat = xg, lon = lo_y,      method = "LOESS"),
  data.frame(lat = xg, lon = ss_y,      method = "SS"),
  data.frame(lat = xg, lon = gpr_y,     method = "GPR"),
  data.frame(lat = xg, lon = atps_y,    method = "ATPS")
)

p_rot <- ggplot() +
  geom_path(
    data = curves_rot,
    aes(x = lat, y = lon, color = method, linetype = method, linewidth = method)
  ) +
  geom_point(
    data = results,
    aes(x = lat, y = lon),
    shape = 21, size = 1.6, stroke = 0.35,
    fill = "white", color = "black"
  ) +
  geom_point(
    data = cons_hwy,
    aes(x = lat, y = lon),
    shape = 21, size = 3.0, stroke = 0.8,
    fill = "gold", color = "black"
  ) +
  geom_text_repel(
    data = cons_hwy,
    aes(x = lat, y = lon, label = name),
    size = 3.0,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.15
  ) +
  coord_fixed() +
  labs(
    title = "Highway 11 (rotated coords)",
    subtitle = sprintf("NW: h=%.4f | ANW/DS-ANW: h=%.4f, Lambda=%.1f | alpha=%.1f",
                       h_nw_opt_1d, h_ANW_opt_1d, Lambda_opt_1d, alpha),
    x = "Latitude",
    y = "Longitude",
    color = "Method", linetype = "Method", linewidth = "Method"
  ) +
  scale_color_manual(values = c(
    "Baseline"       = "black",
    "NW"             = "green4",
    "ANW"            = "firebrick",
    "DS-ANW (M=1)"   = "darkorange3",
    "DS-ANW (M=2)"   = "darkgreen",
    "BS"             = "blue",
    "PS"             = "red",
    "LOESS"          = "orange",
    "SS"             = "purple2",
    "GPR"            = "skyblue",
    "ATPS"           = "goldenrod2"
  )) +
  scale_linetype_manual(values = c(
    "Baseline"       = "solid",
    "NW"             = "dashed",
    "ANW"            = "solid",
    "DS-ANW (M=1)"   = "dashed",
    "DS-ANW (M=2)"   = "dashed",
    "BS"             = "dashed",
    "PS"             = "dashed",
    "LOESS"          = "dashed",
    "SS"             = "dashed",
    "GPR"            = "dashed",
    "ATPS"           = "dashed"
  )) +
  scale_linewidth_manual(values = c(
    "Baseline"       = 0.9,
    "NW"             = 0.8,
    "ANW"            = 1.2,
    "DS-ANW (M=1)"   = 1.0,
    "DS-ANW (M=2)"   = 1.0,
    "BS"             = 0.8,
    "PS"             = 0.8,
    "LOESS"          = 0.8,
    "SS"             = 0.8,
    "GPR"            = 0.8,
    "ATPS"           = 0.8
  )) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.box = "vertical")

print(p_rot)

if (!dir.exists("figs")) dir.create("figs")
pdf("figs/fig_highway11_rotated_allmethods_with_ds.pdf", width=9, height=6)
print(p_rot)
dev.off()
cat("Saved: figs/fig_highway11_rotated_allmethods_with_ds.pdf\n")
