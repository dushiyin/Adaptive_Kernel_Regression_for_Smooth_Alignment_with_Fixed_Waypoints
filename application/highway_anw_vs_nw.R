suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(geosphere)
  library(ggplot2)
  library(ggrepel)
})

shp_path <- "D:/Judy File/UBCO/research/Paper/data/Saskatchewan/HIGHWAY_OFFICIAL.shp"
roads_sf <- st_read(shp_path, quiet = TRUE)

hwy11_sf <- roads_sf %>%
  filter(RTNUMBER1 == 11 | RTNUMBER2 == 11 | RTNUMBER3 == 11 | RTNUMBER4 == 11)

hwy11_ll <- st_transform(hwy11_sf, 4326)  

coords_all <- st_coordinates(hwy11_ll) %>%
  as.data.frame() %>%
  dplyr::rename(lon = X, lat = Y) %>%
  arrange(lat)

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

write.csv(results, "Highway11_major_intersections_coords.csv", row.names = FALSE)
cat("Saved: Highway11_major_intersections_coords.csv\n")

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

Y_star  <- (Y_base - Ybar) / sY
Yc_star <- (Yc     - Ybar) / sY
inv_std <- function(y_star) Ybar + sY * y_star

K1d <- function(u) dnorm(u)

nw_point_star <- function(x0, X, Ystar, h) {
  w   <- K1d((x0 - X) / h)
  den <- sum(w)
  if (den == 0) return(NA_real_)
  sum(w * Ystar) / den
}

anw_point_star <- function(x0, X, Ystar, Xc, Yc_star, h, Lambda) {
  w_data <- K1d((x0 - X) / h)
  w_cons <- Lambda * K1d((x0 - Xc) / h)
  den <- sum(w_data) + sum(w_cons)
  if (den == 0) return(NA_real_)
  (sum(w_data * Ystar) + sum(w_cons * Yc_star)) / den
}

pred_on_grid <- function(xg, FUN) vapply(xg, FUN, numeric(1))

loocv_mse_nw_star <- function(h, X, Ystar) {
  if (h <= 0) return(Inf)
  n  <- length(X)
  se <- numeric(n)
  for (i in 1:n) {
    y_hat <- nw_point_star(X[i], X[-i], Ystar[-i], h)
    se[i] <- (Ystar[i] - y_hat)^2
  }
  mean(se, na.rm = TRUE)
}

Lx        <- diff(range(X_base))
h_grid_1d <- seq(0.1 * Lx, 0.5 * Lx, length.out = 25)

cv_nw <- sapply(h_grid_1d, function(h) loocv_mse_nw_star(h, X = X_base, Ystar = Y_star))
h_nw_opt <- h_grid_1d[which.min(cv_nw)]
cat(sprintf("NW optimal h = %.6f (standardized LOOCV)\n", h_nw_opt))

loocv_mse_anw_star <- function(Lambda, h, X, Ystar, Xc, Yc_star) {
  if (Lambda < 0 || h <= 0) return(Inf)
  n  <- length(X)
  se <- numeric(n)
  for (i in 1:n) {
    y_hat <- anw_point_star(
      x0 = X[i],
      X  = X[-i],
      Ystar = Ystar[-i],
      Xc = Xc,
      Yc_star = Yc_star,
      h  = h,
      Lambda = Lambda
    )
    se[i] <- (Ystar[i] - y_hat)^2
  }
  mean(se, na.rm = TRUE)
}

constraint_sq_star <- function(Lambda, h, X, Ystar, Xc, Yc_star) {
  y_hatC <- vapply(
    Xc,
    function(x0) anw_point_star(x0, X, Ystar, Xc, Yc_star, h, Lambda),
    numeric(1)
  )
  sum((y_hatC - Yc_star)^2)
}

alpha <- 1

loss_anw_star <- function(h, Lambda, X, Ystar, Xc, Yc_star, alpha) {
  loocv_mse_anw_star(Lambda, h, X, Ystar, Xc, Yc_star) +
    alpha * constraint_sq_star(Lambda, h, X, Ystar, Xc, Yc_star)
}

Lambda_grid <- seq(0, 120, by = 0.2)

loss_mat <- matrix(NA_real_, nrow = length(h_grid_1d), ncol = length(Lambda_grid))
for (i in seq_along(h_grid_1d)) {
  for (j in seq_along(Lambda_grid)) {
    loss_mat[i, j] <- loss_anw_star(
      h = h_grid_1d[i],
      Lambda = Lambda_grid[j],
      X = X_base,
      Ystar = Y_star,
      Xc = Xc,
      Yc_star = Yc_star,
      alpha = alpha
    )
  }
}
idx <- which(loss_mat == min(loss_mat, na.rm = TRUE), arr.ind = TRUE)
h_anw_opt  <- h_grid_1d[idx[1, 1]]
Lambda_opt <- Lambda_grid[idx[1, 2]]
cat(sprintf("ANW optimal: h = %.6f, Lambda = %.2f (standardized Loss)\n", h_anw_opt, Lambda_opt))

anw_pred_vec_star <- function(X_train, Y_train_star, X_eval, Xc, Yc_star, h, Lambda) {
  vapply(
    X_eval,
    function(x0) anw_point_star(x0, X_train, Y_train_star, Xc, Yc_star, h, Lambda),
    numeric(1)
  )
}

ds_anw_pred_vec_star <- function(X_train, Y0_star, X_eval, Xc, Yc_star, h, Lambda, M = 1) {
  if (M <= 0) {
    return(anw_pred_vec_star(X_train, Y0_star, X_eval, Xc, Yc_star, h, Lambda))
  }
  Ym <- Y0_star
  for (m in 1:M) {
    fit_train <- anw_pred_vec_star(X_train, Ym, X_train, Xc, Yc_star, h, Lambda)
    Ym <- Y0_star + (Ym - fit_train)
  }
  anw_pred_vec_star(X_train, Ym, X_eval, Xc, Yc_star, h, Lambda)
}

xg <- seq(min(X_base), max(X_base), length.out = 1500)

yhat_nw_star  <- pred_on_grid(xg, function(xx) nw_point_star(xx, X_base, Y_star, h_nw_opt))
yhat_anw_star <- pred_on_grid(xg, function(xx) anw_point_star(xx, X_base, Y_star, Xc, Yc_star, h_anw_opt, Lambda_opt))

yhat_ds1_star <- ds_anw_pred_vec_star(
  X_train = X_base, Y0_star = Y_star, X_eval = xg,
  Xc = Xc, Yc_star = Yc_star,
  h = h_anw_opt, Lambda = Lambda_opt, M = 1
)
yhat_ds2_star <- ds_anw_pred_vec_star(
  X_train = X_base, Y0_star = Y_star, X_eval = xg,
  Xc = Xc, Yc_star = Yc_star,
  h = h_anw_opt, Lambda = Lambda_opt, M = 2
)

yhat_nw   <- inv_std(yhat_nw_star)
yhat_anw  <- inv_std(yhat_anw_star)
yhat_ds1  <- inv_std(yhat_ds1_star)
yhat_ds2  <- inv_std(yhat_ds2_star)

baseline_rot <- coords_all %>%
  transmute(lat = lat, lon = lon, method = "Baseline")

curves_rot <- bind_rows(
  baseline_rot,
  data.frame(lat = xg, lon = yhat_nw,  method = "NW"),
  data.frame(lat = xg, lon = yhat_anw, method = "ANW"),
  data.frame(lat = xg, lon = yhat_ds1, method = "DS-ANW (M=1)"),
  data.frame(lat = xg, lon = yhat_ds2, method = "DS-ANW (M=2)")
)

p_rot <- ggplot() +
  geom_path(
    data = curves_rot,
    aes(x = lat, y = lon, color = method, linetype = method, linewidth = method)
  ) +
  geom_point(
    data = results,
    aes(x = lat, y = lon),
    shape = 21, size = 1.8, stroke = 0.4,
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
    size = 3.1,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.15
  ) +
  coord_fixed() +
  labs(
    title = "Highway 11 (rotated coords): NW vs ANW vs DS-ANW",
    subtitle = sprintf(
      "NW: h=%.4f | ANW/DS-ANW: h=%.4f, Lambda=%.1f | alpha=%.1f",
      h_nw_opt, h_anw_opt, Lambda_opt, alpha
    ),
    x = "Latitude",
    y = "Longitude",
    color = "Method", linetype = "Method", linewidth = "Method"
  ) +
  scale_color_manual(values = c(
    "Baseline"     = "black",
    "NW"           = "blue",
    "ANW"          = "firebrick",
    "DS-ANW (M=1)" = "orange",
    "DS-ANW (M=2)" = "darkgreen"
  )) +
  scale_linetype_manual(values = c(
    "Baseline"     = "solid",
    "NW"           = "dashed",
    "ANW"          = "solid",
    "DS-ANW (M=1)" = "dashed",
    "DS-ANW (M=2)" = "dashed"
  )) +
  scale_linewidth_manual(values = c(
    "Baseline"     = 0.9,
    "NW"           = 0.8,
    "ANW"          = 1.2,
    "DS-ANW (M=1)" = 1.0,
    "DS-ANW (M=2)" = 1.0
  )) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.box = "vertical")

print(p_rot)

if (!dir.exists("figs")) dir.create("figs")
pdf("figs/fig_highway11_nw_vs_anw_vs_ds.pdf", width = 9, height = 6)
print(p_rot)
dev.off()
cat("Saved: figs/fig_highway11_nw_vs_anw_vs_ds.pdf\n")
