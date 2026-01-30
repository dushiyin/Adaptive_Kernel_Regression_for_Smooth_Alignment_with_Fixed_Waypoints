library(dplyr)
library(ggplot2)
library(maps)
library(ggrepel)

library(splines)
library(mgcv)
library(GauPro)

haversine_m <- function(lon1, lat1, lon2, lat2) {
  R <- 6371000
  toRad <- pi/180
  
  dlat <- (lat2 - lat1) * toRad
  dlon <- (lon2 - lon1) * toRad
  
  a <- sin(dlat/2)^2 +
    cos(lat1*toRad) * cos(lat2*toRad) * sin(dlon/2)^2
  
  2 * R * asin(pmin(1, sqrt(a)))
}

lonlat_to_xy_segment <- function(lon, lat) {
  lat0 <- mean(lat)
  toRad <- pi/180
  R <- 6371000
  
  x <- (lon * toRad) * R * cos(lat0 * toRad)
  y <- (lat * toRad) * R
  
  cbind(x, y)
}

polyline_cumlen_m <- function(lon, lat) {
  c(
    0,
    cumsum(
      haversine_m(
        lon[-length(lon)], lat[-length(lat)],
        lon[-1],            lat[-1]
      )
    )
  )
}

project_point_onto_polyline <- function(lon, lat, line_lon, line_lat, s_accum) {
  XY <- lonlat_to_xy_segment(line_lon, line_lat)
  P  <- lonlat_to_xy_segment(lon, lat)
  
  best <- list(dist2 = Inf, s = NA_real_, lon = NA_real_, lat = NA_real_)
  
  for (i in 1:(nrow(XY) - 1)) {
    A <- XY[i, ]
    B <- XY[i + 1, ]
    v <- B - A
    w <- P - A
    
    t <- sum(w * v) / sum(v * v)
    t <- max(0, min(1, t))
    
    lon_q <- line_lon[i] + t * (line_lon[i + 1] - line_lon[i])
    lat_q <- line_lat[i] + t * (line_lat[i + 1] - line_lat[i])
    
    Q  <- A + t * v
    d2 <- sum((P - Q)^2)
    
    if (d2 < best$dist2) {
      seg_len <- haversine_m(
        line_lon[i],     line_lat[i],
        line_lon[i + 1], line_lat[i + 1]
      )
      
      s_here <- s_accum[i] + t * seg_len
      
      best <- list(
        dist2 = d2,
        s     = s_here,
        lon   = lon_q,
        lat   = lat_q
      )
    }
  }
  
  best
}

K1d <- function(u) dnorm(u)

gtfs_dir <- "D:/Judy File/UBCO/research/Paper/data/viarail"

shapes     <- read.csv(file.path(gtfs_dir, "shapes.txt"))
stops      <- read.csv(file.path(gtfs_dir, "stops.txt"))
trips      <- read.csv(file.path(gtfs_dir, "trips.txt"))
stop_times <- read.csv(file.path(gtfs_dir, "stop_times.txt"))

shapes$shape_id    <- as.character(shapes$shape_id)
trips$shape_id     <- as.character(trips$shape_id)
trips$trip_id      <- as.character(trips$trip_id)
stop_times$trip_id <- as.character(stop_times$trip_id)
stop_times$stop_id <- as.character(stop_times$stop_id)
stops$stop_id      <- as.character(stops$stop_id)


shape111 <- shapes %>%
  filter(shape_id == "111") %>%
  arrange(as.numeric(shape_pt_sequence))

stopifnot(nrow(shape111) > 0)

if ("shape_dist_traveled" %in% names(shape111) &&
    all(is.finite(shape111$shape_dist_traveled))) {
  s_baseline <- as.numeric(shape111$shape_dist_traveled)
  s_baseline <- s_baseline - min(s_baseline)
} else {
  s_baseline <- polyline_cumlen_m(
    shape111$shape_pt_lon,
    shape111$shape_pt_lat
  )
}

shape111$s <- s_baseline


trips111 <- trips %>% filter(shape_id == "111")
stopifnot(nrow(trips111) > 0)

main_trip_id <- trips111$trip_id[1]

stops_on_line <- stop_times %>%
  filter(trip_id == main_trip_id) %>%
  arrange(as.numeric(stop_sequence)) %>%
  select(stop_id, stop_sequence) %>%
  distinct(stop_id, .keep_all = TRUE)

stops_line <- stops_on_line %>%
  left_join(stops, by = "stop_id") %>%
  filter(is.finite(stop_lon), is.finite(stop_lat))

proj_list <- lapply(
  1:nrow(stops_line),
  function(i) {
    project_point_onto_polyline(
      lon      = stops_line$stop_lon[i],
      lat      = stops_line$stop_lat[i],
      line_lon = shape111$shape_pt_lon,
      line_lat = shape111$shape_pt_lat,
      s_accum  = shape111$s
    )
  }
)

stops_line$s <- vapply(proj_list, `[[`, numeric(1), "s")
stops_line   <- stops_line %>% arrange(s)

banff <- data.frame(
  name = "Banff",
  lon  = -115.5781,
  lat  = 51.1784
)

projC <- project_point_onto_polyline(
  lon      = banff$lon,
  lat      = banff$lat,
  line_lon = shape111$shape_pt_lon,
  line_lat = shape111$shape_pt_lat,
  s_accum  = shape111$s
)

cons <- data.frame(
  name      = "Banff",
  lon       = banff$lon,
  lat       = banff$lat,
  sC        = projC$s,
  label_lon = -115.0,
  label_lat = 50.2
)

stations9 <- data.frame(
  name=c("Vancouver","Kamloops","Jasper","Edmonton","Saskatoon",
         "Winnipeg","Sioux Lookout","Sudbury Jct.","Toronto"),
  lon=c(-123.0971,-120.3570,-118.0790,-113.4745,-106.6689,
        -97.1384,-91.9163,-80.8575,-79.3807),
  lat=c(49.2734,50.7168,52.8734,53.5426,52.1345,
        49.8951,50.0973,46.5382,43.6457),
  label_lon = c(-120.9, -122.7, -119.1, -113.0, -105.8,
                -97.9, -91.1, -84.0, -81.7),
  label_lat = c(48.2,   51.9,   54.2,   54.5,   53.2,
                51.1,   51.2,   46.3,  43.7)
)

X_data <- stops_line$stop_lon
Y_data <- stops_line$stop_lat

Xc <- cons$lon
Yc <- cons$lat

X_data <- c(X_data, Xc)
Y_data <- c(Y_data, Yc)

nw_1d_point <- function(x0, X, Y, h) {
  w   <- K1d((x0 - X) / h)
  num <- sum(w * Y)
  den <- sum(w)
  if (den == 0) return(NA_real_) else num / den
}

anw_1d_point <- function(x0, X, Y, Xc, Yc, h, lambda) {
  w_data <- K1d((x0 - X) / h)
  w_cons <- lambda * K1d((x0 - Xc) / h)
  num <- sum(w_data * Y) + sum(w_cons * Yc)
  den <- sum(w_data) + sum(w_cons)
  if (den == 0) return(NA_real_) else num / den
}

pred_on_grid_1d <- function(xg, FUN) vapply(xg, FUN, numeric(1))

anw_pred_vec_std <- function(x_train, y_train_std, x_eval, Xc, Yc_std, h, lambda) {
  vapply(
    x_eval,
    function(xx) anw_1d_point(xx, X = x_train, Y = y_train_std, Xc = Xc, Yc = Yc_std, h = h, lambda = lambda),
    numeric(1)
  )
}

ds_anw_pred_vec_std <- function(x_train, y0_std, x_eval, Xc, Yc_std, h, lambda, M = 1) {
  if (M <= 0) {
    return(anw_pred_vec_std(x_train, y0_std, x_eval, Xc, Yc_std, h, lambda))
  }
  ym <- y0_std
  for (m in 1:M) {
    fit_train <- anw_pred_vec_std(x_train, ym, x_train, Xc, Yc_std, h, lambda)
    ym <- y0_std + (ym - fit_train)  
  }
  anw_pred_vec_std(x_train, ym, x_eval, Xc, Yc_std, h, lambda)
}

Y_bar <- mean(Y_data)
sY_sq <- mean((Y_data - Y_bar)^2)
sY    <- sqrt(sY_sq)

Y_std  <- (Y_data - Y_bar) / sY
Yc_std <- (Yc        - Y_bar) / sY

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

Lx        <- diff(range(X_data))
h_grid_1d <- seq(0.02 * Lx, 0.3 * Lx, length.out = 25)

cv_nw_1d <- sapply(
  h_grid_1d,
  function(h) loocv_mse_nw_1d(h, X = X_data, Y = Y_std)
)
h_nw_opt_1d <- h_grid_1d[which.min(cv_nw_1d)]
cat(sprintf("1D NW optimal h: %.4f (longitude degrees)\n", h_nw_opt_1d))


xg <- seq(
  min(c(X_data, Xc)),
  max(c(X_data, Xc)),
  length.out = 1500
)


loocv_mse_anw_1d <- function(lambda, h, X, Y, Xc, Yc) {
  if (lambda < 0 || h <= 0) return(Inf)
  n  <- length(X)
  se <- numeric(n)
  for (i in 1:n) {
    y_hat <- anw_1d_point(
      x0 = X[i],
      X  = X[-i],
      Y  = Y[-i],
      Xc = Xc,
      Yc = Yc,
      h  = h,
      lambda = lambda
    )
    se[i] <- (Y[i] - y_hat)^2
  }
  mean(se, na.rm = TRUE)
}

constraint_sq_error <- function(lambda, h, X, Y, Xc, Yc) {
  y_hatC <- vapply(
    Xc,
    function(x0) {
      anw_1d_point(
        x0,
        X  = X,
        Y  = Y,
        Xc = Xc,
        Yc = Yc,
        h  = h,
        lambda = lambda
      )
    },
    numeric(1)
  )
  sum((y_hatC - Yc)^2)
}

alpha <- 1  

loss_anw <- function(h, lambda, X, Y, Xc, Yc, alpha) {
  if (h <= 0 || lambda < 0) return(Inf)
  mse_loo <- loocv_mse_anw_1d(lambda, h, X, Y, Xc, Yc)
  c_sq    <- constraint_sq_error(lambda, h, X, Y, Xc, Yc)
  mse_loo + alpha * c_sq
}

lambda_grid_1d <- seq(0, 50, by = 0.01)
h_grid_anw     <- h_grid_1d

loss_mat <- matrix(NA_real_, nrow = length(h_grid_anw), ncol = length(lambda_grid_1d))

for (i in seq_along(h_grid_anw)) {
  for (j in seq_along(lambda_grid_1d)) {
    loss_mat[i, j] <- loss_anw(
      h      = h_grid_anw[i],
      lambda = lambda_grid_1d[j],
      X      = X_data,
      Y      = Y_std,
      Xc     = Xc,
      Yc     = Yc_std,
      alpha  = alpha
    )
  }
}

idx_min <- which(loss_mat == min(loss_mat, na.rm = TRUE), arr.ind = TRUE)
h_anw_opt_1d  <- h_grid_anw[idx_min[1, 1]]
lambda_opt_1d <- lambda_grid_1d[idx_min[1, 2]]

cat(sprintf("1D ANW optimal h: %.4f, optimal lambda: %.2f\n",
            h_anw_opt_1d, lambda_opt_1d))

y_hat_nw_1d <- pred_on_grid_1d(xg, function(xx) {
  nw_1d_point(xx, X = X_data, Y = Y_data, h = h_nw_opt_1d)
})

y_hat_anw_1d <- pred_on_grid_1d(xg, function(xx) {
  anw_1d_point(
    xx,
    X  = X_data,
    Y  = Y_data,
    Xc = Xc,
    Yc = Yc,
    h  = h_anw_opt_1d,
    lambda = lambda_opt_1d
  )
})

ds1_star <- ds_anw_pred_vec_std(
  x_train = X_data, y0_std = Y_std, x_eval = xg,
  Xc = Xc, Yc_std = Yc_std,
  h = h_anw_opt_1d, lambda = lambda_opt_1d,
  M = 1
)
ds2_star <- ds_anw_pred_vec_std(
  x_train = X_data, y0_std = Y_std, x_eval = xg,
  Xc = Xc, Yc_std = Yc_std,
  h = h_anw_opt_1d, lambda = lambda_opt_1d,
  M = 2
)

y_hat_ds1_1d <- Y_bar + sY * ds1_star
y_hat_ds2_1d <- Y_bar + sY * ds2_star

B <- bs(X_data)
fit_bs <- lm(Y_data ~ B)
B_new  <- bs(xg)
bs_y   <- predict(fit_bs, newdata = list(B = B_new))
bs_Yc  <- predict(fit_bs, newdata = list(B = bs(Xc)))

fit_ps <- gam(Y_data ~ s(X_data, bs = "ps", k = 20, m = c(2, 2)),
              method = "REML")
ps_y  <- predict(fit_ps, newdata = data.frame(X_data = xg))
ps_Yc <- predict(fit_ps, newdata = data.frame(X_data = Xc))

fit_loess <- loess(Y_data ~ X_data)
LOESS_y   <- predict(fit_loess, data.frame(X_data = xg))
LOESS_Yc  <- predict(fit_loess, data.frame(X_data = Xc))

fit_ss <- smooth.spline(X_data, Y_data)
smooth_spline    <- predict(fit_ss, xg)$y
smooth_spline_Yc <- predict(fit_ss, Xc)$y

fit_gpr <- gpkm(X_data, Y_data)
GPR_y   <- fit_gpr$predict(xg)
GPR_Yc  <- fit_gpr$predict(Xc)

fit_atps <- gam(Y_data ~ s(X_data, bs = "ts", k = 20) +
                  s(X_data, bs = "ts", k = 20, m = 2),
                method = "REML")
ATPS_y  <- predict(fit_atps, data.frame(X_data = xg))
ATPS_Yc <- predict(fit_atps, data.frame(X_data = Xc))

baseline_1d <- data.frame(
  lon    = shape111$shape_pt_lon,
  lat    = shape111$shape_pt_lat,
  method = "Baseline"
)

smooth_nw_1d <- data.frame(lon = xg, lat = y_hat_nw_1d, method = "NW")
smooth_anw_1d <- data.frame(lon = xg, lat = y_hat_anw_1d, method = "ANW")
smooth_ds1_1d <- data.frame(lon = xg, lat = y_hat_ds1_1d, method = "DS-ANW (M=1)")
smooth_ds2_1d <- data.frame(lon = xg, lat = y_hat_ds2_1d, method = "DS-ANW (M=2)")

bsplines <- data.frame(lon = xg, lat = bs_y, method = "BS")
psplines <- data.frame(lon = xg, lat = ps_y, method = "PS")
LOESS <- data.frame(lon = xg, lat = LOESS_y, method = "LOESS")
smoothsplines <- data.frame(lon = xg, lat = smooth_spline, method = "SS")
GPR <- data.frame(lon = xg, lat = GPR_y, method = "GPR")
ATPS <- data.frame(lon = xg, lat = ATPS_y, method = "ATPS")

curves_1d <- bind_rows(
  baseline_1d,
  smooth_nw_1d, smooth_anw_1d, smooth_ds1_1d, smooth_ds2_1d,
  bsplines, psplines, LOESS, smoothsplines, GPR, ATPS
)

world <- map_data("world")

x_rng <- range(c(shape111$shape_pt_lon, cons$lon, stations9$lon), na.rm = TRUE)
y_rng <- range(c(shape111$shape_pt_lat, cons$lat, stations9$lat), na.rm = TRUE)

pad_y <- diff(y_rng) * 0.1
x_lim <- c(-123, -80)

p1d <- ggplot() +
  geom_polygon(
    data = subset(world, region %in% c("Canada", "USA")),
    aes(long, lat, group = group),
    fill  = "gray95", color = "gray70", linewidth = 0.3
  ) +
  geom_path(
    data = curves_1d,
    aes(x = lon, y = lat, color = method, linetype = method, linewidth = method)
  ) +
  geom_point(
    data = stops_line,
    aes(x = stop_lon, y = stop_lat),
    shape = 21, size = 1.5, stroke = 0.3,
    fill = "lightblue", color = "black", alpha = 0.6
  ) +
  geom_point(
    data = cons,
    aes(x = lon, y = lat),
    shape = 21, size = 3.5, stroke = 1,
    fill = "gold", color = "black"
  ) +
  geom_text(
    data = cons,
    aes(x = label_lon, y = label_lat, label = name),
    fontface = "bold",
    size = 4
  ) +
  geom_point(
    data = stations9,
    aes(x = lon, y = lat),
    shape = 21, size = 3, stroke = 0.8,
    fill = "lightblue", color = "black"
  ) +
  geom_text(
    data = stations9,
    aes(x = label_lon, y = label_lat, label = name),
    size = 3.2,
    fontface = "bold"
  ) +
  coord_fixed(
    xlim = x_lim,
    ylim = y_rng + c(-pad_y, pad_y)
  ) +
  labs(
    title = "Railway",
    subtitle = sprintf(
      "NW: h = %.4f | ANW/DS-ANW: h = %.4f, Lambda = %.2f | alpha = %.1f",
      h_nw_opt_1d, h_anw_opt_1d, lambda_opt_1d, alpha
    ),
    x = "Longitude", y = "Latitude",
    color = "Method", linetype = "Method", linewidth = "Method"
  ) +
  scale_color_manual(values = c(
    "Baseline"     = "black",
    "NW"           = "green",
    "ANW"          = "firebrick",
    "DS-ANW (M=1)" = "darkorange3",
    "DS-ANW (M=2)" = "darkgreen",
    "BS"           = "blue",
    "PS"           = "red",
    "LOESS"        = "orange",
    "SS"           = "purple2",
    "GPR"          = "skyblue",
    "ATPS"         = "gold"
  )) +
  scale_linetype_manual(values = c(
    "Baseline"     = "solid",
    "NW"           = "dashed",
    "ANW"          = "dashed",
    "DS-ANW (M=1)" = "dashed",
    "DS-ANW (M=2)" = "dashed",
    "BS"           = "dashed",
    "PS"           = "dashed",
    "LOESS"        = "dashed",
    "SS"           = "dashed",
    "GPR"          = "dashed",
    "ATPS"         = "dashed"
  )) +
  scale_linewidth_manual(values = c(
    "Baseline"     = 0.8,
    "NW"           = 0.8,
    "ANW"          = 1.2,
    "DS-ANW (M=1)" = 1.1,
    "DS-ANW (M=2)" = 1.1,
    "BS"           = 0.8,
    "PS"           = 0.8,
    "LOESS"        = 0.8,
    "SS"           = 0.8,
    "GPR"          = 0.8,
    "ATPS"         = 0.8
  )) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box      = "vertical"
  )

pdf("figs/fig_railway_1d_allmethods_with_ds.pdf", width = 9, height = 5)
print(p1d)
dev.off()

rmse <- function(truth, pred) sqrt(mean((truth - pred)^2, na.rm = TRUE))

constraint_rmse <- function(xC, yC, f_pred) {
  y_hatC <- vapply(xC, f_pred, numeric(1))
  sqrt(mean((y_hatC - yC)^2, na.rm = TRUE))
}

smooth_index <- function(y) {
  d2 <- diff(y, differences = 2)
  mean(d2^2, na.rm = TRUE)
}

phi <- function(x) pmax(0, x - 1)

baseline_sorted <- baseline_1d[order(baseline_1d$lon), ]
truth_on_grid <- approx(
  x = baseline_sorted$lon,
  y = baseline_sorted$lat,
  xout = xg,
  rule = 2
)$y

rmse_nw   <- rmse(truth_on_grid, y_hat_nw_1d)
rmse_anw  <- rmse(truth_on_grid, y_hat_anw_1d)
rmse_ds1  <- rmse(truth_on_grid, y_hat_ds1_1d)
rmse_ds2  <- rmse(truth_on_grid, y_hat_ds2_1d)
rmse_bs   <- rmse(truth_on_grid, bs_y)
rmse_ps   <- rmse(truth_on_grid, ps_y)
rmse_lo   <- rmse(truth_on_grid, LOESS_y)
rmse_ss   <- rmse(truth_on_grid, smooth_spline)
rmse_gpr  <- rmse(truth_on_grid, GPR_y)
rmse_atps <- rmse(truth_on_grid, ATPS_y)

ce_nw <- constraint_rmse(
  Xc, Yc,
  function(x0) nw_1d_point(x0, X_data, Y_data, h_nw_opt_1d)
)

ce_anw <- constraint_rmse(
  Xc, Yc,
  function(x0) anw_1d_point(x0, X_data, Y_data, Xc, Yc, h_anw_opt_1d, lambda_opt_1d)
)

ds1_c_star <- ds_anw_pred_vec_std(
  x_train = X_data, y0_std = Y_std, x_eval = Xc,
  Xc = Xc, Yc_std = Yc_std,
  h = h_anw_opt_1d, lambda = lambda_opt_1d,
  M = 1
)
ds2_c_star <- ds_anw_pred_vec_std(
  x_train = X_data, y0_std = Y_std, x_eval = Xc,
  Xc = Xc, Yc_std = Yc_std,
  h = h_anw_opt_1d, lambda = lambda_opt_1d,
  M = 2
)

ce_ds1 <- sqrt(mean(( (Y_bar + sY * ds1_c_star) - Yc )^2, na.rm = TRUE))
ce_ds2 <- sqrt(mean(( (Y_bar + sY * ds2_c_star) - Yc )^2, na.rm = TRUE))

ce_bs    <- sqrt(mean((Yc - bs_Yc)^2))
ce_ps    <- sqrt(mean((Yc - ps_Yc)^2))
ce_lo    <- sqrt(mean((Yc - LOESS_Yc)^2))
ce_ss    <- sqrt(mean((Yc - smooth_spline_Yc)^2))
ce_gpr   <- sqrt(mean((Yc - GPR_Yc)^2))
ce_atps  <- sqrt(mean((Yc - ATPS_Yc)^2))

sm_nw   <- smooth_index(y_hat_nw_1d)
sm_anw  <- smooth_index(y_hat_anw_1d)
sm_ds1  <- smooth_index(y_hat_ds1_1d)
sm_ds2  <- smooth_index(y_hat_ds2_1d)
sm_bs   <- smooth_index(bs_y)
sm_ps   <- smooth_index(ps_y)
sm_lo   <- smooth_index(LOESS_y)
sm_ss   <- smooth_index(smooth_spline)
sm_gpr  <- smooth_index(GPR_y)
sm_atps <- smooth_index(ATPS_y)

tau_c <- 0.1
tau_g <- 0.06
tau_s <- median(c(sm_nw, sm_anw, sm_ds1, sm_ds2), na.rm = TRUE)

w_c <- 0.5; w_g <- 0.3; w_s <- 0.2

css <- function(rmse_v, ce_v, sm_v) {
  w_c * phi(ce_v   / tau_c) +
    w_g * phi(rmse_v / tau_g) +
    w_s * phi(sm_v  / tau_s)
}

tab_all <- data.frame(
  Model         = c("NW", "ANW", "DS-ANW (M=1)", "DS-ANW (M=2)", "BS", "PS", "LOESS", "SS", "GPR", "ATPS"),
  RMSE          = c(rmse_nw, rmse_anw, rmse_ds1, rmse_ds2, rmse_bs, rmse_ps, rmse_lo, rmse_ss, rmse_gpr, rmse_atps),
  ConstraintErr = c(ce_nw,   ce_anw,   ce_ds1,   ce_ds2,   ce_bs,   ce_ps,   ce_lo,   ce_ss,   ce_gpr,   ce_atps),
  Smoothness    = c(sm_nw,   sm_anw,   sm_ds1,   sm_ds2,   sm_bs,   sm_ps,   sm_lo,   sm_ss,   sm_gpr,   sm_atps)
)

tab_all$CSS <- mapply(css, tab_all$RMSE, tab_all$ConstraintErr, tab_all$Smoothness)

print(tab_all)


cat("\nSaved figure:\n- figs/fig_railway_1d_allmethods_with_ds.pdf\n")
cat("If you uncomment write.csv, it will also save the summary table.\n")
