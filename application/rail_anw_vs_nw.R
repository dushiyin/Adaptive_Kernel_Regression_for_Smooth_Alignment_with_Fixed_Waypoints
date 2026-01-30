if (!require("dplyr"))   install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("maps"))    install.packages("maps")
if (!require("ggrepel")) install.packages("ggrepel")

library(dplyr)
library(ggplot2)
library(maps)
library(ggrepel)

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
stops_line    <- stops_line %>% arrange(s)

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
  lon=c(-123.0971,-120.3570,-118.0790,-113.4745,-106.6689,-97.1384,-91.9163,-80.8575,-79.3807),
  lat=c(49.2734,50.7168,52.8734,53.5426,52.1345,49.8951,50.0973,46.5382,43.6457),
  label_lon = c(-120.9, -122.7, -119.1, -113.0, -105.8, -97.9, -91.1, -84.0, -81.7),
  label_lat = c(48.2,   51.9,   54.2,   54.5,   53.2,   51.1,   51.2,   46.3,  43.7)
)

X_data <- c(stops_line$stop_lon, cons$lon)  
Y_data <- c(stops_line$stop_lat, cons$lat)  

Xc <- cons$lon   
Yc <- cons$lat 

Y_bar <- mean(Y_data)
s_Y   <- sqrt(mean((Y_data - Y_bar)^2))
if (!is.finite(s_Y) || s_Y <= 0) stop("s_Y is not positive; cannot standardize Y.")

Y_star  <- (Y_data - Y_bar) / s_Y
Yc_star <- (Yc - Y_bar) / s_Y

nw_1d_point <- function(x0, X, Ystar, h) {
  w   <- K1d((x0 - X) / h)
  den <- sum(w)
  if (den == 0) return(NA_real_)
  sum(w * Ystar) / den
}

anw_1d_point <- function(x0, X, Ystar, Xc, Yc_star, h, Lambda) {
  w_data <- K1d((x0 - X) / h)
  w_cons <- Lambda * K1d((x0 - Xc) / h)
  den <- sum(w_data) + sum(w_cons)
  if (den == 0) return(NA_real_)
  (sum(w_data * Ystar) + sum(w_cons * Yc_star)) / den
}

pred_on_grid_1d <- function(xg, FUN) vapply(xg, FUN, numeric(1))

anw_pred_vec_1d <- function(X_train, Ystar_train, X_eval, Xc, Yc_star, h, Lambda) {
  vapply(
    X_eval,
    function(x0) anw_1d_point(x0, X_train, Ystar_train, Xc, Yc_star, h, Lambda),
    numeric(1)
  )
}

ds_anw_pred_vec_1d <- function(X_train, Y0_star, X_eval, Xc, Yc_star, h, Lambda, M = 1) {
  if (M == 0) {
    return(anw_pred_vec_1d(X_train, Y0_star, X_eval, Xc, Yc_star, h, Lambda))
  }
  Ym <- Y0_star
  for (m in 1:M) {
    fit_on_train <- anw_pred_vec_1d(X_train, Ym, X_train, Xc, Yc_star, h, Lambda)
    Ym <- Y0_star + (Ym - fit_on_train) 
  }
  anw_pred_vec_1d(X_train, Ym, X_eval, Xc, Yc_star, h, Lambda)
}

make_folds <- function(n, K = 5, seed = 2025) {
  set.seed(seed)
  fold_id <- sample(rep(1:K, length.out = n))
  fold_id
}

cv_mse_nw_kfold <- function(h, X, Ystar, K = 5, seed = 2025) {
  if (h <= 0) return(Inf)
  n <- length(X)
  fold_id <- make_folds(n, K, seed)
  
  se_all <- c()
  for (k in 1:K) {
    idx_te <- which(fold_id == k)
    idx_tr <- setdiff(1:n, idx_te)
    
    Xtr <- X[idx_tr]; Ytr <- Ystar[idx_tr]
    Xte <- X[idx_te]; Yte <- Ystar[idx_te]
    
    yhat <- vapply(Xte, function(x0) nw_1d_point(x0, Xtr, Ytr, h), numeric(1))
    se   <- (Yte - yhat)^2
    se_all <- c(se_all, se)
  }
  mean(se_all, na.rm = TRUE)
}

cv_mse_anw_kfold <- function(Lambda, h, X, Ystar, Xc, Yc_star, K = 5, seed = 2025) {
  if (Lambda < 0 || h <= 0) return(Inf)
  n <- length(X)
  fold_id <- make_folds(n, K, seed)
  
  se_all <- c()
  for (k in 1:K) {
    idx_te <- which(fold_id == k)
    idx_tr <- setdiff(1:n, idx_te)
    
    Xtr <- X[idx_tr]; Ytr <- Ystar[idx_tr]
    Xte <- X[idx_te]; Yte <- Ystar[idx_te]
    
    yhat <- vapply(
      Xte,
      function(x0) anw_1d_point(x0, Xtr, Ytr, Xc, Yc_star, h, Lambda),
      numeric(1)
    )
    se <- (Yte - yhat)^2
    se_all <- c(se_all, se)
  }
  mean(se_all, na.rm = TRUE)
}

waypoint_sq_error_star <- function(Lambda, h, X, Ystar, Xc, Yc_star) {
  yhatC <- vapply(
    Xc,
    function(x0) anw_1d_point(x0, X, Ystar, Xc, Yc_star, h, Lambda),
    numeric(1)
  )
  sum((yhatC - Yc_star)^2)
}

Lx        <- diff(range(X_data))
h_grid_1d <- seq(0.02 * Lx, 0.3 * Lx, length.out = 25)

Kfold <- 5
seed_cv <- 2025

cv_nw_1d <- sapply(
  h_grid_1d,
  function(h) cv_mse_nw_kfold(h, X = X_data, Ystar = Y_star, K = Kfold, seed = seed_cv)
)
h_nw_opt_1d <- h_grid_1d[which.min(cv_nw_1d)]
cat(sprintf("1D NW optimal h: %.4f (longitude degrees)\n", h_nw_opt_1d))

xg <- seq(
  min(c(X_data, Xc)),
  max(c(X_data, Xc)),
  length.out = 1500
)

alpha <- 1  

loss_anw <- function(h, Lambda, X, Ystar, Xc, Yc_star, alpha, K = 5, seed = 2025) {
  if (h <= 0 || Lambda < 0) return(Inf)
  cv_star <- cv_mse_anw_kfold(Lambda, h, X, Ystar, Xc, Yc_star, K = K, seed = seed)
  wp_star <- waypoint_sq_error_star(Lambda, h, X, Ystar, Xc, Yc_star)
  cv_star + alpha * wp_star
}

Lambda_grid_1d <- seq(4, 50, by = 1)
h_grid_anw     <- h_grid_1d

loss_mat <- matrix(NA_real_, nrow = length(h_grid_anw), ncol = length(Lambda_grid_1d))

for (i in seq_along(h_grid_anw)) {
  for (j in seq_along(Lambda_grid_1d)) {
    loss_mat[i, j] <- loss_anw(
      h       = h_grid_anw[i],
      Lambda  = Lambda_grid_1d[j],
      X       = X_data,
      Ystar   = Y_star,
      Xc      = Xc,
      Yc_star = Yc_star,
      alpha   = alpha,
      K       = Kfold,
      seed    = seed_cv
    )
  }
}

idx_min <- which(loss_mat == min(loss_mat, na.rm = TRUE), arr.ind = TRUE)
h_anw_opt_1d  <- h_grid_anw[idx_min[1, 1]]
Lambda_opt_1d <- Lambda_grid_1d[idx_min[1, 2]]

cat(sprintf("1D ANW optimal h: %.4f, optimal Lambda: %.2f | alpha = %.1f | K = %d\n",
            h_anw_opt_1d, Lambda_opt_1d, alpha, Kfold))

y_hat_nw_star <- pred_on_grid_1d(xg, function(xx) {
  nw_1d_point(xx, X = X_data, Ystar = Y_star, h = h_nw_opt_1d)
})
y_hat_nw_1d <- Y_bar + s_Y * y_hat_nw_star

y_hat_anw_star <- pred_on_grid_1d(xg, function(xx) {
  anw_1d_point(
    xx,
    X       = X_data,
    Ystar   = Y_star,
    Xc      = Xc,
    Yc_star = Yc_star,
    h       = h_anw_opt_1d,
    Lambda  = Lambda_opt_1d
  )
})
y_hat_anw_1d <- Y_bar + s_Y * y_hat_anw_star

y_hat_ds1_star <- ds_anw_pred_vec_1d(
  X_train = X_data, Y0_star = Y_star, X_eval = xg,
  Xc = Xc, Yc_star = Yc_star,
  h = h_anw_opt_1d, Lambda = Lambda_opt_1d, M = 1
)
y_hat_ds1_1d <- Y_bar + s_Y * y_hat_ds1_star

y_hat_ds2_star <- ds_anw_pred_vec_1d(
  X_train = X_data, Y0_star = Y_star, X_eval = xg,
  Xc = Xc, Yc_star = Yc_star,
  h = h_anw_opt_1d, Lambda = Lambda_opt_1d, M = 2
)
y_hat_ds2_1d <- Y_bar + s_Y * y_hat_ds2_star

smooth_nw_1d <- data.frame(
  lon    = xg,
  lat    = y_hat_nw_1d,
  method = "NW"
)

smooth_anw_1d <- data.frame(
  lon    = xg,
  lat    = y_hat_anw_1d,
  method = "ANW"
)

smooth_ds1_1d <- data.frame(
  lon    = xg,
  lat    = y_hat_ds1_1d,
  method = "DS-ANW (M=1)"
)

smooth_ds2_1d <- data.frame(
  lon    = xg,
  lat    = y_hat_ds2_1d,
  method = "DS-ANW (M=2)"
)

baseline_1d <- data.frame(
  lon    = shape111$shape_pt_lon,
  lat    = shape111$shape_pt_lat,
  method = "Baseline"
)

curves_1d <- bind_rows(baseline_1d, smooth_nw_1d, smooth_anw_1d, smooth_ds1_1d, smooth_ds2_1d)

world <- map_data("world")

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
    title = "ANW vs NW (with DS-ANW)",
    subtitle = sprintf(
      "NW: h = %.4f | ANW: h = %.4f, Lambda = %.1f | DS-ANW uses same (h, Lambda) | alpha = %.1f",
      h_nw_opt_1d, h_anw_opt_1d, Lambda_opt_1d, alpha
    ),
    x = "Longitude", y = "Latitude",
    color = "Method", linetype = "Method", linewidth = "Method"
  ) +
  scale_color_manual(values = c(
    "Baseline"    = "black",
    "NW"          = "green",
    "ANW"         = "firebrick",
    "DS-ANW (M=1)"= "dodgerblue4",
    "DS-ANW (M=2)"= "purple4"
  )) +
  scale_linetype_manual(values = c(
    "Baseline"    = "solid",
    "NW"          = "dashed",
    "ANW"         = "dashed",
    "DS-ANW (M=1)"= "dotdash",
    "DS-ANW (M=2)"= "twodash"
  )) +
  scale_linewidth_manual(values = c(
    "Baseline"    = 0.8,
    "NW"          = 0.8,
    "ANW"         = 1.2,
    "DS-ANW (M=1)"= 1.0,
    "DS-ANW (M=2)"= 1.0
  )) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box      = "vertical"
  )

if (!dir.exists("figs")) dir.create("figs")

pdf(
  file   = "figs/fig_1d_railway_anw_vs_nw_ds.pdf",
  width  = 10,
  height = 4
)
print(p1d)
dev.off()

cat("Saved: figs/fig_1d_railway_anw_vs_nw_ds.pdf\n")
