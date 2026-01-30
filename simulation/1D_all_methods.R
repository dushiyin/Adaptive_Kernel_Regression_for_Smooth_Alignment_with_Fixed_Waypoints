set.seed(2025)
n       <- 500
noiseSD <- 0.15
x_min   <- 0
x_max   <- 10
h_fix   <- 0.25
ngrid   <- 400

if(!dir.exists("data")) dir.create("data")
if(!dir.exists("figs")) dir.create("figs")

m_true <- function(x) sin(x) + 0.3*cos(2*x)

x <- sort(runif(n, x_min, x_max))
y <- m_true(x) + rnorm(n, mean = 0, sd = noiseSD)

xg <- seq(x_min, x_max, length.out = ngrid)
truth_on_grid <- m_true(xg)

chosen_idx_fixed <- c(216, 304)

C_fixed <- data.frame(
  xC  = x[chosen_idx_fixed],
  yC  = y[chosen_idx_fixed],
  idx = chosen_idx_fixed
)

pick_two_close_to_baseline <- function(
    x, y,
    baseline = c("spline","llr","nw"),
    h = 0.25,
    top_p = 0.3,      
    min_dx = 0.0,    
    seed = 20251212
){
  baseline <- match.arg(baseline)
  set.seed(seed)
  
  if (baseline == "spline") {
    yhat <- as.numeric(predict(smooth.spline(x, y), x)$y)
  } else if (baseline == "llr") {
    K <- function(u) dnorm(u)
    llr_point <- function(x0, x, y, h) {
      w <- K((x - x0)/h)
      X <- cbind(1, x-x0)
      beta <- solve(t(X*w)%*%X + diag(1e-8,2), t(X*w)%*%y)
      as.numeric(beta[1])
    }
    yhat <- vapply(x, function(xx) llr_point(xx, x, y, h), numeric(1))
  } else { 
    K <- function(u) dnorm(u)
    nw_point <- function(x0, x, y, h) {
      w <- K((x0 - x) / h)
      if (sum(w)==0) return(mean(y))
      sum(w*y)/sum(w)
    }
    yhat <- vapply(x, function(xx) nw_point(xx, x, y, h), numeric(1))
  }
  
  res <- abs(y - yhat)
  k <- max(2, floor(length(x) * top_p))
  pool <- order(res)[1:k]   
  
  i1 <- sample(pool, 1)
  pool2 <- pool[abs(x[pool] - x[i1]) >= min_dx]
  if (length(pool2) == 0) pool2 <- setdiff(pool, i1)
  i2 <- sample(pool2, 1)
  
  sort(c(i1, i2))
}

chosen_idx_rand <- pick_two_close_to_baseline(
  x, y,
  baseline = "spline",  
  h = h_fix,
  top_p = 0.2,
  min_dx = 0.0,
  seed = 20251212
)

C_rand <- data.frame(
  xC  = x[chosen_idx_rand],
  yC  = y[chosen_idx_rand],
  idx = chosen_idx_rand
)

pick_constraints_internal_observed <- function(
    x, y, m, pilot = c("spline","llr"), h_llr = 0.25, min_dx = 0.4
){
  pilot <- match.arg(pilot)
  
  if (pilot == "spline") {
    y_hat <- as.numeric(predict(smooth.spline(x, y), x)$y)
  } else {
    K_loc <- function(u) dnorm(u)
    llr_point_loc <- function(x0, x, y, h) {
      w <- K_loc((x - x0)/h)
      X <- cbind(1, x - x0)
      beta <- solve(t(X * w) %*% X + diag(1e-8,2), t(X * w) %*% y)
      as.numeric(beta[1])
    }
    y_hat <- vapply(x, function(z) llr_point_loc(z, x, y, h_llr), numeric(1))
  }
  
  res <- abs(y - y_hat)
  res_cap <- quantile(res, 0.96)
  res_trim <- pmin(res, res_cap)
  ord <- order(res_trim, decreasing = TRUE)
  
  selected <- integer(0)
  for (i in ord) {
    if (length(selected) >= m) break
    if (length(selected) == 0 || all(abs(x[i] - x[selected]) >= min_dx))
      selected <- c(selected, i)
  }
  
  if (length(selected) < m) {
    need <- setdiff(ord, selected)
    selected <- c(selected, head(need, m - length(selected)))
  }
  
  data.frame(
    xC = x[selected],
    yC = y[selected],
    idx = selected,
    resid = res[selected]
  )
}

y_mean <- mean(y)
y_sd   <- sd(y)
y_std  <- (y - y_mean) / y_sd

K <- function(u) dnorm(u)

nw_point <- function(x0, x, y, h) {
  w <- K((x0 - x) / h)
  if (sum(w)==0) return(0)
  sum(w*y)/sum(w)
}

llr_point <- function(x0, x, y, h) {
  w <- K((x - x0)/h)
  X <- cbind(1, x-x0)
  beta <- solve(t(X*w)%*%X + diag(1e-8,2), t(X*w)%*%y)
  as.numeric(beta[1])
}

spline_predict_vec <- function(x_train, y_train, x_out) {
  fit <- smooth.spline(x_train, y_train)
  as.numeric(predict(fit, x_out)$y)
}
spline_predict_point <- function(x_train, y_train, x0) {
  as.numeric(predict(smooth.spline(x_train, y_train), x0)$y)
}

bspline_predict_vec <- function(X_data, Y_data, x_out) {
  library(splines)
  B <- bs(X_data)
  fit <- lm(Y_data ~ B)
  B_new <- bs(x_out)
  as.numeric(predict(fit, newdata=list(B=B_new)) )
}

bspline_predict_point <- function(X_data, Y_data, x0) {
  library(splines)
  B <- bs(X_data)
  fit <- lm(Y_data ~ B)
  as.numeric(predict(fit, newdata=list(B=bs(x0)) ))
}

pspline_predict_vec <- function(X_data, Y_data, x_out) {
  library(mgcv)
  fit=gam(Y_data ~ s(X_data,bs="ps",k=20, m=c(2,2)), method = "REML" )
  as.numeric(predict(fit, newdata=data.frame(X_data=x_out)))
}

pspline_predict_point <- function(X_data, Y_data, x0) {
  library(mgcv)
  fit=gam(Y_data ~ s(X_data,bs="ps",k=20, m=c(2,2)), method = "REML" )
  as.numeric(predict(fit, newdata=data.frame(X_data=x0)))
}

LOESS_predict_vec <- function(X_data, Y_data, x_out) {
  fit <- loess(Y_data ~ X_data,span=0.7,degree=2)
  as.numeric(predict(fit, data.frame(X_data=x_out)))
}

LOESS_predict_point <- function(X_data, Y_data, x0) {
  fit <- loess(Y_data ~ X_data,span=0.7,degree=2)
  as.numeric(predict(fit, data.frame(X_data=x0)))
}

GPR_predict_vec <-function(X_data, Y_data, x_out) {
  library(GauPro)
  fit <- gpkm(X_data, Y_data)
  as.numeric(fit$predict(x_out))
}

GPR_predict_point <-function(X_data, Y_data, x0) {
  library(GauPro)
  fit <- gpkm(X_data, Y_data)
  as.numeric(fit$predict(x0))
}

ATPS_predict_vec <-function(X_data, Y_data, x_out) {
  library(mgcv)
  fit <- gam(Y_data ~ s(X_data, bs="ts", k=20) + s(X_data, bs="ts", k=20, m=2) ,
             method="REML")
  as.numeric(predict(fit, data.frame(X_data=x_out)))
}

ATPS_predict_point <-function(X_data, Y_data, x0) {
  library(mgcv)
  fit <- gam(Y_data ~ s(X_data, bs="ts", k=20) + s(X_data, bs="ts", k=20, m=2) ,
             method="REML")
  as.numeric(predict(fit, data.frame(X_data=x0)))
}

anw_point <- function(x0, x, y, xC, yC, h, lambda) {
  w  <- K((x0 - x)/h)
  wC <- lambda * K((x0 - xC)/h)
  num <- sum(w*y) + sum(wC*yC)
  den <- sum(w) + sum(wC)
  if (den==0) return(0)
  num/den
}

pred_on_grid <- function(xg, fun) vapply(xg, fun, numeric(1))

anw_pred_vec <- function(x_train, y_train, x_eval, xC, yC, h, lambda) {
  vapply(x_eval, function(xx)
    anw_point(xx, x_train, y_train, xC, yC, h, lambda), numeric(1))
}

ds_anw_pred_vec <- function(x_train, y0, x_eval, xC, yC, h, lambda, M = 1) {
  
  if (M == 0) {
    return(anw_pred_vec(x_train, y0, x_eval, xC, yC, h, lambda))
  }
  ym <- y0
  for (m in 1:M) {
    fit_on_train <- anw_pred_vec(x_train, ym, x_train, xC, yC, h, lambda)
    ym <- y0 + (ym - fit_on_train) 
  }
  anw_pred_vec(x_train, ym, x_eval, xC, yC, h, lambda)
}

rmse <- function(truth, pred) sqrt(mean((truth-pred)^2))
max_abs_err <- function(truth,pred) max(abs(truth-pred))

constraint_rmse <- function(xC, yC, pred_fun_point) {
  yhatC <- vapply(xC, pred_fun_point, numeric(1))
  sqrt(mean((yC-yhatC)^2))
}

smooth_index <- function(yhat) {
  d2 <- diff(yhat, differences=2)
  mean(d2^2)
}

cv_mse_anw <- function(h, lambda, x, y_std, xC, yC_std, Kfold=5) {
  if (h<=0 || lambda<0) return(Inf)
  n <- length(x)
  folds <- cut(seq_len(n), breaks=Kfold, labels=FALSE)
  errs <- numeric(Kfold)
  for (k in 1:Kfold) {
    tr <- which(folds!=k); te <- which(folds==k)
    yhat <- vapply(x[te], function(xx)
      anw_point(xx, x[tr], y_std[tr], xC, yC_std, h, lambda), numeric(1))
    errs[k] <- mean((y_std[te]-yhat)^2)
  }
  mean(errs)
}

constraint_sq_error_std <- function(h, lambda, x, y_std, xC, yC_std) {
  yhatC_std <- vapply(xC, function(xx)
    anw_point(xx, x, y_std, xC, yC_std, h, lambda), numeric(1))
  mean((yhatC_std-yC_std)^2)
}

loss_anw <- function(h, lambda, x, y_std, xC, yC_std, alpha=1) {
  if (h<=0 || lambda<0) return(Inf)
  cv_mse_anw(h,lambda,x,y_std,xC,yC_std)+
    alpha*constraint_sq_error_std(h,lambda,x,y_std,xC,yC_std)
}

run_main_case <- function(case_name, C_df,
                          h_grid = seq(0.15,0.45,length.out=15),
                          lambda_grid = seq(0,100,1),
                          alpha = 1,
                          tau_c = 0.1, tau_g = 0.05,
                          w_c = 0.5, w_g = 0.3, w_s = 0.2,
                          fig_file = NULL,
                          csv_file = NULL){
  
  
  C_df$yC_std <- (C_df$yC - y_mean) / y_sd
  
  
  loss_mat <- matrix(NA_real_,
                     nrow=length(h_grid),
                     ncol=length(lambda_grid))
  
  for (i in seq_along(h_grid)) {
    for (j in seq_along(lambda_grid)) {
      loss_mat[i,j] <- loss_anw(h_grid[i], lambda_grid[j],
                                x, y_std, C_df$xC, C_df$yC_std, alpha=alpha)
    }
  }
  
  idx <- which(loss_mat==min(loss_mat,na.rm=TRUE), arr.ind=TRUE)
  h_opt  <- h_grid[idx[1,1]]
  L_opt  <- lambda_grid[idx[1,2]]
  
  cat(sprintf("[%s] ANW optimal h=%.4f, lambda=%.2f\n", case_name, h_opt, L_opt))
  
  
  pred_nw  <- pred_on_grid(xg, function(xx) nw_point(xx,x,y,h_fix))
  pred_llr <- pred_on_grid(xg, function(xx) llr_point(xx,x,y,h_fix))
  pred_sp  <- spline_predict_vec(x,y,xg)
  pred_bs  <- bspline_predict_vec(x,y,xg)
  pred_ps  <- pspline_predict_vec(x,y,xg)
  pred_LOESS <- LOESS_predict_vec(x,y,xg)
  pred_gpr <- GPR_predict_vec(x,y,xg)
  pred_atps <- ATPS_predict_vec(x,y,xg)
  
  pred_anw_std <- pred_on_grid(
    xg, function(xx) anw_point(xx,x,y_std,C_df$xC,C_df$yC_std,h_opt,L_opt)
  )
  pred_anw <- y_mean + y_sd*pred_anw_std
  
  
  pred_ds1_std <- ds_anw_pred_vec(
    x_train = x, y0 = y_std, x_eval = xg,
    xC = C_df$xC, yC = C_df$yC_std,
    h = h_opt, lambda = L_opt, M = 1
  )
  pred_ds1 <- y_mean + y_sd * pred_ds1_std
  
  
  pred_ds2_std <- ds_anw_pred_vec(
    x_train = x, y0 = y_std, x_eval = xg,
    xC = C_df$xC, yC = C_df$yC_std,
    h = h_opt, lambda = L_opt, M = 2
  )
  pred_ds2 <- y_mean + y_sd * pred_ds2_std
  
  
  tab_main <- data.frame(
    Case  = case_name,
    Model = c("NW","LLR","Spline","ANW","DS-ANW (M=1)","DS-ANW (M=2)","BS","PS","LOESS","GPR","ATPS"),
    RMSE  = c(rmse(truth_on_grid,pred_nw),
              rmse(truth_on_grid,pred_llr),
              rmse(truth_on_grid,pred_sp),
              rmse(truth_on_grid,pred_anw),
              rmse(truth_on_grid,pred_ds1),
              rmse(truth_on_grid,pred_ds2),
              rmse(truth_on_grid,pred_bs),
              rmse(truth_on_grid,pred_ps),
              rmse(truth_on_grid[-c(1,length(truth_on_grid))],pred_LOESS[-c(1,length(truth_on_grid))]),
              rmse(truth_on_grid,pred_gpr),
              rmse(truth_on_grid,pred_atps)
    ),
    Max   = c(max_abs_err(truth_on_grid,pred_nw),
              max_abs_err(truth_on_grid,pred_llr),
              max_abs_err(truth_on_grid,pred_sp),
              max_abs_err(truth_on_grid,pred_anw),
              max_abs_err(truth_on_grid,pred_ds1),
              max_abs_err(truth_on_grid,pred_ds2),
              max_abs_err(truth_on_grid,pred_bs),
              max_abs_err(truth_on_grid,pred_ps),
              max_abs_err(truth_on_grid[-c(1,length(truth_on_grid))],pred_LOESS[-c(1,length(truth_on_grid))]),
              max_abs_err(truth_on_grid,pred_gpr),
              max_abs_err(truth_on_grid,pred_atps)
    ),
    ConstraintErr = c(
      constraint_rmse(C_df$xC,C_df$yC, function(xx) nw_point(xx,x,y,h_fix)),
      constraint_rmse(C_df$xC,C_df$yC, function(xx) llr_point(xx,x,y,h_fix)),
      constraint_rmse(C_df$xC,C_df$yC, function(xx) spline_predict_point(x,y,xx)),
      constraint_rmse(C_df$xC,C_df$yC, function(xx){
        yhat_std <- anw_point(xx,x,y_std,C_df$xC,C_df$yC_std,h_opt,L_opt)
        y_mean + y_sd*yhat_std
      }),
      
      constraint_rmse(C_df$xC,C_df$yC, function(xx){
        yhat_std <- ds_anw_pred_vec(
          x_train = x, y0 = y_std, x_eval = xx,
          xC = C_df$xC, yC = C_df$yC_std,
          h = h_opt, lambda = L_opt, M = 1
        )
        y_mean + y_sd*yhat_std
      }),
      
      constraint_rmse(C_df$xC,C_df$yC, function(xx){
        yhat_std <- ds_anw_pred_vec(
          x_train = x, y0 = y_std, x_eval = xx,
          xC = C_df$xC, yC = C_df$yC_std,
          h = h_opt, lambda = L_opt, M = 2
        )
        y_mean + y_sd*yhat_std
      }),
      constraint_rmse(C_df$xC,C_df$yC, function(xx) bspline_predict_point(x,y,xx)),
      constraint_rmse(C_df$xC,C_df$yC, function(xx) pspline_predict_point(x,y,xx)),
      constraint_rmse(C_df$xC,C_df$yC, function(xx) LOESS_predict_point(x,y,xx)),
      constraint_rmse(C_df$xC,C_df$yC, function(xx) GPR_predict_point(x,y,xx)),
      constraint_rmse(C_df$xC,C_df$yC, function(xx) ATPS_predict_point(x,y,xx))
    ),
    Smoothness = c(smooth_index(pred_nw),
                   smooth_index(pred_llr),
                   smooth_index(pred_sp),
                   smooth_index(pred_anw),
                   smooth_index(pred_ds1),
                   smooth_index(pred_ds2),
                   smooth_index(pred_bs),
                   smooth_index(pred_ps),
                   smooth_index(pred_LOESS[-c(1,length(truth_on_grid))]),
                   smooth_index(pred_gpr),
                   smooth_index(pred_atps)
    )
  )
  
  
  phi <- function(z) pmax(0, z-1)
  tau_s <- median(tab_main$Smoothness)
  
  tab_main$CSS <- with(tab_main,
                       w_c * phi(ConstraintErr/tau_c) +
                         w_g * phi(RMSE/tau_g) +
                         w_s * phi(Smoothness/tau_s))
  
  
  if (!is.null(csv_file)) {
    write.csv(tab_main, csv_file, row.names=FALSE)
  }
  

  if (!is.null(fig_file)) {
    pdf(fig_file, width=8, height=6)
    par(mar=c(4,4,2,1))
    plot(xg,truth_on_grid,type="l",col="black",lwd=3,
         xlab="x",ylab="value",
         main=paste0("1D Spatial Reconstruction (", case_name, ")"))
    lines(xg,pred_nw,col="dodgerblue",lwd=2)
    
    
    lines(xg,pred_anw,col="firebrick",lwd=2,lty=4)
    
    
    lines(xg,pred_ds1,col="darkgreen",lwd=2,lty=10)
    lines(xg,pred_ds2,col="darkgreen",lwd=2,lty=11)
    
    lines(xg,pred_bs,col="purple",lwd=2,lty=5)
    lines(xg,pred_ps,col="orange",lwd=2,lty=6)
    lines(xg,pred_LOESS,col="blue",lwd=2,lty=7)
    lines(xg,pred_gpr,col="gold",lwd=2,lty=8)
    lines(xg,pred_atps,col="skyblue",lwd=2,lty=9)
    
    points(x,y,pch=16,col=rgb(0,0,0,0.25),cex=0.6)
    points(C_df$xC,C_df$yC,pch=16,col="red",cex=1.5)
    
    
    anw_lab  <- bquote(ANW~"(h="*.(sprintf("%.3f", h_opt))*", "*lambda*"="*.(sprintf("%.1f", L_opt))*")")
    ds1_lab  <- bquote(DS-ANW~"(M=1)")
    ds2_lab  <- bquote(DS-ANW~"(M=2)")
    
    legend("bottomright",
           legend=c("Truth","NW", anw_lab, ds1_lab, ds2_lab,
                    "BS", "PS","LOESS", "GPR", "ATPS","Control points","Samples"),
           col=c("black","dodgerblue","firebrick","darkgreen","darkgreen",
                 "purple","orange","blue","gold","skyblue","red",rgb(0,0,0,0.25)),
           lty=c(1,1,4,10,11,5,6,7,8,9,NA,NA),
           pch=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,16,16),
           lwd=c(3,2,2,2,2,2,2,2,2,2,NA,NA),
           bty="n")
    
    grid()
    dev.off()
  }
  
  list(
    tab = tab_main,
    h_opt = h_opt,
    lambda_opt = L_opt
  )
}

res_fixed <- run_main_case(
  case_name = "Fixed control points",
  C_df      = C_fixed,
  fig_file  = "figs/fig_1d_main_fixed_com.pdf",
  csv_file  = "data/exp1d_summary_fixed_com.csv"
)

res_rand <- run_main_case(
  case_name = "Random control points",
  C_df      = C_rand,
  fig_file  = "figs/fig_1d_main_random_com.pdf",
  csv_file  = "data/exp1d_summary_random_com.csv"
)

tab_two <- rbind(res_fixed$tab, res_rand$tab)
write.csv(tab_two, "data/exp1d_summary_two_cases_com.csv", row.names=FALSE)

cat("\n=== Fixed + Random===\n")
print(tab_two)

m_grid      <- c(1,2,3)
lambda_grid <- c(10,100,1000)

sens_rows <- list()

for (m_now in m_grid) {
  C_now <- pick_constraints_internal_observed(
    x, y, m=m_now, pilot="spline", h_llr=h_fix, min_dx=0.4
  )
  for (L in lambda_grid) {
    pred_now <- pred_on_grid(
      xg, function(xx) anw_point(xx,x,y,C_now$xC,C_now$yC,h_fix,L)
    )
    sens_rows[[length(sens_rows)+1]] <- data.frame(
      m=m_now, lambda=L,
      RMSE = rmse(truth_on_grid,pred_now),
      ConstraintErr = constraint_rmse(C_now$xC,C_now$yC,
                                      function(xx) anw_point(xx,x,y,C_now$xC,C_now$yC,h_fix,L))
    )
  }
}

tab_sens <- do.call(rbind, sens_rows)
write.csv(tab_sens, "data/exp1d_sensitivity_com.csv", row.names=FALSE)

cat("\n=== sensitivity completed ===\n")
print(tab_sens)

cat("\n print：\n",
    "- data/exp1d_summary_fixed_com.csv\n",
    "- data/exp1d_summary_random_com.csv\n",
    "- data/exp1d_summary_two_cases_com.csv\n",
    "- data/exp1d_sensitivity_com.csv\n",
    "- figs/fig_1d_main_fixed_com.pdf\n",
    "- figs/fig_1d_main_random_com.pdf\n", sep="")
