### Load packages
library(dplyr) 
library(tidyverse)
library(lubridate)
library(stringr)
library(zoo)
library(ggplot2)
library(urbnmapr)
library(devtools)
library(readxl)
library(spdep)
library(sp)
library(huge)
library(INLA)
library(HMMpa)
library(invgamma)
library(brinla)
library(reshape2)
library(patchwork)
library(jsonlite)
library(geosphere)
library(urbnmapr)
library(RAQSAPI)
library(con2aqi)
library(pscl)
library(randtoolbox)
library(lhs)
library(scales)

### Load INLA workspace 
load("~/work/P2-Causal-Analysis-NB-v2-6.18.RData")

set.seed(10)

### Treatment/PS models

kgr_model1_treatment = function(dataset,rho_rbf,rho_periodic,sigma2,link=1){
  
  ###Computer covariance
  K_df = dataset %>% select(time,cum_monthly_total_pm25_6) %>% st_drop_geometry()
  K_cov = cov_similarity_kernel(K_df,time_span = 72,rho_rbf = rho_rbf,rho_periodic = rho_periodic,sigma2 = sigma2)
  
  #Calculate trace norm of gram matrix
  K_cov_weight = norm((1/72)*K_cov,type = "F")
  
  #Calculate proposed kernel
  covGP1 = kronecker((K_cov/72),(H^2/58))
  
  #Need to ensure precision matrix is not computationally singular i.e det > 0
  covGP_jittered = desingularize(covGP1,threshold = 1e-5,increment = 0.2)
  covGP1 = covGP_jittered[[1]]
  num_jitters = covGP_jittered[[2]]
  
  inv_covGP1 = solve(covGP1)
  
  ###Fit INLA model 
  kgr_formula1 = (cum_monthly_total_pm25_6 + 0.1) ~ -1 + Intercept_1 + Intercept_2 + Intercept_3 + Intercept_4 + Intercept_5 + Intercept_6 + Intercept_7 + Intercept_8 + Intercept_9 + Intercept_10 + Intercept_11 + Intercept_12 + Intercept_13 + Intercept_14 + Intercept_15 + Intercept_16 + Intercept_17 + Intercept_18 + Intercept_19 + Intercept_20 + Intercept_21 + Intercept_22 + Intercept_23 + Intercept_24 + Intercept_25 + Intercept_26 + Intercept_27 + Intercept_28 + Intercept_29 + Intercept_30 + Intercept_31 + Intercept_32 + Intercept_33 + Intercept_34 + Intercept_35 + Intercept_36 + Intercept_37 + Intercept_38 + Intercept_39 + Intercept_40 + Intercept_41 + Intercept_42 + Intercept_43 + Intercept_44 + Intercept_45 + Intercept_46 + Intercept_47 + Intercept_48 + Intercept_49 + Intercept_50 + Intercept_51 + Intercept_52 + Intercept_53 + Intercept_54 + Intercept_55 + Intercept_56 + Intercept_57 + Intercept_58 + month + aqi_resids + total_precip + avg_air_temp + avg_dew_point + avg_wind_speed + cum_monthly_total_pm25_6_lag1 + cum_monthly_total_pm25_6_lag2 + cum_monthly_total_pm25_6_lag3 + cum_monthly_total_pm25_6_lag4 + cum_monthly_total_pm25_6_lag5 + cum_monthly_total_pm25_6_lag6 + f(id2, model = "generic0", Cmatrix = inv_covGP1)  
  model = inla(formula = kgr_formula1,family = "gamma",data = dataset, num.threads = 50,
               control.compute = list(dic=TRUE,waic=TRUE,
                                      return.marginals.predictor=TRUE),
               control.inla = list(strategy = "simplified.laplace"),
               control.predictor = list(compute = TRUE, link = link))
  
  ###Extract relevant information and store in the list
  # model_summary <- model$summary.fixed
  # bri_hyperpar_summary <- bri.hyperpar.summary(model)
  model_DIC <- model$dic$dic
  model_WAIC <- model$waic$waic
  # preds_model <- model$summary.fitted.values
  # preds_model <- cbind(dataset$NAME, dataset$time, dataset$ID, dataset$monthly_total_pm25, preds_model)
  # colnames(preds_model) <- c("county", "time", "ID", "monthly_total_pm25", "mean", "sd", 
  #                            "0.025quant", "0.5quant", "0.975quant", "mode")
  
  # #Exponentiating parameter to get better interpretation of estimates 
  # multeff <- exp(model$summary.fixed$mean)
  # names(multeff) <- model$names.fixed
  # 
  # #Plot of each parameters' posterior density 
  # mf <- melt(model$marginals.fixed)
  # cf <- spread(mf,Var2,value)
  # names(cf)[2] <- 'parameter'
  # param_plot = ggplot(cf,aes(x=x,y=y)) + geom_line()+facet_wrap(~ parameter, 
  #                                                               scales="free") + geom_vline(xintercept=0) + ylab("density")
  # 
  # #Plot of precision of random effect (main hyperparameter of interest)
  # sden <- data.frame(bri.hyper.sd(model$marginals.hyperpar[[1]]))
  # hyperparam_plot = ggplot(sden,aes(x,y)) + geom_line() + ylab("density") + 
  #   xlab("linear predictor")
  
  #Store the results in the list
  kgr_model1_results = list(
    #covmatrix = covGP1,
    #prec = inv_covGP1,
    #model_summary = model_summary,
    #bri_hyperpar_summary = bri_hyperpar_summary,
    #exp_effects = multeff,
    #param_plot = param_plot,
    #hyperparam_plot = hyperparam_plot,
    model_DIC = model_DIC,
    model_WAIC = model_WAIC
    #fitted_values = preds_model
  )
  
  return(kgr_model1_results)
}

### Define grid for 3 parameters in K_time
hyper_grid = randomLHS(100,3)
colnames(hyper_grid) = c("rho_rbf","rho_periodic","sigma2")
hyper_grid[,1:2] = hyper_grid[,1:2]*100
hyper_grid[,3] = hyper_grid[,3]*100

head(hyper_grid,20)

### Perform the grid search
apply_function <- function(row, index) {
  rho_rbf <- row[1]
  rho_periodic <- row[2]
  sigma2 <- row[3]
  
  result <- tryCatch(
    {
      kgr_model1_treatment(dataset = inla_df, rho_rbf = rho_rbf,
                           rho_periodic = rho_periodic, sigma2 = sigma2, link=1)
    },
    error = function(e) {
      # Handle the error by printing a message and returning NULL
      cat("Error in model", index, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  # If result is not NULL, process the result
  if (!is.null(result)) {
    
    kgr_model1_treatment_running_results <<- c(kgr_model1_treatment_running_results, result$model_WAIC)
    # Print and save WAIC every 10th model
    if (index %% 10 == 0) {
      cat("Model", index, "WAIC:", result$model_WAIC, "\n")
      saveRDS(kgr_model1_treatment_running_results, paste0("kgr_model1_treatment_running_results_WAIC.rds"))
    }
  }
  
  # Return the result (which could be NULL if there was an error)
  return(result)
}

kgr_model1_treatment_running_results <- vector()  # Initialize an empty vector to store WAIC results

kgr_model1_treatment_results_list <- lapply(seq_len(nrow(hyper_grid)), function(i) {
  apply_function(hyper_grid[i, ], i)
})

# Remove NULL results from the list if needed
kgr_model1_treatment_results_list <- Filter(Negate(is.null), kgr_model1_treatment_results_list)

#Extracting WAIC values 
kgr_model1_treatment_results_WAIC = c()

for (i in 1:length(kgr_model1_treatment_results_list)){
  kgr_model1_treatment_results_WAIC[i] = pred_data = kgr_model1_treatment_results_list[[i]]$model_WAIC
}


hist(kgr_model1_treatment_results_WAIC)
top5 = head(sort(kgr_model1_treatment_results_WAIC))
top5
top5_idx = which(kgr_model1_treatment_results_WAIC <= top5[5])

kgr_model1_treatment_results_WAIC = cbind(hyper_grid,kgr_model1_treatment_results_WAIC)
colnames(kgr_model1_treatment_results_WAIC) = c(colnames(hyper_grid),"WAIC")
kgr_model1_treatment_results_WAIC[top5_idx,]


saveRDS(kgr_model1_treatment_results_WAIC, paste0("kgr_model1_treatment_results_WAIC.rds"))


#-------------------------------------------------------------------------------------------------------


kgr_model2_treatment = function(dataset,rho_rbf,rho_periodic,sigma2,link=1){
  
  #Calculating gram matrix K_time
  K_time = time_kernel(time_span = length(unique(inla_df$time)),rho_rbf = rho_rbf, 
                       rho_periodic = rho_periodic, sigma2 = sigma2)
  
  #Calculate proposed kernel
  covGP2 = kronecker((K_time/72),(H^2/58))
  
  #Need to ensure precision matrix is not computationally singular i.e det > 0
  covGP_jittered = desingularize(covGP2,threshold = 1e-5,increment = 0.5)
  covGP2 = covGP_jittered[[1]]
  num_jitters = covGP_jittered[[2]]
  
  inv_covGP2 = solve(covGP2)
  
  ###Fit INLA model 
  kgr_formula2 = (cum_monthly_total_pm25_6 + 0.1) ~ -1 + Intercept_1 + Intercept_2 + Intercept_3 + Intercept_4 + Intercept_5 + Intercept_6 + Intercept_7 + Intercept_8 + Intercept_9 + Intercept_10 + Intercept_11 + Intercept_12 + Intercept_13 + Intercept_14 + Intercept_15 + Intercept_16 + Intercept_17 + Intercept_18 + Intercept_19 + Intercept_20 + Intercept_21 + Intercept_22 + Intercept_23 + Intercept_24 + Intercept_25 + Intercept_26 + Intercept_27 + Intercept_28 + Intercept_29 + Intercept_30 + Intercept_31 + Intercept_32 + Intercept_33 + Intercept_34 + Intercept_35 + Intercept_36 + Intercept_37 + Intercept_38 + Intercept_39 + Intercept_40 + Intercept_41 + Intercept_42 + Intercept_43 + Intercept_44 + Intercept_45 + Intercept_46 + Intercept_47 + Intercept_48 + Intercept_49 + Intercept_50 + Intercept_51 + Intercept_52 + Intercept_53 + Intercept_54 + Intercept_55 + Intercept_56 + Intercept_57 + Intercept_58 + month + aqi_resids + total_precip + avg_air_temp + avg_dew_point + avg_wind_speed + cum_monthly_total_pm25_6_lag1 + cum_monthly_total_pm25_6_lag2 + cum_monthly_total_pm25_6_lag3 + cum_monthly_total_pm25_6_lag4 + cum_monthly_total_pm25_6_lag5 + cum_monthly_total_pm25_6_lag6 + f(id2, model = "generic0", Cmatrix = inv_covGP2)  
  model = inla(formula = kgr_formula2,family = "gamma",data = dataset, num.threads = 50,
               control.compute = list(dic=TRUE,waic=TRUE,
                                      return.marginals.predictor=TRUE),
               control.inla = list(strategy = "simplified.laplace"),
               control.predictor = list(compute = TRUE, link = link))
  
  ###Extract relevant information and store in the list
  # model_summary <- model$summary.fixed
  # bri_hyperpar_summary <- bri.hyperpar.summary(model)
  model_DIC <- model$dic$dic
  model_WAIC <- model$waic$waic
  # preds_model <- model$summary.fitted.values
  # preds_model <- cbind(dataset$NAME, dataset$time, dataset$ID, dataset$monthly_total_pm25, preds_model)
  # colnames(preds_model) <- c("county", "time", "ID", "monthly_total_pm25", "mean", "sd", 
  #                            "0.025quant", "0.5quant", "0.975quant", "mode")
  
  # #Exponentiating parameter to get better interpretation of estimates 
  # multeff <- exp(model$summary.fixed$mean)
  # names(multeff) <- model$names.fixed
  # 
  # #Plot of each parameters' posterior density 
  # mf <- melt(model$marginals.fixed)
  # cf <- spread(mf,Var2,value)
  # names(cf)[2] <- 'parameter'
  # param_plot = ggplot(cf,aes(x=x,y=y)) + geom_line()+facet_wrap(~ parameter, 
  #                                                               scales="free") + geom_vline(xintercept=0) + ylab("density")
  # 
  # #Plot of precision of random effect (main hyperparameter of interest)
  # sden <- data.frame(bri.hyper.sd(model$marginals.hyperpar[[1]]))
  # hyperparam_plot = ggplot(sden,aes(x,y)) + geom_line() + ylab("density") + 
  #   xlab("linear predictor")
  
  #Store the results in the list
  kgr_model2_results = list(
    #covmatrix = covGP1,
    #prec = inv_covGP1,
    #model_summary = model_summary,
    #bri_hyperpar_summary = bri_hyperpar_summary,
    #exp_effects = multeff,
    #param_plot = param_plot,
    #hyperparam_plot = hyperparam_plot,
    model_DIC = model_DIC,
    model_WAIC = model_WAIC
    #fitted_values = preds_model
  )
  
  return(kgr_model2_results)
}

### Define grid for 3 parameters in K_time
hyper_grid = randomLHS(100,3)
colnames(hyper_grid) = c("rho_rbf","rho_periodic","sigma2")
hyper_grid[,1:2] = hyper_grid[,1:2]*100
hyper_grid[,3] = hyper_grid[,3]*100

head(hyper_grid,20)

### Perform the grid search
apply_function <- function(row, index) {
  rho_rbf <- row[1]
  rho_periodic <- row[2]
  sigma2 <- row[3]
  
  result <- tryCatch(
    {
      kgr_model2_treatment(dataset = inla_df, rho_rbf = rho_rbf,
                           rho_periodic = rho_periodic, sigma2 = sigma2, link=1)
    },
    error = function(e) {
      # Handle the error by printing a message and returning NULL
      cat("Error in model", index, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  # If result is not NULL, process the result
  if (!is.null(result)) {
    
    kgr_model2_treatment_running_results <<- c(kgr_model2_treatment_running_results, result$model_WAIC)
    # Print and save WAIC every 10th model
    if (index %% 10 == 0) {
      cat("Model", index, "WAIC:", result$model_WAIC, "\n")
      saveRDS(kgr_model2_treatment_running_results, paste0("kgr_model2_treatment_running_results_WAIC.rds"))
    }
  }
  
  # Return the result (which could be NULL if there was an error)
  return(result)
}

kgr_model2_treatment_running_results <- vector()  # Initialize an empty vector to store WAIC results

kgr_model2_treatment_results_list <- lapply(seq_len(nrow(hyper_grid)), function(i) {
  apply_function(hyper_grid[i, ], i)
})

# Remove NULL results from the list if needed
kgr_model2_treatment_results_list <- Filter(Negate(is.null), kgr_model2_treatment_results_list)

#Extracting WAIC values 
kgr_model2_treatment_results_WAIC = c()

for (i in 1:length(kgr_model2_treatment_results_list)){
  kgr_model2_treatment_results_WAIC[i] = pred_data = kgr_model2_treatment_results_list[[i]]$model_WAIC
}


hist(kgr_model2_treatment_results_WAIC)
top5 = head(sort(kgr_model2_treatment_results_WAIC))
top5
top5_idx = which(kgr_model2_treatment_results_WAIC <= top5[5])

kgr_model2_treatment_results_WAIC = cbind(hyper_grid,kgr_model2_treatment_results_WAIC)
colnames(kgr_model2_treatment_results_WAIC) = c(colnames(hyper_grid),"WAIC")
kgr_model2_treatment_results_WAIC[top5_idx,]


saveRDS(kgr_model2_treatment_results_WAIC, paste0("kgr_model2_treatment_results_WAIC.rds"))


#-------------------------------------------------------------------------------------------------------


### Outcome model (NB)
kgr_model1_response = function(dataset,rho_rbf,rho_periodic,sigma2,link=1){
  
  ###Computer covariance
  K_df = dataset %>% select(time,cum_monthly_total_pm25_6) %>% st_drop_geometry()
  K_cov = cov_similarity_kernel(K_df,time_span = 72,rho_rbf = rho_rbf,rho_periodic = rho_periodic,sigma2 = sigma2)
  
  #Calculate trace norm of gram matrix
  K_cov_weight = norm((1/72)*K_cov,type = "F")
  
  #Calculate proposed kernel
  covGP1 = kronecker((K_cov/72),(H^2/58))
  
  #Need to ensure precision matrix is not computationally singular i.e det > 0
  covGP_jittered = desingularize(covGP1,threshold = 1e-5,increment = 0.2)
  covGP1 = covGP_jittered[[1]]
  num_jitters = covGP_jittered[[2]]
  
  inv_covGP1 = solve(covGP1)
  
  # ###Fit INLA model 
  kgr_formula1 = deaths ~ -1 + offset(log(Total_Pop)) + factor(month) + aqi_resids + cum_total_precip_6 + avg_avg_air_temp_6 + avg_avg_dew_point_6 + max_avg_wind_speed_6 + Cluster*ns(cum_monthly_total_pm25_6,df = 3) + stabilized_IPW_kgr + f(id2,model = "generic0",Cmatrix = inv_covGP1)
  
  
  model = inla(formula = kgr_formula1, family = "nbinomial", data = inla_cf_df, num.threads = 50,
               control.compute = list(dic=TRUE,waic=TRUE, return.marginals.predictor=TRUE),
               #control.family = list(variant = 1), E = Total_Pop,
               control.family = list(hyper = list(theta = list(prior = "pc.mgamma", param = 1))),
               control.inla = list(strategy = "simplified.laplace"),
               control.predictor = list(compute = TRUE, link = 1))
  
  ###Extract relevant information and store in the list
  # model_summary <- model$summary.fixed
  # bri_hyperpar_summary <- bri.hyperpar.summary(model)
  model_DIC <- model$dic$dic
  model_WAIC <- model$waic$waic
  # preds_model <- model$summary.fitted.values
  # preds_model <- cbind(dataset$NAME, dataset$time, dataset$ID, dataset$deaths, preds_model)
  # colnames(preds_model) <- c("county", "time", "ID", "deaths", "mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
  
  # #Exponentiating parameter to get better interpretation of estimates 
  # multeff <- exp(model$summary.fixed$mean)
  # names(multeff) <- model$names.fixed
  # 
  # #Plot of each parameters' posterior density 
  # mf <- melt(model$marginals.fixed)
  # cf <- spread(mf,Var2,value)
  # names(cf)[2] <- 'parameter'
  # param_plot = ggplot(cf,aes(x=x,y=y)) + geom_line()+facet_wrap(~ parameter, 
  #                                                               scales="free") + geom_vline(xintercept=0) + ylab("density")
  # 
  # #Plot of precision of random effect (main hyperparameter of interest)
  # sden <- data.frame(bri.hyper.sd(model$marginals.hyperpar[[1]]))
  # hyperparam_plot = ggplot(sden,aes(x,y)) + geom_line() + ylab("density") + 
  #   xlab("linear predictor")
  
  #Store the results in the list
  kgr_model1_results = list(
    #covmatrix = covGP1,
    #prec = inv_covGP1,
    #model_summary = model_summary,
    #bri_hyperpar_summary = bri_hyperpar_summary,
    #exp_effects = multeff,
    #param_plot = param_plot,
    #hyperparam_plot = hyperparam_plot,
    model_DIC = model_DIC,
    model_WAIC = model_WAIC
    #fitted_values = preds_model
  )
  
  return(kgr_model1_results)
}

### Perform the grid search
apply_function <- function(row, index) {
  rho_rbf <- row[1]
  rho_periodic <- row[2]
  sigma2 <- row[3]
  
  result <- tryCatch(
    {
      kgr_model1_response(dataset = inla_df, rho_rbf = rho_rbf,
                          rho_periodic = rho_periodic, sigma2= sigma2, link=1)
    },
    error = function(e) {
      # Handle the error by printing a message and returning NULL
      cat("Error in model", index, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  # If result is not NULL, process the result
  if (!is.null(result)) {
    
    kgr_model1_response_running_results <<- c(kgr_model1_response_running_results, result$model_WAIC)
    # Print and save WAIC every 10th model
    if (index %% 10 == 0) {
      cat("Model", index, "WAIC:", result$model_WAIC, "\n")
      saveRDS(kgr_model1_response_running_results, paste0("kgr_model1_response_running_results_WAIC.rds"))
    }
  }
  
  # Return the result (which could be NULL if there was an error)
  return(result)
}

kgr_model1_response_running_results <- vector()  # Initialize an empty vector to store WAIC results

kgr_model1_response_results_list <- lapply(seq_len(nrow(hyper_grid)), function(i) {
  apply_function(hyper_grid[i, ], i)
})

# Remove NULL results from the list if needed
kgr_model1_response_results_list <- Filter(Negate(is.null), kgr_model1_response_results_list)

#Extracting WAIC values 
kgr_model1_response_results_WAIC = c()

for (i in 1:length(kgr_model1_response_results_list)){
  kgr_model1_response_results_WAIC[i] = pred_data = kgr_model1_response_results_list[[i]]$model_WAIC
}



hist(kgr_model1_response_results_WAIC)
top5 = head(sort(kgr_model1_response_results_WAIC))
top5
top5_idx = which(kgr_model1_response_results_WAIC <= top5[5])

kgr_model1_response_results_WAIC = cbind(hyper_grid,kgr_model1_response_results_WAIC)
colnames(kgr_model1_response_results_WAIC) = c(colnames(hyper_grid),"WAIC")
kgr_model1_response_results_WAIC[top5_idx,]


saveRDS(kgr_model1_response_results_WAIC, paste0("kgr_model1_response_results_WAIC.rds"))


#-------------------------------------------------------------------------------------------------------


### Outcome model (NB)
kgr_model2_response = function(dataset,rho_rbf,rho_periodic,sigma2,link=1){
  
  #Calculating gram matrix K_time
  K_time = time_kernel(time_span = length(unique(inla_df$time)),rho_rbf = rho_rbf, 
                       rho_periodic = rho_periodic, sigma2 = sigma2)
  
  #Calculate proposed kernel
  covGP2 = kronecker((K_time/72),(H^2/58))
  
  #Need to ensure precision matrix is not computationally singular i.e det > 0
  covGP_jittered = desingularize(covGP2,threshold = 1e-5,increment = 0.5)
  covGP2 = covGP_jittered[[1]]
  num_jitters = covGP_jittered[[2]]
  
  inv_covGP2 = solve(covGP2)
  
  # ###Fit INLA model 
  kgr_formula2 = deaths ~ -1 + offset(log(Total_Pop)) + factor(month) + aqi_resids + cum_total_precip_6 + avg_avg_air_temp_6 + avg_avg_dew_point_6 + max_avg_wind_speed_6 + Cluster*ns(cum_monthly_total_pm25_6,df = 3) + stabilized_IPW_kgr + f(id2,model = "generic0",Cmatrix = inv_covGP2)
  
  
  model = inla(formula = kgr_formula2, family = "nbinomial", data = inla_cf_df, num.threads = 50,
               control.compute = list(dic=TRUE,waic=TRUE, return.marginals.predictor=TRUE),
               #control.family = list(variant = 1), E = Total_Pop,
               control.family = list(hyper = list(theta = list(prior = "pc.mgamma", param = 1))),
               control.inla = list(strategy = "simplified.laplace"),
               control.predictor = list(compute = TRUE, link = 1))
  
  ###Extract relevant information and store in the list
  # model_summary <- model$summary.fixed
  # bri_hyperpar_summary <- bri.hyperpar.summary(model)
  model_DIC <- model$dic$dic
  model_WAIC <- model$waic$waic
  # preds_model <- model$summary.fitted.values
  # preds_model <- cbind(dataset$NAME, dataset$time, dataset$ID, dataset$deaths, preds_model)
  # colnames(preds_model) <- c("county", "time", "ID", "deaths", "mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
  
  # #Exponentiating parameter to get better interpretation of estimates 
  # multeff <- exp(model$summary.fixed$mean)
  # names(multeff) <- model$names.fixed
  # 
  # #Plot of each parameters' posterior density 
  # mf <- melt(model$marginals.fixed)
  # cf <- spread(mf,Var2,value)
  # names(cf)[2] <- 'parameter'
  # param_plot = ggplot(cf,aes(x=x,y=y)) + geom_line()+facet_wrap(~ parameter, 
  #                                                               scales="free") + geom_vline(xintercept=0) + ylab("density")
  # 
  # #Plot of precision of random effect (main hyperparameter of interest)
  # sden <- data.frame(bri.hyper.sd(model$marginals.hyperpar[[1]]))
  # hyperparam_plot = ggplot(sden,aes(x,y)) + geom_line() + ylab("density") + 
  #   xlab("linear predictor")
  
  #Store the results in the list
  kgr_model2_results = list(
    #covmatrix = covGP1,
    #prec = inv_covGP1,
    #model_summary = model_summary,
    #bri_hyperpar_summary = bri_hyperpar_summary,
    #exp_effects = multeff,
    #param_plot = param_plot,
    #hyperparam_plot = hyperparam_plot,
    model_DIC = model_DIC,
    model_WAIC = model_WAIC
    #fitted_values = preds_model
  )
  
  return(kgr_model2_results)
}

### Perform the grid search
apply_function <- function(row, index) {
  rho_rbf <- row[1]
  rho_periodic <- row[2]
  sigma2 <- row[3]
  
  result <- tryCatch(
    {
      kgr_model2_response(dataset = inla_df, rho_rbf = rho_rbf,
                          rho_periodic = rho_periodic, sigma2= sigma2, link=1)
    },
    error = function(e) {
      # Handle the error by printing a message and returning NULL
      cat("Error in model", index, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  # If result is not NULL, process the result
  if (!is.null(result)) {
    
    kgr_model2_response_running_results <<- c(kgr_model2_response_running_results, result$model_WAIC)
    # Print and save WAIC every 10th model
    if (index %% 10 == 0) {
      cat("Model", index, "WAIC:", result$model_WAIC, "\n")
      saveRDS(kgr_model2_response_running_results, paste0("kgr_model2_response_running_results_WAIC.rds"))
    }
  }
  
  # Return the result (which could be NULL if there was an error)
  return(result)
}

kgr_model2_response_running_results <- vector()  # Initialize an empty vector to store WAIC results

kgr_model2_response_results_list <- lapply(seq_len(nrow(hyper_grid)), function(i) {
  apply_function(hyper_grid[i, ], i)
})

# Remove NULL results from the list if needed
kgr_model2_response_results_list <- Filter(Negate(is.null), kgr_model2_response_results_list)

#Extracting WAIC values 
kgr_model2_response_results_WAIC = c()

for (i in 1:length(kgr_model2_response_results_list)){
  kgr_model2_response_results_WAIC[i] = pred_data = kgr_model2_response_results_list[[i]]$model_WAIC
}



hist(kgr_model2_response_results_WAIC)
top5 = head(sort(kgr_model2_response_results_WAIC))
top5
top5_idx = which(kgr_model2_response_results_WAIC <= top5[5])

kgr_model2_response_results_WAIC = cbind(hyper_grid,kgr_model2_response_results_WAIC)
colnames(kgr_model2_response_results_WAIC) = c(colnames(hyper_grid),"WAIC")
kgr_model2_response_results_WAIC[top5_idx,]


saveRDS(kgr_model2_response_results_WAIC, paste0("kgr_model2_response_results_WAIC.rds"))

