library(sandwich)
library(lmtest)
library(gamlr)
library(dplyr)


data <- all_full_report_stations_94_07_hourly
data <- subset(data, window365 != 0)
data$month2 <- NULL
data$month <- NULL

data$hp_opl1 <- data$opl1 * data$high_poverty
data$hd_opl1 <- data$opl1 * data$high_drate
data$cl_opl1 <- data$opl1 * data$station_near
data$rh_opl1 <- data$opl1 * data$hour_rush
data$time_opl1 <- data$opl1 * data$time
data$time2_opl1 <- data$opl1 * data$time_2
data$time3_opl1 <- data$opl1 * data$time_3
data$dow_hour <- data$dow * data$Hour

data_sub <- subset(data, window365 == 1)

gr_vars <- names(data_sub)[grep("^gr", names(data_sub))]
Wind_vars <- names(data_sub)[grep("^Wind", names(data_sub))]
Temp_vars <- names(data_sub)[grep("^Temp", names(data_sub))]
hum_vars <- names(data_sub)[grep("^hum", names(data_sub))]
mon_vars <- names(data_sub)[grep("^mon", names(data_sub))]

run_regression <- function(dependent_var, independent_vars, cluster_var, data) {
  formula_string <- paste(dependent_var, "~", paste(independent_vars, collapse = " + "))
  formula <- as.formula(formula_string)
  model <- lm(formula, data = data)
  cluster_se <- vcovCL(model, cluster = data[[cluster_var]])
  result <- coeftest(model, vcov = cluster_se)
  return(list(model = model, result = result))
}

independent_vars_poverty <- c(
  "high_poverty", "hp_opl1", gr_vars,
  "as.factor(opl1):time", "as.factor(opl1):time_2", "as.factor(opl1):time_3",
  Wind_vars, Temp_vars, hum_vars, mon_vars,
  "as.factor(dow):as.factor(Hour)", "opl1", "as.factor(dow)", "as.factor(Hour)"
)

regression1_lNOx <- run_regression("lNOx", independent_vars_poverty, "weekid", data_sub)
regression1_lCO <- run_regression("lCO", independent_vars_poverty, "weekid", data_sub)

print(regression1_lCO$result)
length(coef(regression1_lCO$result))-1

og_CO_estimate <- regression1_lCO$result["hp_opl1", "Estimate"] 
og_NOx_estimate <- regression1_lNOx$result["hp_opl1", "Estimate"]


##Placebo cutoffs
##THIS TAKES A LONG TIME TO RUN##
set.seed(12345678)
B <- 200
placebo_hour <- rep(0, B)
placebo_betas_CO <- rep(0, B)
placebo_betas_NOx <- rep(0, B)
for (b in 1:B)
{
  placebo <- sample.int(nrow(data_sub), size = 1, replace=FALSE)
  placebo_hour[b] <- placebo
  data_sub$oplextra <- ifelse(1:nrow(data_sub) > placebo, 1, 0)
  
  data_sub$hp_oplextra <- data_sub$oplextra * data_sub$high_poverty
  
  independent_vars_placebo <- c(
    "high_poverty", "hp_oplextra", gr_vars,
    "as.factor(oplextra):time", "as.factor(oplextra):time_2", "as.factor(oplextra):time_3",
    Wind_vars, Temp_vars, hum_vars, mon_vars,
    "as.factor(dow):as.factor(Hour)", "oplextra", "as.factor(dow)", "as.factor(Hour)"
  )
  
  placebo_regression_lCO <- run_regression("lCO", independent_vars_placebo, "weekid", data_sub)
  placebo_regression_lNOx <- run_regression("lNOx", independent_vars_placebo, "weekid", data_sub)
  
  placebo_betas_CO[b]<- placebo_regression_lCO$result["hp_oplextra", "Estimate"]
  placebo_betas_NOx[b]<- placebo_regression_lNOx$result["hp_oplextra", "Estimate"]
}

og_CO_estimate > quantile(abs(placebo_betas_CO),.95)
og_NOx_estimate > quantile(abs(placebo_betas_NOx),.95)

#No of bins based on Sturge’s Rule

placebo_betas_CO <- abs(placebo_betas_CO)
placebo_betas_NOx<- abs(placebo_betas_NOx)

hist(placebo_betas_CO, main = "CO estimates (n=200)", 
     xlab = "Estimates", ylab = "Density", col = "grey",border = "black",freq=FALSE,breaks=19)
lines(density(placebo_betas_CO), col = "red", lwd = 2)
legend("topright", legend = c("95th percentile","Original Estimate"), col = c("blue","purple"), lwd = 2)
abline(v=quantile(abs(placebo_betas_CO),.95), col="blue")
abline(v=og_CO_estimate, col="purple")

hist(placebo_betas_NOx, main = "NOx estimates (n=200)", 
     xlab = "Estimates", ylab = "Density", col = "grey",border = "black",freq=FALSE,breaks=19,ylim = c(0, 10))
lines(density(placebo_betas_NOx), col = "red", lwd = 2)
legend("topright", legend = c("95th percentile","Original Estimate"), col = c("blue","purple"), lwd = 2)
abline(v=quantile(abs(placebo_betas_NOx),.95), col="blue")
abline(v=og_NOx_estimate, col="purple")

## Our results show no significant heterogeneous effects between high poverty and low poverty areas


## Lasso
## Check to see what variables are good predictors of pollutant. 
## Gamlr is auto stantardized
##Don't want to drop an indicator variable
##Lambda is picked using AICc cause we dont want to cross validate with time series or RD. 

vars <- c(
  "high_poverty", "hp_opl1", gr_vars, "time_opl1","time2_opl1","time3_opl1","dow_hour",
  Wind_vars, Temp_vars, hum_vars, mon_vars,"opl1", "dow", "Hour","lNOx","lCO"
  )


data_sub_lasso0 <- subset(data_sub, is.finite(lNOx))
data_sub_lasso1 <- data_sub_lasso0[,vars]
data_sub_lasso <- data_sub_lasso1[complete.cases(data_sub_lasso1), ]

data_clean <- data_sub_lasso %>% select(-"lNOx", -"lCO")

data_clean$opl1 <- factor(data_clean$opl1)
data_clean$Hour <- factor(data_clean$Hour)
data_clean$dow <- factor(data_clean$dow)

data_clean_matrix <- model.matrix(~ . - 1, data = data_clean)
lasso_model <- gamlr(data_clean_matrix, data_sub_lasso$lNOx , verb=FALSE)
best_coefs <- coef(lasso_model, which = "AICc")
print(best_coefs)

lasso_model <- gamlr(data_clean_matrix, data_sub_lasso$lCO , verb=FALSE)
best_coefs <- coef(lasso_model, which = "AICc")
print(best_coefs)
