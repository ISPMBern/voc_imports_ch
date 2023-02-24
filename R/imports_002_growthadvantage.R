#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: estimate growth advantage of new variants
#' author: "Martina Reichmuth"
#' date: "'01/12/2021"
#' ---

parameters_output <-c()
for (n in 1:2) {
  if(n==1){
    # only data where first time >75% of variant
    variant_data <- alpha_ch[alpha_ch$date %in% seq(min(alpha_ch$date),min(alpha_ch$date[alpha_ch$proportion>0.75]),1),]
    variant_name <- "Alpha"
    time_window_Re <- as_date(c("2020-11-01", "2021-01-31")) # early growth phase of the variant: November 2020 to 31 January 2021; range: (0.8-0.98)
  }
  if(n==2){
    # only data where first time >75% of variant
    variant_data <- delta_ch[delta_ch$date %in% seq(min(delta_ch$date),min(delta_ch$date[delta_ch$proportion>0.75]),1),]
    variant_name <- "Delta" 
    time_window_Re <- as_date(c("2021-04-01", "2021-06-30")) # early growth phase of the variant: April 2021 to 30 June 2021; range: (0.62-1.64)
    
  }
  
  cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
  time_window <- c(min(variant_data$date),max(variant_data$date))
  time_window1 <- c(min(variant_data1$date),max(variant_data1$date))
  ## Sample size for parameter sets and bootstrapping
  n_sim <- 1e4
  # Generation time
  gen_mean <- 5.2 #Ganyani et al 2020; https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257
  gen_sd <- 1.72 #Ganyani et al 2020; https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257
  gen_var <- gen_sd^2
  gamma_rate <- gen_mean/gen_var #gamma_shape/generation_time
  gamma_shape <- gen_mean^2/gen_var
  generation <- rgamma(n_sim, shape = gamma_shape, rate = gamma_rate)
  generation <- ifelse(generation > 0, generation, rnorm(1, gen_mean, gen_sd))
  
  # effective reproduction number
  R_ch <- cov_ch$median_R_mean[cov_ch$date %in% seq(time_window_Re[1],time_window_Re[2],1)]
  paste0("range: (", min(R_ch),"-", max(R_ch), ")")
  R_ch <- sample(R_ch, n_sim, replace = TRUE)
  
  ## Estimate logistic growth rate of variants
  # Logistic growth function
  logistic <- function(p_0, rho, x) {
    p <- 1/((1 - p_0)/p_0*exp(- rho*x) + 1)
    return(p)
  }
  # Negative log-likelihood function
  nll <- function(p_0, rho) {
    p_0 <- plogis(p_0)
    rho <- rho
    times <- as.numeric(variant_data$date - min(variant_data$date))
    p <- logistic(p_0, rho, times)
    ll <- sum(dbinom(variant_data$variant, variant_data$total, p, log = TRUE))
    return(-ll)
  }
  # Fit logistic growth curve (using parameter transformation)
  free <- c(p_0 = qlogis(0.01), rho = 0.1)
  fit <- mle2(nll, start = as.list(free))#, method = "Nelder-Mead")
  fit.coef <- coef(fit)
  rho_ch <- fit.coef[2]
  fit.ci <- confint(fit)
  fit.varcov <- vcov(fit)
  sim_coef <- data.frame(rmvnorm(n_sim, mean = fit.coef, sigma = fit.varcov))
  rho_sample_ch <- sim_coef[, 2]
  sim_coef$p_0 <- plogis(sim_coef$p_0)
  
  logistic_out <- c()
  logistic_out$logistic_fit <- logistic(mean(sim_coef$p_0),rho_ch, as.numeric(seq(time_window[1],time_window[2],1))-as.numeric(time_window[1]))
  
  logistic_out$date <- seq(time_window[1],time_window[2],1)
  variant_data$logistic_fit <-logistic_out$logistic_fit[logistic_out$date %in% variant_data$date]
  
  p<- ggplot(data= variant_data)+
    geom_errorbar(aes(x = date, ymin=lower, ymax=upper), color= col_9[9], alpha=0.5,width=.1) +
    geom_point(aes(x = date, y=proportion), color= col_9[9])+
    geom_line(aes(x=date, y=logistic_fit), size=2)+
    theme_minimal()+
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0,1))+
    theme(plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 10),
          title = element_text(size = 10))+
    scale_x_date(date_breaks = "1 month", 
                 date_labels = "%b", limits = c(time_window1[1]-1,time_window1[2]+1))+
    labs(tag="",subtitle = bquote(), x = "", y =paste0("Proportion of ",variant_name,"\n"))
  
  ## safe parameters
  output <- data.frame(matrix(NA, nrow = 3, ncol = 6))
  names(output) <- c("parameter", "region", "median", "lower", "upper", "variant")
  # Rho (logistic growth rate)
  output[1, ] <- c("rho", "CH", round(rho_ch, 3), round(quantile(rho_sample_ch, probs = c(0.025, 0.975)), 3), NA)
  # Kappa (increased transmissibility)
  output[2, ] <- c("kappa", "CH", round(quantile(rho_sample_ch*generation/R_ch, probs = c(0.5, 0.025, 0.975)), 2), NA)
  # Tau ( generation time)
  output[3, ] <- c("tau", "CH", round(quantile(rho_sample_ch/(1/generation - rho_sample_ch), probs = c(0.5, 0.025, 0.975)), 2), NA)
  
  output$variant <- variant_name
  parameters_output <- rbind(parameters_output, output)
  
  
  
  # Prediction
  fit <- glm(cbind(variant, total - variant) ~ date, data = variant_data, family = "binomial")
  fit.ci <- confint(fit)
  interval = 0.95
  interval <- qnorm(1 - (1 - interval)/2)
  prediction_times <- data.frame(date = as_date(time_window[1]:time_window[2]))
  prediction <- predict(fit, newdata = prediction_times, type = "link", se.fit = TRUE)
  prediction.lower <- plogis(prediction$fit - interval*prediction$se.fit)
  prediction.upper <- plogis(prediction$fit + interval*prediction$se.fit)
  prediction <- plogis(prediction$fit)
  variant_data[,c("prediction", "prediction.lower","prediction.upper")] <- c(prediction, prediction.lower,prediction.upper)
  
  if(n %in% c(1)){
    plot_alpha_logistic_fit <- p
    alpha_ch<- merge(variant_data[,c("date","prediction", "prediction.lower","prediction.upper")], variant_data1, by="date", all.y = TRUE)
  }
  if(n %in% c(2)){
    plot_delta_logistic_fit <- p
    delta_ch <- merge(variant_data[,c("date","prediction", "prediction.lower","prediction.upper")], variant_data1, by="date", all.y = TRUE)
  }
}
plot_logistic_fit <- ggarrange(plot_alpha_logistic_fit, plot_delta_logistic_fit,
                               ncol = 2, nrow = 1)
ggsave(plot_logistic_fit, filename = paste0("./data/figures/SF2_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 4, width = 8,  bg = "transparent")
