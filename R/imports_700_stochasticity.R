

#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: stochasticity of imports
#' author: "Martina Reichmuth"
#' date: "'26/07/2022"
#' ---
library(MASS)
library(lubridate)
library(MCMCglmm)
library(RDS)
library(doParallel)
library(foreach)
library(reshape2)
library(dplyr)

# Set seed
set.seed(60321)
registerDoParallel(cores=4)
getDoParWorkers()
n_sim <- 1e4

# Generation time
gen_mean <- 5.2 #Ganyani et al 2020; https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257
gen_sd <- 1.72 #Ganyani et al 2020; https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257
gen_var <- gen_sd^2
gamma_rate <- gen_mean/gen_var #gamma_shape/generation_time
gamma_shape <- gen_mean^2/gen_var
generation <- rgamma(n_sim, shape = gamma_shape, rate = gamma_rate)
generation <- ifelse(generation > 0, generation, rnorm(1, gen_mean, gen_sd))

# dispersion
dispersion_parameters <- as.vector(rtnorm(1e5,mean=0.51,lower=0.49,upper=0.52))
dispersion_parameters <- dispersion_parameters[order(dispersion_parameters)]

Rv_low <- 1.05
Rv_high <- 1.15
Rv_all <- as.vector(runif(n_sim,min(Rv_low),max(Rv_high)))

import_scenarios <- c(1,10,100)

stoch_imports <- function(Re, dispersion, import){
  cases_time <- c()
  foreach(Ri=1:length(Re)) %dopar% {
    R <- Re[Ri]
    cases_time <-secondary_t <- rep(1,import)
    k <- sample(dispersion, size=1)
    secondary <- 1
    while(length(cases_time)<1e5 & sum(secondary>0)>0) {#sum(secondary_t<max_time+0.5) >0 
      secondary <- rnbinom(length(secondary_t), size = k, mu = R)
      if(sum(secondary>0)>0){
      secondary_t <- rep(secondary_t, secondary)
      secondary_t <- secondary_t + round(rgamma(length(secondary_t), shape = gamma_shape, rate = gamma_rate))
      cases_time <- c(cases_time,secondary_t)
      }
    }
    model_outputs<- c(import,R,k,length(cases_time),cases_time)
    return(model_outputs)
  }
}

model_outputs <- c()
model_outputs <- lapply(1:length(import_scenarios), function(x) stoch_imports(Rv_all, dispersion_parameters, import_scenarios[x]))
model_outputs <- unlist(model_outputs, recursive = FALSE)
models_output_inc<- lapply(1:length(model_outputs), function(x) tabulate(as.numeric(model_outputs[[x]][-c(1:4)])))
models_output_cuminc <- lapply(1:length(models_output_inc), function(x) cumsum(models_output_inc[[x]]))


stochastic_model_outputs_variation <- c()
stochastic_model_outputs_variation <- data.frame(array(NA, dim = c(0, 6)))

for(i in 1:3){
  output_cuminc <- models_output_cuminc[((i-1)*n_sim+1):(i*n_sim)]
  output_inc <- models_output_inc[((i-1)*n_sim+1):(i*n_sim)]
  cum_inc <- sort(unique(unlist(output_cuminc)))
  cum_inc <- cum_inc[cum_inc<=5*1e3]
  output_cuminc <- as.data.frame(do.call(cbind, lapply(output_cuminc, `length<-`, max(lengths(output_cuminc)))))
  output_inc <- as.data.frame(do.call(cbind, lapply(output_inc, `length<-`, max(lengths(output_cuminc)))))
 
  output_cuminc$date <- c(1:length(output_cuminc[,1]))
  output_cuminc<- melt(output_cuminc, id.vars=c("date"))
  
  output_inc$date <- c(1:length(output_inc[,1]))
  output_inc<- melt(output_inc, id.vars=c("date"))
  output_cuminc<- cbind(output_cuminc, output_inc)
  output_cuminc <- output_cuminc[!output_cuminc[,6] %in% 0,]
  
  output_cuminc <- output_cuminc[!is.na(output_cuminc$value),]
  output_cuminc <- output_cuminc[,1:3]
  
  variation_outputs_range <-sapply(cum_inc, function(j) sum(diff(sort(output_cuminc$date[output_cuminc$value %in% j]))))
  variation_outputs_sd <-sapply(cum_inc, function(j) sd(output_cuminc$date[output_cuminc$value %in% j]))
  variation_outputs_mean <-sapply(cum_inc, function(j) mean(output_cuminc$date[output_cuminc$value %in% j]))
   variation_outputs_num <-sapply(cum_inc, function(j) sum(output_cuminc$value %in% j))
  
  variation_outputs <- as.data.frame(cbind(variation_outputs_range,variation_outputs_mean, variation_outputs_sd,cum_inc,variation_outputs_num))
  variation_outputs$seeds <- import_scenarios[i]
  colnames(stochastic_model_outputs_variation) <- colnames(variation_outputs) <- c("range","mean", "sd", "value","num_sim", "seeds")
  
  stochastic_model_outputs_variation <- rbind(variation_outputs,stochastic_model_outputs_variation)
}

days_max_plot<- 150
days_min_plot <- 20
for (i in 1:3) {
  output_inc <- models_output_inc[((i-1)*n_sim+1):(i*n_sim)]
  output_inc <- as.data.frame(do.call(rbind, lapply(output_inc, `length<-`, max(lengths(output_inc)))))
  output_inc$date <- c(1:length(output_inc[,1]))
  max_time<- length(output_inc)
  models_output_CI <- data.frame(array(as.numeric(0), dim = c(max_time, 6)))
  for(n in 1:max_time){
    models_output_CI[n,2:6] <- quantile(as.numeric(na.omit(output_inc[,n])), probs = c(.025,.25,.5,.75,.975))
    models_output_CI[n,7] <- length(na.omit(output_inc[,n]))
    
    }
  models_output_CI[,1] <- c(1:max_time)
  colnames(models_output_CI) <- c("date","low_95","low_50","median","up_50","up_95", "num_sim")
  if(i==1){col_l <- col_9[2]}
  if(i==2){col_l <- col_9[3]}
  if(i==3){col_l <- col_9[4]}
  max_yaxis <-400
  incidence_seeds_plot <- ggplot()+
    theme_minimal()+
    scale_y_continuous(labels = yscaling)+
    geom_ribbon(data=models_output_CI,aes(x=date, ymin=low_95,ymax=up_95), fill=col_l, alpha=0.4)+
    geom_ribbon(data=models_output_CI,aes(x=date, ymin=low_50,ymax=up_50), fill=col_l, alpha=0.4)+
    geom_line(data=models_output_CI, aes( x=(date), y=median), col=col_l)+#[models_output$variable%in%c(20001:20011),]
    theme(legend.position="none")+
    geom_hline(yintercept=3500, linetype="dashed", color = col_9[5])+
    coord_cartesian(ylim = c(1, max_yaxis),xlim = c(days_min_plot,days_max_plot)) +
    labs( x = "Time (in days)", y =bquote("Incidence"))
  if(i==1){incidence_seeds1_plot <- incidence_seeds_plot+ labs(tab="1 import")}
  if(i==2){incidence_seeds10_plot <- incidence_seeds_plot+ labs(tab="10 imports")}
  if(i==3){incidence_seeds100_plot <- incidence_seeds_plot+ labs(tab="100 imports")}
}
incidence_seeds_plot <- ggarrange(incidence_seeds1_plot, incidence_seeds10_plot,incidence_seeds100_plot,
                                  ncol = 3, nrow = 1,labels = c("","",""))

for (i in 1:3) {
  output_inc <- models_output_cuminc[((i-1)*n_sim+1):(i*n_sim)]
  output_inc <- as.data.frame(do.call(rbind, lapply(output_inc, `length<-`, max(lengths(output_inc)))))
  output_inc$date <- c(1:length(output_inc[,1]))
  max_time<- length(output_inc)
  models_output_CI <- data.frame(array(as.numeric(0), dim = c(max_time, 6)))
  for(n in 1:max_time){
    models_output_CI[n,2:6] <- quantile(as.numeric(na.omit(output_inc[,n])), probs = c(.025,.25,.5,.75,.975))
  }
  models_output_CI[,1] <- c(1:max_time)
  colnames(models_output_CI) <- c("date","low_95","low_50","median","up_50","up_95")
  if(i==1){col_l <- col_9[2]}
  if(i==2){col_l <- col_9[3]}
  if(i==3){col_l <- col_9[4]}
  max_yaxis <-5000
  cumincidence_seeds_plot <- ggplot()+
    theme_minimal()+
    scale_y_continuous(labels = yscaling)+
    geom_ribbon(data=models_output_CI,aes(x=date, ymin=low_95,ymax=up_95), fill=col_l, alpha=0.4)+
    geom_ribbon(data=models_output_CI,aes(x=date, ymin=low_50,ymax=up_50), fill=col_l, alpha=0.4)+
    geom_line(data=models_output_CI, aes( x=(date), y=median), col=col_l)+#[models_output$variable%in%c(20001:20011),]
    theme(legend.position="none")+
    coord_cartesian(ylim = c(1, max_yaxis),xlim = c(days_min_plot,days_max_plot)) +
    labs( x = "Time (in days)", y =bquote("Cumulative incidence"))
  if(i==1){incidence_seeds1_plot <- cumincidence_seeds_plot+ labs(tab="1 import")}
  if(i==2){incidence_seeds10_plot <- cumincidence_seeds_plot+ labs(tab="10 imports")}
  if(i==3){incidence_seeds100_plot <- cumincidence_seeds_plot+ labs(tab="100 imports")}
}
cumincidence_seeds_plot <- ggarrange(incidence_seeds1_plot, incidence_seeds10_plot,incidence_seeds100_plot,
                                     ncol = 3, nrow = 1,labels = c("","",""))

stochastic_model_outputs_variation$seeds <-  factor(stochastic_model_outputs_variation$seeds, levels=c("1", "10", "100"))
variation_cumincidence_sd <- ggplot()+
  theme_minimal()+
  geom_point(data=stochastic_model_outputs_variation, aes( x=(value), y=(sd), color=as.character(seeds), group=as.character(seeds)),alpha=0.3)+
  geom_vline(xintercept=3500, linetype="dashed", color = col_9[5])+
  coord_cartesian(ylim = c(days_min_plot,days_max_plot)) +
  scale_color_manual(name="Number of imports",
                     values = (col_9[2:4]),
                     labels = c("1", "10", "100"))+
  labs( x = "Cumulative incidence", y =bquote("Variation in days (sd)"))

variation_cumincidence_range <- ggplot()+
  theme_minimal()+
  geom_point(data=stochastic_model_outputs_variation, aes( x=(value), y=(range), color=as.character(seeds), group=as.character(seeds)),alpha=0.3)+
  geom_vline(xintercept=3500, linetype="dashed", color = col_9[5])+
  scale_color_manual(name="Number of imports",
                     values = (col_9[2:4]),
                     labels = c("1", "10", "100"))+
  labs( x = "Cumulative incidence", y =bquote("Variation in days (range)"))


stochastic_seeds_plot <- ggarrange(incidence_seeds_plot, cumincidence_seeds_plot,variation_cumincidence_sd,
                                     ncol = 1, nrow = 3,labels = c("A","B","C"))
ggsave(stochastic_seeds_plot, filename = paste0("./data/figures/SF6.png"), height =14, width = 10,  bg = "transparent")



