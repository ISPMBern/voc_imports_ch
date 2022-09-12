

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
registerDoParallel(cores=10)
getDoParWorkers()
# Set seed
set.seed(60321)
z = as.numeric(commandArgs(trailingOnly=TRUE))
#imports_007_stochasticity
n_sim <- 1e5

# Re (range? start voc and previous dif!)
#R_e = (1 - p)*R_w + p*R_v
#R_v = (1 + kappa)*R_w
#R_w = R_e/(1 + p*kappa)

# Generation time
gen_mean <- 5.2 #Ganyani et al 2020; https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257
gen_sd <- 1.72 #Ganyani et al 2020; https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257
gen_var <- gen_sd^2
gamma_rate <- gen_mean/gen_var #gamma_shape/generation_time
gamma_shape <- gen_mean^2/gen_var
generation <- rgamma(n_sim, shape = gamma_shape, rate = gamma_rate)
generation <- ifelse(generation > 0, generation, rnorm(1, gen_mean, gen_sd))

# dispersion
dispersion_parameters <- as.vector(rtnorm(1e5,mean=0.51,lower=0.49,upper=0.52))#on 2021-09-18#on 2021-05-20
dispersion_parameters <- dispersion_parameters[order(dispersion_parameters)]

max_time <- 160

#Swiss epidemic (FOPH)
swiss_cov <- read.csv("swiss_cov.csv")

#imports
# Alpha imports from Emma Hodcroft
imports_alpha <-read.csv("liberalClusters-Alpha.csv")[,-1]
imports_alpha$MinDateSwissChildren[grepl("2020-01-05", imports_alpha$MinDateSwissChildren)] <- ("2020-12-26")
imports_alpha$approach <- ifelse(imports_alpha$HeadOfClade!="True" & is.na(imports_alpha$Clade)| imports_alpha$HeadOfClade=="True" & !is.na(imports_alpha$Clade)  & imports_alpha$Clade>-1,"conservative", "liberal" )

# Delta imports from Emma Hodcroft
imports_delta <-read.csv("liberalClusters-Delta.csv")[,-1]
imports_delta <- imports_delta[!grepl("2021-02-13", imports_delta$MinDateSwissChildren),]
imports_delta$approach <- ifelse(imports_delta$HeadOfClade!="True" & is.na(imports_delta$Clade)| imports_delta$HeadOfClade=="True" & !is.na(imports_delta$Clade)  & imports_delta$Clade>-1,"conservative", "liberal" )

#imports_files <- list(imports_alpha, imports_delta)

#for (z in 1:6) {
  if(z %in% c(1,2,3)){
    imports_all <- imports_alpha
    zi <- z
    variant<- "Alpha"
    beta <- 0.43
  }
  if(z %in% c(4,5,6)){
    zi <- z-3
    imports_all <- imports_delta
    variant<- "Delta"
    beta <- 0.53
    
  }
  
  Rv <- (1 + beta)* na.omit(swiss_cov$median_R_mean[as_date(swiss_cov$date) %in% seq(as_date(min(imports_all$MinDateSwissChildren))-3,as_date(min(imports_all$MinDateSwissChildren))+3,by=1)])
  Rv_low <- (1 + beta)* na.omit(swiss_cov$median_R_lowHPD[as_date(swiss_cov$date) %in% seq(as_date(min(imports_all$MinDateSwissChildren))-3,as_date(min(imports_all$MinDateSwissChildren))+3,by=1)])
  Rv_high <- (1 + beta)* na.omit(swiss_cov$median_R_highHPD[as_date(swiss_cov$date) %in% seq(as_date(min(imports_all$MinDateSwissChildren))-3,as_date(min(imports_all$MinDateSwissChildren))+3,by=1)])
  
Rv_all <- as.vector(runif(n_sim,min(Rv_low),max(Rv_high)))

imports_afterfirst <- as.numeric(as_date(imports_all$MinDateSwissChildren)) - as.numeric(as_date(min(imports_all$MinDateSwissChildren)))
imports_afterfirst <- tabulate(findInterval(imports_afterfirst, vec=seq(0,length.out=max(imports_afterfirst))), nbins=max(imports_afterfirst))

first_import_date <- as_date(min(imports_all$MinDateSwissChildren))-7
last_import_date <- max(as_date(imports_all$MinDateSwissChildren))-7
length_imports <- length(first_import_date:last_import_date)


import_t <- rep(c(1:length(imports_afterfirst)),imports_afterfirst)
cases_d_runs <- data.frame(array(0, dim = c(max_time,1)))
import_scenarios <- list(import_t[1],import_t[1:10],import_t[1:100])


model_outputs <- data.frame(array(as.numeric(0), dim = c(n_sim, 7+max_time)))
colnames(model_outputs) <- c("seeds","imports","Re","dispersion_parameter","-","cum_cases","final_incidence", 1:max_time)



stoch_imports <- function(Re, dispersion, imports){
  #Ri <- Re #<- Rv
  #for(Ri in (1:length(Re))){
    #print(i)
  #foreach(i=1:length(import_t)) %dopar% {#foreach(Ri=1:length(Re)) %dopar% {
  foreach(Ri=1:length(Re)) %dopar% {#lapply(Re,function(R) {
  
    R <- Re[Ri]
    secondary_t <- imports#import_t[1:i]#import_scenarios[[zi]]#
    #R <- Re[Ri]
    k <- sample(dispersion, size=1)
    cases_d_runs[,1] <- tabulate(findInterval(imports, vec=seq(1,length.out=max(imports_afterfirst))), nbins=max_time)
    
    while(length(secondary_t<max_time+0.5) >0 & sum(cases_d_runs[,1])<1e7) {
      secondary <- rnbinom(length(secondary_t), size = k, mu = R)
      secondary_t <- rep(secondary_t, secondary)
      
      secondary_t <- secondary_t + round(rgamma(length(secondary_t), shape = gamma_shape, rate = gamma_rate))
      secondary_t <- secondary_t[secondary_t<max_time+0.5]
      cases_d_runs[match(names(table(secondary_t[secondary_t>0&secondary_t<(max_time+0.5)])),rownames(cases_d_runs)),1] <-  cases_d_runs[match(names(table(secondary_t[secondary_t>0&secondary_t<(max_time+0.5)])),rownames(cases_d_runs)),1] + table(secondary_t[secondary_t>0&secondary_t<(max_time+0.5)])
    }
    model_outputs[Ri,]<- c(length(import_scenarios[[zi]]),length(import_scenarios[[zi]]),R,k,NA,sum(cases_d_runs),cases_d_runs[max_time,1], cases_d_runs[,1])
    #rm(cases_d_runs)
    return(model_outputs[Ri,])
    }#)
}

model_outputs <- stoch_imports(Rv_all, dispersion_parameters, import_scenarios[[zi]])
model_outputs <- as.data.frame(do.call(rbind, lapply(model_outputs, `length<-`, max(lengths(model_outputs)))))

model_outputs$variant <- variant

# account for variation due to stochasticity
colnames(model_outputs)[is.na(colnames(model_outputs))] <- c(1:max_time)
models_output_num <- as.data.frame(t(model_outputs[,c(8:(length(model_outputs)-1))]))
models_output_num$date <- c(1:max_time)
models_output_num<- melt(models_output_num, id.vars=c("date"))
models_output_num$variant <- variant
models_output_num$seeds <- unique(model_outputs[,1])

models_output_num$time <- models_output_num$date
models_output_num$date <-  models_output_num$date+as_date(first_import_date)

models_output_num[,c("variable","cum_inc")]<- models_output_num %>% group_by(variable) %>% summarise(cum_inc=cumsum(value))
max_cumincidence_perseed<- models_output_num %>% group_by(seeds) %>% summarise(max_cumincidence_perseed = max(cum_inc))
variance_outputs <- data.frame(array(NA, dim = c(sum(max_cumincidence_perseed$max_cumincidence_perseed)+3, 5)))
colnames(variance_outputs) <- c("value","variance","sd","range", "seeds")

#for(k in 1:3){
  i <- c(1,10,100)[zi]
  subset_variance_outputs <- models_output_num[models_output_num$seeds==i,]
  #l <- sum(max_cumincidence_perseed$max_cumincidence_perseed[0:(k -1)])
  for(j in sort(as.numeric(unique(subset_variance_outputs$cum_inc)))){
    #s<- sort(as.numeric(unique(subset_variance_outputs$cum_inc)))[j]
    time_j<- subset_variance_outputs[subset_variance_outputs$cum_inc%in%j,]
    time_j <- time_j$time[!duplicated(time_j$variable)]
    time_j<- time_j - min(time_j)
    variance_outputs$variance[j] <- var(time_j)
    variance_outputs$sd[j] <- sd(time_j)
    variance_outputs$range[j] <- max(time_j)
    variance_outputs$value[j] <- j
    variance_outputs$seeds[j] <- i
  }
#}

variance_outputs<- variance_outputs[variance_outputs$seeds!=0,]
variance_outputs<- variance_outputs[!is.na(variance_outputs$seeds),]


if(z %in% c(1,2,3)){
saveRDS(model_outputs, paste0("import_voc_alpha_105_",zi,"_",Sys.Date(),".rds"))
write.csv(variance_outputs, paste0("import_variance_voc_alpha_105_",zi,"_",Sys.Date(),".csv"))
}
if(z %in% c(4,5,6)){
  saveRDS(model_outputs, paste0("import_voc_delta_105_",zi,"_",Sys.Date(),".rds"))
  write.csv(variance_outputs, paste0("import_variance_voc_delta_105_",zi,"_",Sys.Date(),".csv"))
}

#}


