#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: text for manuscript
#' author: "Martina Reichmuth"
#' date: "'29/04/2022"
#' ---

for (n in 1:2) {
  if(n %in% c(1)){
    variant_data <- alpha_ch
    variant <- "alpha"
    imports <- imports_alpha
    sim_voc <- sim_alpha
    scenario_closure <- scenarios_closure[scenarios_closure$variant=="Alpha",]
    scenario_closure_afterconcern <- scenarios_closure_afterconcern[scenarios_closure_afterconcern$variant=="Alpha",]
    scenarios_surveillance <- scenarios_surveillance[scenarios_surveillance$variant=="Alpha",]

  }
  if(n %in% c(2)){
    variant_data <- delta_ch
    variant <- "delta"
    imports <- imports_delta
    sim_voc <- sim_delta
    scenario_closure <- scenarios_closure[scenarios_closure$variant=="Delta",]
    scenario_closure_afterconcern <- scenarios_closure_afterconcern[scenarios_closure_afterconcern$variant=="Delta",]
    scenarios_surveillance <- scenarios_surveillance[scenarios_surveillance$variant=="Delta",]
    
  }
  print(variant)
  
  cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
  cov_ch1 <- subset(swiss_cov, as_date(date) %in% seq(min(alpha_ch$date),max(delta_ch$date),1))
  
  cov_seq_ch <- subset(variant_data, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
  print(paste0("cases reported: ", sum(cov_ch$cases_num)))
  print(paste0("Average test positivity: ",   100*mean(swiss_cov$tests_pos_num[swiss_cov$date%in% c(min(alpha_ch$date):max(delta_ch$date))]/swiss_cov$tests_num[swiss_cov$date%in% c(min(alpha_ch$date):max(delta_ch$date))])))
  print(paste0("Average test positivity: ",   100*quantile(swiss_cov$tests_pos_num[swiss_cov$date%in% c(min(alpha_ch$date):max(delta_ch$date))]/swiss_cov$tests_num[swiss_cov$date%in% c(min(alpha_ch$date):max(delta_ch$date))])))
  
  print(paste0("50% VoC:  ", min(cov_seq_ch$date[cov_seq_ch$proportion>=0.5])))
  
  print(paste0("proportion sequenced:  ", round(sum(cov_seq_ch$variant)/sum(cov_seq_ch$total)*100,2)))
  print(paste0("proportion sequenced:  ", mean(swiss_cov$prop_seq_cases[swiss_cov$date%in% c(min(alpha_ch$date):max(delta_ch$date))])))
  print(paste0("variant proportion:  ", quantile(round((cov_seq_ch$variant)/(cov_seq_ch$total)*100,2))))
  print(paste0("Re: ", round(mean(cov_ch1$median_R_mean),2)))
  print(paste0("Re: ", quantile(cov_ch1$median_R_mean)))
  
  print(paste0("number of imports: ", length(imports[,1])))
  print(paste0("number of imports: ", length(imports[imports$approach=="conservative",1])))
  
  print(paste0("dates of imports of liberal approach: from ",min(imports$MinDateSwissChildren), " to ",max(imports$MaxDateSwissChildren)))
  print(paste0("dates of imports of conservative approach: from ",min(imports$MinDateSwissChildren[imports$approach=="conservative"]), " to ",max(imports$MaxDateSwissChildren[imports$approach=="conservative"])))
  
  print(paste0("simulated reported cases: ",sum(sim_voc$T_total[sim_voc$approach=="liberal"])))
  print(paste0("simulated reported cases: ",sum(sim_voc$T_total[sim_voc$approach=="conservative"])))
  
  print(paste0("simulated reported VoC cases: ",sum(sim_voc$T2[sim_voc$approach=="liberal"])))
  print(paste0("simulated reported VoC cases: ",sum(sim_voc$T2[sim_voc$approach=="conservative"])))
  
  print(paste0("simulated proportion VoC cases: ",round(sum(sim_voc$T2[sim_voc$approach=="liberal"])/sum(sim_voc$T_total[sim_voc$approach=="liberal"])*100,2)))
  print(paste0("simulated proportion VoC cases: ",round(sum(sim_voc$T2[sim_voc$approach=="conservative"])/sum(sim_voc$T_total[sim_voc$approach=="conservative"])*100,2)))
  
  
  print(paste0("50% Alpha (prediction model): ",paste(min(delta_ch$date[delta_ch$prediction>=0.5]), "(",min(delta_ch$date[delta_ch$prediction.upper>=0.5]), " to ", min(delta_ch$date[delta_ch$prediction.lower>=0.5]),")")))
  print(paste0("50% Delta (prediction model): ",paste(min(alpha_ch$date[alpha_ch$prediction>=0.5]), "(",min(alpha_ch$date[alpha_ch$prediction.upper>=0.5]), " to ", min(alpha_ch$date[alpha_ch$prediction.lower>=0.5]),")")))
  
  print(paste0("50% VoC (dynamic transmission model): ",min(sim_voc$date[sim_voc$approach=="liberal"&sim_voc$proportion>=0.5])))
  print(paste0("50% VoC (dynamic transmission model): ",min(sim_voc$date[sim_voc$approach=="conservative"&sim_voc$proportion>=0.5])))
  
  
  print(paste0("counterfactual scenario 1 month border closed results in a delay of: ",scenarios_closure$delay_50[scenarios_closure$closed_borders==30]))
  print(paste0("counterfactual scenario 1 month border closed after alert results in a delay of: ",scenarios_closure_afterconcern$delay_50[scenarios_closure_afterconcern$closed_borders==30]))
  print(paste0("counterfactual scenario 2 weeks border closed results in a delay of: ",scenarios_closure$delay_50[scenarios_closure$closed_borders==14]))
  print(paste0("counterfactual scenario 2 weeks border closed after alert results in a delay of: ",scenarios_closure_afterconcern$delay_50[scenarios_closure_afterconcern$closed_borders==14]))
  
  print(paste0("counterfactual scenario improving surveillance by 25% results in a delay of: ",paste0(quantile(scenarios_surveillance$time_50[scenarios_surveillance$surveillance >0.24 &scenarios_surveillance$surveillance <0.26]))))
  print(paste0("counterfactual scenario improving surveillance by 25% results in a delay of: ",paste0(quantile(scenarios_surveillance$time_50[scenarios_surveillance$surveillance >0.49 &scenarios_surveillance$surveillance <0.51]))))
  print(paste0("counterfactual scenario improving surveillance by 75% results in a delay of: ",paste0(quantile(scenarios_surveillance$time_50[scenarios_surveillance$surveillance >0.74 &scenarios_surveillance$surveillance <0.76]))))
  
  samplex <- (scenarios_surveillance$delay_50[scenarios_surveillance$surveillance >0.74 &scenarios_surveillance$surveillance <0.76])
  samplex_boot <- rnorm(n_sim, samplex)
 plot(samplex)
 plot(samplex_boot)
  
  sd_surveillance_delay <- sd(samplex_boot)
  mean_surveillance_delay <- mean(samplex_boot)  
  n_surveillance_delay <- length(samplex_boot)
  mean_surveillance_delay
  mean_surveillance_delay+1.96 *(sd_surveillance_delay/sqrt(n_surveillance_delay))
  mean_surveillance_delay-1.96 *(sd_surveillance_delay/sqrt(n_surveillance_delay))

  
  
  time_window <- c(min(variant_data$date),max(variant_data$date))
  period <- seq(min(variant_data$date),max(variant_data$date),1)
  times_length <- length(period)
  
  
}

