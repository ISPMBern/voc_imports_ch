#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: text for manuscript
#' author: "Martina Reichmuth"
#' date: "'29/04/2022"
#' ---

for (n in 1:2) {
  if(n %in% c(1)){
    variant_data <- alpha_ch
    variant <- "Alpha"
    imports <- imports_alpha
    sim_voc <- sim_alpha
    scenario_closure <- scenarios_closure[scenarios_closure$variant=="Alpha",]
    scenario_closure_afterconcern <- scenarios_closure_afterconcern[scenarios_closure_afterconcern$variant=="Alpha",]
    scenarios_surveillance_variant <- scenarios_surveillance[scenarios_surveillance$variant=="Alpha",]

  }
  if(n %in% c(2)){
    variant_data <- delta_ch
    variant <- "Delta"
    imports <- imports_delta
    sim_voc <- sim_delta
    scenario_closure <- scenarios_closure[scenarios_closure$variant=="Delta",]
    scenario_closure_afterconcern <- scenarios_closure_afterconcern[scenarios_closure_afterconcern$variant=="Delta",]
    scenarios_surveillance_variant <- scenarios_surveillance[scenarios_surveillance$variant=="Delta",]
    
  }
  print(variant)
  cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(sim_voc$date),max(sim_voc$date),1))
  cov_ch_all <- subset(swiss_cov, as_date(date) %in% seq(min(sim_alpha$date),max(sim_delta$date),1))
  seq_ch_variant <- subset(seq_ch, as_date(date) %in% seq(min(sim_voc$date),max(sim_voc$date),1))
  
  print(paste0("reported cases:  ", round(sum(cov_ch$cases_num))))
  print(paste0("Average test positivity: ",   100*mean(swiss_cov$tests_pos_num[swiss_cov$date%in% c(min(alpha_ch$date):max(delta_ch$date))]/swiss_cov$tests_num[swiss_cov$date%in% c(min(alpha_ch$date):max(delta_ch$date))])))
  
  print(paste0("First sequence/import:  ", min(imports$MinDateSwissChildren)))
  print(paste0("proportion sequenced:  ", mean(swiss_cov$prop_seq_cases[swiss_cov$date%in% c(min(alpha_ch$date):max(delta_ch$date))])))
  print(paste0("variant sequenced:  ", round(sum(variant_data$variant))))
  print(paste0("variant proportion:  ", round(sum(variant_data$variant)/sum(seq_ch_variant$value)*100,2)))
 
  print(paste0("number of imports (liberal approach): ", length(imports[,1])))
  print(paste0("number of imports (conservative approach): ", length(imports[imports$approach=="conservative",1])))
  print(paste0("dates of imports of liberal approach: from ",min(imports$MinDateSwissChildren), " to ",max(imports$MaxDateSwissChildren)))
  print(paste0("dates of imports of conservative approach: from ",min(imports$MinDateSwissChildren[imports$approach=="conservative"]), " to ",max(imports$MaxDateSwissChildren[imports$approach=="conservative"])))
  
  print(paste0("simulated reported cases (liberal approach): ",sum(sim_voc$T_total[sim_voc$approach=="liberal"])))
  print(paste0("simulated reported cases (conservative approach): ",sum(sim_voc$T_total[sim_voc$approach=="conservative"])))
  
  print(paste0("simulated reported VoC cases (liberal approach): ",sum(sim_voc$T2[sim_voc$approach=="liberal"])))
  print(paste0("simulated reported VoC cases (conservative approach): ",sum(sim_voc$T2[sim_voc$approach=="conservative"])))
  
  print(paste0("simulated proportion VoC cases (liberal approach): ",round(sum(sim_voc$T2[sim_voc$approach=="liberal"])/sum(sim_voc$T_total[sim_voc$approach=="liberal"])*100,2)))
  print(paste0("simulated proportion VoC cases (conservative approach): ",round(sum(sim_voc$T2[sim_voc$approach=="conservative"])/sum(sim_voc$T_total[sim_voc$approach=="conservative"])*100,2)))
  
  print(paste0("sequenced cases:  ", round(sum(cov_ch$num_seq))))
  print(paste0("Variant sequenced (%):  ", round(sum(seq_ch_variant$who_variants %in% variant)/sum(seq_ch_variant$value)*100,2)))
  
  print(paste0("50% VoC (dynamic transmission model): ",min(sim_voc$date[sim_voc$approach=="liberal"&sim_voc$proportion>=0.5])))
  print(paste0("50% VoC (dynamic transmission model): ",min(sim_voc$date[sim_voc$approach=="conservative"&sim_voc$proportion>=0.5])))
  
  print(paste0("50% VoC:  ", min(variant_data$date[variant_data$proportion >=0.5])))
  print(paste0("Model lacks (liberal approach): ",min(variant_data$date[variant_data$proportion >=0.5])-min(sim_voc$date[sim_voc$approach=="liberal"&sim_voc$proportion>=0.5])))
  print(paste0("Model lacks (conservative approach): ",min(variant_data$date[variant_data$proportion >=0.5])-min(sim_voc$date[sim_voc$approach=="conservative"&sim_voc$proportion>=0.5])))
  
  print(paste0("Delay increases from 30 to 60 complete border closure by: ",scenarios_closure$delay_50[scenarios_closure$closed_borders%in%60]-scenarios_closure$delay_50[scenarios_closure$closed_borders%in%30]))
  print(paste0("Delay increases from 30 to 60 complete border closure after warning by: ",scenarios_closure_afterconcern$delay_50[scenarios_closure_afterconcern$closed_borders%in%60]-scenarios_closure_afterconcern$delay_50[scenarios_closure_afterconcern$closed_borders%in%30]))
  print(paste0("Delay increases from 50 to 100 complete border closure by: ",scenarios_closure$delay_50[scenarios_closure$closed_borders%in%100]-scenarios_closure$delay_50[scenarios_closure$closed_borders%in%50]))
  
  #print(paste0("Re: ", round(mean(cov_ch_all$median_R_mean),2)))
  #print(paste0("Re: ", quantile(cov_ch_all$median_R_mean)))
  
  print(paste0("50% Delta (prediction model): ",paste(min(na.omit(delta_ch$date[delta_ch$prediction>=0.5])), "(",min(na.omit(delta_ch$date[delta_ch$prediction.upper>=0.5])), " to ", min(na.omit(delta_ch$date[delta_ch$prediction.lower>=0.5])),")")))
  print(paste0("50% Alpha (prediction model): ",paste(min(na.omit(alpha_ch$date[alpha_ch$prediction>=0.5])), "(",min(na.omit(alpha_ch$date[alpha_ch$prediction.upper>=0.5])), " to ", min(na.omit(alpha_ch$date[alpha_ch$prediction.lower>=0.5])),")")))
  
  
  print(paste0("counterfactual scenario 1 month border closed results in a delay of: ",scenarios_closure$delay_50[scenarios_closure$closed_borders==30]))
  print(paste0("counterfactual scenario 1 month border closed after alert results in a delay of: ",scenarios_closure_afterconcern$delay_50[scenarios_closure_afterconcern$closed_borders==30]))
  print(paste0("counterfactual scenario 2 weeks border closed results in a delay of: ",scenarios_closure$delay_50[scenarios_closure$closed_borders==14]))
  print(paste0("counterfactual scenario 2 weeks border closed after alert results in a delay of: ",scenarios_closure_afterconcern$delay_50[scenarios_closure_afterconcern$closed_borders==14]))
  
  for(imports in import_scenarios){
    print(paste0((imports)," imports"))
    cum_inc <- 3500
    range_st  <- stochastic_model_outputs_variation$range[stochastic_model_outputs_variation$value %in% cum_inc & stochastic_model_outputs_variation$seeds %in% imports]
    sd_st <- stochastic_model_outputs_variation$sd[stochastic_model_outputs_variation$value %in% cum_inc & stochastic_model_outputs_variation$seeds %in% imports]
    print(paste0("Variation in time  ",round(range_st), " (sd: ",round(sd_st),")"))
  }
  
  
  print(paste0("counterfactual scenario improving surveillance by 25% results in a delay of: ",paste0(quantile(scenarios_surveillance_variant$delay_50[scenarios_surveillance_variant$surveillance >0.24 &scenarios_surveillance_variant$surveillance <0.26],probs = c(.025,.5,.975)))))
  print(paste0("counterfactual scenario improving surveillance by 50% results in a delay of: ",paste0(quantile(scenarios_surveillance_variant$delay_50[scenarios_surveillance_variant$surveillance >0.49 &scenarios_surveillance_variant$surveillance <0.51],probs = c(.025,.5,.975)))))
  print(paste0("counterfactual scenario improving surveillance by 75% results in a delay of: ",paste0(quantile(scenarios_surveillance_variant$delay_50[scenarios_surveillance_variant$surveillance >0.74 &scenarios_surveillance_variant$surveillance <0.76],probs = c(.025,.5,.975)))))
  
  samplex <- (scenarios_surveillance_variant$delay_50[scenarios_surveillance_variant$surveillance >0.74 &scenarios_surveillance_variant$surveillance <0.76])
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

