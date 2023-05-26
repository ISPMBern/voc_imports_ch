#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: model spread of new variants in Switzerland
#' author: "Martina Reichmuth"
#' date: "'01/12/2021"
#' ---


### Load data
#Rows that have no entry for 'Clade' nor 'HeadOfClade' (end of file) are 'separate' nodes like b) below - they don't have any Swiss sequences above them. These are very likely introductions.
#Rows that have an entry for 'Clade' and 'HeadOfClade' are the first node in a chain, like node 1554 in a) below. They don't have any Swiss sequences above them. These are also likely introductions.
#Rows that have an entry for 'Clade' but not 'HeadOfClade' are the other nodes in a chain, like the other nodes in a) below. They can either be interpreted as further introductions, or not, depending on how you want to count.
# Alpha imports from Emma Hodcroft
alpha_files <- list.files(paste0(path_name,"data/Alpha/"))
conservative_alpha <- alpha_files[grepl("conservClusters_dates", alpha_files)]
conservative_alpha <- conservative_alpha[order(as.numeric(gsub("[^0-9]+","\\1",conservative_alpha)))]
liberal_alpha <- alpha_files[grepl("liberalClusters_dates", alpha_files)]
liberal_alpha <- liberal_alpha[order(as.numeric(gsub("[^0-9]+","\\1",liberal_alpha)))]
# Delta imports from Emma Hodcroft
delta_files <- list.files(paste0(path_name,"data/Delta/"))
conservative_delta <- delta_files[grepl("conservClusters_dates", delta_files)]
conservative_delta <- conservative_delta[order(as.numeric(gsub("[^0-9]+","\\1",conservative_delta)))]
liberal_delta <- delta_files[grepl("liberalClusters_dates", delta_files)]
liberal_delta <- liberal_delta[order(as.numeric(gsub("[^0-9]+","\\1",liberal_delta)))]


subsampling_i <- length(liberal_alpha)
lagdetection_l <- 7#c(3:13)
proportion_plot_Re_t_lags <- c()
simulations <-c()
imports_files <- c()
imports_estimates<- matrix(ncol=4, nrow=subsampling_i)
for(l in 1:length(lagdetection_l)){
for(i in 1:subsampling_i){
  imports_conservative_alpha <- read.csv(paste0(path_name,"data/Alpha/",conservative_alpha[i]))[,-1]
  imports_alpha <- read.csv(paste0(path_name,"data/Alpha/",liberal_alpha[i]))[,-1]
  imports_alpha$approach <- ifelse(imports_alpha$HeadOfClade!="True" & is.na(imports_alpha$Clade)| imports_alpha$HeadOfClade=="True" & !is.na(imports_alpha$Clade)  & imports_alpha$Clade>-1,"conservative", "liberal" )
  imports_alpha$MinDateSwissChildren[grepl("2020-01-05", imports_alpha$MinDateSwissChildren)] <- ("2020-12-26")
  imports_conservative_alpha$MinDateSwissChildren[grepl("2020-01-05", imports_conservative_alpha$MinDateSwissChildren)] <- ("2020-12-26")
  imports_alpha$MinDateNonSwissChildren <- as_date(imports_alpha$MinDateNonSwissChildren)
  
  imports_conservative_delta <- read.csv(paste0(path_name,"data/Delta/",conservative_delta[i]))[,-1]
  imports_delta <- read.csv(paste0(path_name,"data/Delta/",liberal_delta[i]))[,-1]
  imports_delta$approach <- ifelse(imports_delta$HeadOfClade!="True" & is.na(imports_delta$Clade)| imports_delta$HeadOfClade=="True" & !is.na(imports_delta$Clade)  & imports_delta$Clade>-1,"conservative", "liberal" )
  imports_delta$MinDateNonSwissChildren <- as_date(imports_delta$MinDateNonSwissChildren)
  
imports_files[[i]] <- list(imports_conservative_alpha, imports_alpha,imports_conservative_delta, imports_delta)
  
for (n in 1:length(imports_files[[i]])) {
  imports <- imports_files[[i]][[n]]
  if(n %in% c(1,2)){
    variant_data <- alpha_ch
    variant <- "alpha"
  }
  if(n %in% c(3,4)){
    variant_data <- delta_ch
    variant <- "delta"
  }
  cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
  time_window <- c(min(variant_data$date),max(variant_data$date))
  period <- seq(min(variant_data$date),max(variant_data$date),1)
  times_length <- length(period)
  
  imports <- (as.numeric(as_date(imports$MinDateSwissChildren)) - as.numeric(min(variant_data$date)))
  imports_sum <- length(imports)
  imports <-  imports-lagdetection_l[l] # delay to getting tested
  imports <- hist(imports, 0:times_length, plot = FALSE)$counts
  
  sim <- imports_300_SEIR(imports,variant)
  
  sim$variant <- variant
  sim$lagdetection <- lagdetection_l[l]
  
  if(n %in% c(1)){
    sim$approach <- "conservative"
    imports_estimates[is.na(imports_estimates[,1]),1][1] <- imports_sum
  }
   if(n %in% c(2)){
    sim$approach <- "liberal"
    imports_estimates[is.na(imports_estimates[,2]),2][1] <- imports_sum
   }
  if(n %in% c(3)){
    sim$approach <- "conservative"
    imports_estimates[is.na(imports_estimates[,3]),3][1] <- imports_sum
  }
  if(n %in% c(4)){
    sim$approach <- "liberal"
    imports_estimates[is.na(imports_estimates[,4]),4][1] <- imports_sum
  }
  simulations[[length(simulations)+1]] <- sim
  sim <- NULL
}

}
}
simulations_subsamples <- simulations

for (n in 1:length(imports_estimates[1,])) {
print(paste0( round(mean(imports_estimates[,n])), " (range: ",round(quantile(imports_estimates[,n],c(0))),"-",round(quantile(imports_estimates[,n],c(1))),")")) #round(mean(imports_estimates[,n])-qnorm(0.975)*sd(imports_estimates[,n])/sqrt(subsampling_i)),"-",round(mean(imports_estimates[,n])+qnorm(0.975)*sd(imports_estimates[,n])/sqrt(subsampling_i)),")"))
}

for(l in length(lagdetection_l)) {
for(m in 1:2) {
  if(m %in% c(1)){
    labels_variant <- c("All cases","Alpha cases", "other cases")
    titel_variant <- bquote("Proportion of Alpha\n")
    variant_name <- "Alpha cases" 
    variant_data <- alpha_ch
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
    time_window <- c(min(variant_data$date),max(variant_data$date))
    period <- seq(min(variant_data$date),max(variant_data$date),1)
    variant <- "alpha"
  }
  if(m %in% c(2)){
    labels_variant <- c("All cases","Delta cases", "other cases")
    titel_variant <- bquote("Proportion of Delta\n")
    variant_name <- "Delta cases" 
    variant_data <- delta_ch
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
    time_window <- c(min(variant_data$date),max(variant_data$date))
    period <- seq(min(variant_data$date),max(variant_data$date),1)
    variant <- "delta"
  }
  
  sim <- as.data.frame(do.call(rbind, lapply(simulations_subsamples, `length<-`, max(lengths(simulations_subsamples)))))
  sim <- sim[sim$variant==variant,]
  sim$Imports[is.na(sim$Imports)] <- 0
  
  sim_adapt <- c()
  
  for(approach in c("conservative","liberal")){
    s <- sim[sim$approach==approach,]
    sims <- c()
    sims$time <- unique(sim$time)
    sims$date <- unique(sim$date)
    sims <- as.data.frame(sims) 
    sims[,c("T2","T2_lower","T2_upper","T_total","T_total_lower","T_total_upper", "Imports","Imports_lower","Imports_upper")] <- NA
    for(d in sims$time){
      
      sims[d,c("T2","T2_lower","T2_upper")]<- quantile(s$T2[s$time==d], c(0.5,0,1))#c(0.5,0.025,0.975))
      
      sims[d,c("T_total","T_total_lower","T_total_upper")] <- quantile(s$T_total[s$time==d], c(0.5,0,1))
      
      sims[d,c("Imports","Imports_lower","Imports_upper")] <- quantile(s$Imports[s$time==d],c(0.5,0,1))
      
     }
    sims$variant <- variant
    sims$approach <- approach

    sim_adapt <- rbind(sims,sim_adapt)
  }
  
  sim_adapt$Imports_upper[sim_adapt$Imports==0] <- NA
  sim_adapt$Imports_lower[sim_adapt$Imports==0] <- NA
  sim_adapt$Imports[sim_adapt$Imports==0] <- NA
  sim_adapt$approach <- factor(sim_adapt$approach, levels =  c("liberal","conservative"))
  
  time_window_imports <- max(sim$date[sim$Imports>0])
  
  sim <- as.data.frame(do.call(rbind, lapply(simulations_baseline, `length<-`, max(lengths(simulations_baseline)))))
  sim <- sim[sim$variant==variant,]
  sim$Imports[is.na(sim$Imports)] <- 0
  sim$Imports[sim$Imports==0] <- NA
  sim$approach <- factor(sim$approach, levels =  c("liberal","conservative"))
  
  
  # Plot
  imports_hist_plot <- ggplot()+
   geom_point(sim_adapt, mapping=aes(x=date, y=Imports, fill=approach,color=approach),alpha=0.5)+
    geom_errorbar(sim_adapt, mapping=aes(x=date, ymin=Imports_lower,ymax=Imports_upper,color=approach))+
    scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0), limits = c(time_window[1]-1,time_window_imports+1))+
    coord_cartesian(xlim=c(time_window[1]-1,time_window_imports+1))+
    scale_y_continuous(limits=c(0,100))+
    theme_minimal()+
    scale_fill_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    scale_color_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    theme(legend.position=c(.2,.7),
          plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 10),
          title = element_text(size = 10))+
    labs(tag="",subtitle = bquote(), x = "", y ="Number of imports \n")


  proportion_plot_Re_t <- ggplot(data= sim_adapt)+
    geom_errorbar(data= variant_data, aes(x = date, ymin=lower, ymax=upper), color= col_9[9], alpha=0.5,width=.1) +
    geom_point(data= variant_data,aes(x = date, y=proportion), color= col_9[9])+
    geom_ribbon(data= sim_adapt, aes(x = date, y = T2/T_total, ymin = T2_lower/T_total_lower,ymax = T2_upper/T_total_upper, fill=approach), alpha=0.8)+
    geom_line(data= sim, aes(x = date, y = T2/T_total, color=approach))+
    scale_color_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    scale_fill_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    theme_minimal()+
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0,1))+
    theme(plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 10),
          title = element_text(size = 10))+
    scale_x_date(date_breaks = "1 month", 
                 date_labels = "%b", limits = c(time_window[1]-1,time_window[2]+1))+
    labs(tag="",subtitle = bquote(), x = "", y =paste0(titel_variant,"\n"))

  
  cases_day_plot_Re_t <-  ggplot(data= sim_adapt)+
    geom_line(data = cov_ch, aes(x = date, y = weigthed_cases),color=col_9[9], size=1)+
    geom_ribbon(data= sim_adapt, aes(x = date, y = T2, ymin = T2_lower,ymax = T2_upper, fill=approach))+
    scale_color_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    scale_fill_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    theme_minimal()+
    scale_y_continuous(limits = c(0,8100))+
    scale_x_date(date_breaks = "1 month", 
                 date_labels = "%b", limits = c(time_window[1]-1,time_window[2]+1))+
    theme_minimal()+
    theme(legend.position = "none",
          plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 10),
          title = element_text(size = 10))+
    labs(tag="",subtitle = "", x = "", y =bquote("Number of reported \ncases per day"))
  
  if(m %in% c(1)){
    plot_model_ch_alpha <- ggarrange(imports_hist_plot, cases_day_plot_Re_t,proportion_plot_Re_t, 
                                     labels = c("A", "C", "E"),
                                     ncol = 1, nrow = 3)

 }
  if(m %in% c(2)){
    plot_model_ch_delta <- ggarrange(imports_hist_plot, cases_day_plot_Re_t,proportion_plot_Re_t, 
                                     labels = c("B", "D", "F"),
                                     ncol = 1, nrow = 3)
    
  }
  
}
plot_model_ch <- ggarrange(plot_model_ch_alpha, plot_model_ch_delta,
                           ncol = 2, nrow = 1)
  ggsave(plot_model_ch, filename = paste0("./data/figures/SF3.pdf"), height = 7, width = 7,  bg = "transparent")
  ggsave(plot_model_ch, filename = paste0("./data/figures/SF3.png"), height = 7, width = 7,  bg = "transparent")
}
  
  